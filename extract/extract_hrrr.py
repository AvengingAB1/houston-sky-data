#!/usr/bin/env python3
"""
extract_hrrr.py - NOAA HRRR 3 km (CONUS) for the Houston region, one pass, two outputs:
  hrrr_latest.json      low/mid/high cloud cover + surface visibility in the app's cloud-grid schema
                        (same as dmi_latest.json / meps_latest.json: row 0 = north, UTC times, run
                        'YYMMDD_HHMM', div 20 = 5 % steps, u8-b64 frames, vis in 50 m steps)
  elements_latest.json  10 m wind (+ arrows) and hourly precipitation overlays (overlay.py bands)

Source: NOAA HRRR on AWS Open Data (s3://noaa-hrrr-bdp-pds, public, no key). One GRIB2 file per
forecast hour for all of CONUS (~150 MB); only the 8 messages we need are fetched with HTTP Range
requests via the wgrib2 .idx (~8 MB per hour). The synoptic runs (00/06/12/18Z) reach +48 h; the
newest one whose +48 h file is published is used (about 2 h after its init time).
HRRR cloud layers: LCDC low (surface..642 hPa), MCDC middle (642..350 hPa), HCDC high (above 350 hPa).
10 m winds are grid-relative on the Lambert grid and are rotated to true north here.

Exit codes: 0 written, 3 newest complete run already published (--skip-if-run), 4 source unreachable.

Usage: python extract_hrrr.py [--run YYYY-MM-DDTHH] [--out hrrr_latest.json]
                              [--elements-out elements_latest.json] [--skip-if-run YYMMDD_HHMM] [--hours 48]
"""
from __future__ import annotations
import argparse
import base64
import datetime as dt
import gzip
import json
import math
import os
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import eccodes
from scipy.spatial import cKDTree

import overlay as O

_ECC = threading.Lock()                  # ecCodes is not thread-safe: downloads run in parallel, decoding one at a time

S3 = 'https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.%s/conus/hrrr.t%02dz.wrfsfcf%02d.grib2'
# (short name, level) -> our key. APCP is matched separately (the hourly bucket of this hour).
FIELDS = {('LCDC', 'low cloud layer'): 'low', ('MCDC', 'middle cloud layer'): 'mid',
          ('HCDC', 'high cloud layer'): 'high', ('VIS', 'surface'): 'vis',
          ('UGRD', '10 m above ground'): 'u', ('VGRD', '10 m above ground'): 'v'}
HOURS = 48
# Output grid: the Houston map box with ~300 km up-sun margin (sunset W/WNW, sunrise over the Gulf).
BOX = dict(lat0=27.0, lat1=33.0, lon0=-99.5, lon1=-91.5)
DLAT, DLON = 0.03, 0.034                                     # ~3.3 x 3.3 km at 30 N (HRRR 3 km)
IMG_W = 640                                                  # elements overlay width (px)
DIV, NODATA = 20, 255
VIS_STEP, VIS_CLEAR = 50, 101
UA = {'User-Agent': 'sunsetglow-data/1.0 github.com/avengingab1/houston-sky-data'}
CREDIT = 'NOAA HRRR 3 km (NOAA Open Data Dissemination, AWS)'


class SourceUnreachable(Exception):
    pass


def url(run, h):
    return S3 % (run.strftime('%Y%m%d'), run.hour, h)


def http(u, rng=None, tries=4, timeout=120):
    h = dict(UA)
    if rng:
        h['Range'] = 'bytes=%d-%s' % (rng[0], '' if rng[1] is None else rng[1])
    err = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                raise
            err = e
        except Exception as e:                               # timeouts, resets
            err = e
        time.sleep(3 * (k + 1))
    raise err


def exists(u):
    try:
        http(u, rng=(0, 15), tries=2, timeout=30)
        return True
    except urllib.error.HTTPError:
        return False


def newest_complete_run(hours=HOURS, now=None):
    """Newest synoptic run (00/06/12/18Z) whose +hours file and .idx are published."""
    now = now or dt.datetime.now(dt.timezone.utc)
    reached = False
    for back in range(0, 6):
        t = now - dt.timedelta(hours=6 * back)
        run = t.replace(hour=t.hour - t.hour % 6, minute=0, second=0, microsecond=0)
        try:
            ok = exists(url(run, hours) + '.idx') and exists(url(run, hours))
            reached = True
        except Exception:
            continue
        if ok:
            return run
    if not reached:
        raise SourceUnreachable('noaa-hrrr-bdp-pds not answering')
    return None


def decode(msg, latlon=False):
    with _ECC:
        return _decode(msg, latlon)


def _decode(msg, latlon=False):
    g = eccodes.codes_new_from_message(msg)
    try:
        vals = eccodes.codes_get_values(g).astype('float32')
        if eccodes.codes_get(g, 'bitmapPresent'):
            vals[vals == eccodes.codes_get(g, 'missingValue')] = np.nan
        meta = {}
        if latlon:
            meta['lat'] = eccodes.codes_get_array(g, 'latitudes')
            meta['lon'] = eccodes.codes_get_array(g, 'longitudes')
            meta['LoV'] = eccodes.codes_get(g, 'LoVInDegrees')
            meta['Latin1'] = eccodes.codes_get(g, 'Latin1InDegrees')
            meta['Latin2'] = eccodes.codes_get(g, 'Latin2InDegrees')
        try:
            meta['gridRel'] = bool(eccodes.codes_get(g, 'uvRelativeToGrid'))
        except Exception:
            meta['gridRel'] = None
        return vals, meta
    finally:
        eccodes.codes_release(g)


def fetch_hour(run, h, latlon=False):
    """{key: 1D field} for one forecast hour (+ grid metadata if asked)."""
    u = url(run, h)
    rows = [l.split(':') for l in http(u + '.idx').decode().strip().splitlines()]
    offs = [int(r[1]) for r in rows]
    want = {}
    for i, r in enumerate(rows):
        key = FIELDS.get((r[3], r[4]))
        if key is None and r[3] == 'APCP' and r[4] == 'surface' and r[5] == '%d-%d hour acc fcst' % (h - 1, h):
            key = 'apcp'
        if key and key not in want:
            want[key] = (offs[i], offs[i + 1] - 1 if i + 1 < len(offs) else None)
    out, meta = {}, {}
    for key, rng in want.items():
        vals, m = decode(http(u, rng=rng), latlon=latlon and 'lat' not in meta)
        out[key] = vals
        if key == 'u':
            meta['gridRel'] = m['gridRel']
        if 'lat' in m:
            meta.update(m)
    need = set(FIELDS.values()) | ({'apcp'} if h > 0 else set())
    missing = need - set(out)
    if missing:
        raise RuntimeError('hour +%d: fields %s not in %s' % (h, sorted(missing), u))
    return out, meta


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype=np.uint8).tobytes()).decode('ascii')


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def build(run, hours=HOURS):
    t0 = time.time()
    first, meta = fetch_hour(run, 0, latlon=True)
    lat, lon = meta['lat'], np.where(meta['lon'] > 180, meta['lon'] - 360, meta['lon'])
    near = (lat > BOX['lat0'] - 0.5) & (lat < BOX['lat1'] + 0.5) & (lon > BOX['lon0'] - 0.7) & (lon < BOX['lon1'] + 0.7)
    src = np.flatnonzero(near)
    tree = cKDTree(unit(lat[src], lon[src]))

    def sampler(qlat, qlon, max_km):
        dist, k = tree.query(unit(qlat, qlon).reshape(-1, 3))
        return src[k], (dist * 6371.0 > max_km)

    # cloud grid (regular lat/lon)
    out_lats = np.round(np.arange(BOX['lat1'], BOX['lat0'] - 1e-9, -DLAT), 5)   # N -> S
    out_lons = np.round(np.arange(BOX['lon0'], BOX['lon1'] + 1e-9, DLON), 5)    # W -> E
    glat, glon = np.meshgrid(out_lats, out_lons, indexing='ij')
    gpick, gfar = sampler(glat, glon, 2.5)
    gshape = glat.shape
    # overlay pixels (Mercator rows)
    ebox = dict(south=BOX['lat0'], north=BOX['lat1'], west=BOX['lon0'], east=BOX['lon1'])
    img_h = O.image_height(ebox, IMG_W)
    plat, plon = O.pixel_latlon(ebox, IMG_W, img_h)
    ppick, pfar = sampler(plat, plon, 3.0)
    pshape = plat.shape
    pfar = pfar.reshape(pshape)
    # grid-relative -> true-north wind rotation (Lambert conformal): angle = n * (lon - LoV)
    la1, la2 = math.radians(meta['Latin1']), math.radians(meta['Latin2'])
    n = math.sin(la1) if abs(la1 - la2) < 1e-6 else (
        math.log(math.cos(la1) / math.cos(la2)) / math.log(math.tan(math.pi / 4 + la2 / 2) / math.tan(math.pi / 4 + la1 / 2)))
    lov = meta['LoV'] - 360 if meta['LoV'] > 180 else meta['LoV']
    ang = np.radians(n * (((lon[ppick] - lov + 180) % 360) - 180)).reshape(pshape)
    cosA, sinA = np.cos(ang), np.sin(ang)

    def sample(fields, h):
        res = {}
        for k in ('low', 'mid', 'high', 'vis'):
            v = fields[k][gpick]
            v[gfar] = np.nan
            res[k] = v.reshape(gshape)
        u, v = fields['u'][ppick].reshape(pshape), fields['v'][ppick].reshape(pshape)
        if meta.get('gridRel', True) is not False:
            u, v = cosA * u + sinA * v, -sinA * u + cosA * v
        spd = np.hypot(u, v)
        dirn = (np.degrees(np.arctan2(-u, -v)) + 360) % 360              # direction the wind blows FROM
        spd[pfar] = np.nan
        res['spd'], res['dir'] = spd, dirn
        if 'apcp' in fields:
            p = fields['apcp'][ppick].reshape(pshape)
            p[pfar] = np.nan
            res['apcp'] = p
        return res

    got = {0: sample(first, 0)}
    del first, lat, lon
    with ThreadPoolExecutor(8) as ex:
        for h, s in zip(range(1, hours + 1), ex.map(lambda h: sample(fetch_hour(run, h)[0], h), range(1, hours + 1))):
            got[h] = s
    print('  fetched %d hours in %.0fs' % (len(got), time.time() - t0))

    def q_cloud(x):
        return np.where(np.isnan(x), NODATA, np.clip(np.rint(x / 100.0 * DIV), 0, DIV)).astype(np.uint8)

    def q_vis(x):
        return np.where(np.isnan(x), NODATA, np.minimum(np.rint(x / VIS_STEP), VIS_CLEAR)).astype(np.uint8)

    hrs = sorted(got)
    times = [(run + dt.timedelta(hours=h)).strftime('%Y-%m-%dT%H:%MZ') for h in hrs]
    rid = run.strftime('%y%m%d_%H%M')
    gen = dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    clouds = {
        'source': 'noaa_hrrr_3km', 'region': 'houston',
        'scale': 'fraction_div', 'div': DIV, 'nodata': NODATA,
        'run': rid, 'generated': gen,
        'lat0': float(out_lats[0]), 'lat1': float(out_lats[-1]),
        'lon0': float(out_lons[0]), 'lon1': float(out_lons[-1]),
        'nlat': len(out_lats), 'nlon': len(out_lons), 'dlat': DLAT, 'dlon': DLON,
        'frames': hrs, 'times': times, 'cloud_enc': 'u8-b64',
        'low': [b64(q_cloud(got[h]['low'])) for h in hrs],
        'mid': [b64(q_cloud(got[h]['mid'])) for h in hrs],
        'high': [b64(q_cloud(got[h]['high'])) for h in hrs],
        'vis': [b64(q_vis(got[h]['vis'])) for h in hrs],
        'vis_enc': 'u8-b64', 'vis_step': VIS_STEP, 'vis_clear': VIS_CLEAR,
        'credit': CREDIT,
    }
    wl, pl = O.legends(len(hrs))
    elements = {
        'source': 'noaa_hrrr_3km_elements', 'region': 'houston', 'run': rid, 'generated': gen,
        'north': BOX['lat1'], 'south': BOX['lat0'], 'west': BOX['lon0'], 'east': BOX['lon1'],
        'img_w': IMG_W, 'img_h': img_h, 'fmt': 'webp', 'mercator_rows': True,
        'frames': hrs, 'times': times,
        'wind': [O.wind_png(got[h]['spd'], got[h]['dir']) for h in hrs],
        'precip': [O.precip_png(got[h]['apcp']) if 'apcp' in got[h] else None for h in hrs],
        'wind_legends': wl, 'precip_legends': pl,
        'credit': CREDIT,
    }
    return clouds, elements


def write(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(obj, f, separators=(',', ':'))
    os.replace(tmp, path)
    raw = open(path, 'rb').read()
    return len(raw) / 1e6, len(gzip.compress(raw, 6)) / 1e6


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', help='YYYY-MM-DDTHH (default: newest complete synoptic run)')
    ap.add_argument('--out', default='hrrr_latest.json')
    ap.add_argument('--elements-out', default='elements_latest.json')
    ap.add_argument('--skip-if-run')
    ap.add_argument('--hours', type=int, default=HOURS)
    a = ap.parse_args(argv)
    try:
        run = (dt.datetime.strptime(a.run, '%Y-%m-%dT%H').replace(tzinfo=dt.timezone.utc)
               if a.run else newest_complete_run(a.hours))
    except SourceUnreachable as e:
        print('HRRR unreachable (%s); keeping the published run' % e)
        return 4
    if run is None:
        raise SystemExit('no complete HRRR synoptic run in the last 36 h')
    rid = run.strftime('%y%m%d_%H%M')
    if a.skip_if_run and rid == a.skip_if_run:
        print('run %s already published; nothing to do' % rid)
        return 3
    print('run=%s' % rid)
    clouds, elements = build(run, a.hours)
    mb, gz = write(a.out, clouds)
    print('wrote %s: run=%s %dx%d x %d frames, %.1f MB (gzip %.1f MB)' % (
        a.out, rid, clouds['nlat'], clouds['nlon'], len(clouds['times']), mb, gz))
    if a.elements_out:
        mb, gz = write(a.elements_out, elements)
        print('wrote %s: %dx%d px x %d frames, %.1f MB' % (a.elements_out, elements['img_w'], elements['img_h'], len(elements['times']), mb))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
