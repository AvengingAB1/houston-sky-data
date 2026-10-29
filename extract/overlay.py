"""overlay.py - shared drawing for the Wind/Precip overlay JSON (elements_latest.json): fixed colour
bands (one legend serves every frame), wind arrows, Web-Mercator-spaced image rows so Leaflet's
imageOverlay registers exactly. Same bands and look as the Norway feeds (extract_meps_elements.py)."""
from __future__ import annotations
import base64
import io
import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ALPHA = 190
WIND_EDGES = [0, 3, 5, 8, 11, 14, 17, 21, 25, 29]
WIND_RGB = [(70, 110, 170), (60, 150, 190), (70, 180, 140), (120, 200, 90), (220, 220, 70),
            (240, 170, 50), (235, 110, 40), (215, 50, 50), (170, 30, 110), (120, 40, 160)]
PRECIP_EDGES = [0.1, 0.3, 0.5, 1, 2, 5, 10]
PRECIP_RGB = [(170, 210, 245), (120, 180, 240), (70, 140, 230), (40, 95, 205),
              (60, 50, 180), (140, 40, 170), (210, 40, 120)]


def _ym(lat):
    return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))


def image_height(box, width):
    """Pixel height that keeps the Mercator aspect of box for a given width."""
    return int(round(width * (_ym(box['north']) - _ym(box['south'])) / math.radians(box['east'] - box['west'])))


def pixel_latlon(box, w, h):
    yN, yS = _ym(box['north']), _ym(box['south'])
    rows = [(yN - (r + 0.5) / h * (yN - yS)) for r in range(h)]
    lats = np.degrees(np.arctan(np.sinh(np.array(rows))))
    lons = box['west'] + (np.arange(w) + 0.5) / w * (box['east'] - box['west'])
    return np.meshgrid(lats, lons, indexing='ij')


def colourise(v, edges, rgb, transparent_below=True):
    idx = np.digitize(v, edges) - 1
    out = np.zeros(v.shape + (4,), np.uint8)
    pal = np.array(rgb, np.uint8)
    ok = idx >= 0 if transparent_below else np.ones(v.shape, bool)
    out[..., :3] = pal[np.clip(idx, 0, len(rgb) - 1)]
    out[..., 3] = np.where(ok & ~np.isnan(v), ALPHA, 0)
    return out


def draw_arrows(img, spd, dirn, step=40):
    d = ImageDraw.Draw(img)
    h, w = spd.shape
    for r in range(step // 2, h, step):
        for c in range(step // 2, w, step):
            s, frm = spd[r, c], dirn[r, c]
            if np.isnan(s) or np.isnan(frm) or s < 1:
                continue
            to = math.radians(frm + 180.0)            # meteorological FROM -> arrow points TO
            L = 9 + min(s, 25) * 0.35
            dx, dy = math.sin(to) * L, -math.cos(to) * L
            x0, y0, x1, y1 = c - dx / 2, r - dy / 2, c + dx / 2, r + dy / 2
            d.line([(x0, y0), (x1, y1)], fill=(20, 20, 30, 230), width=2)
            for a in (150, -150):
                ha = to + math.radians(a)
                d.line([(x1, y1), (x1 + math.sin(ha) * 4, y1 - math.cos(ha) * 4)], fill=(20, 20, 30, 230), width=2)


def webp_b64(img):
    buf = io.BytesIO()
    img.save(buf, format='WEBP', lossless=True, method=4)
    return base64.b64encode(buf.getvalue()).decode('ascii')


def legend_b64(edges, rgb, fmt):
    cell, h = 34, 14
    img = Image.new('RGBA', (cell * len(rgb), h * 2), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    for i, col in enumerate(rgb):
        d.rectangle([i * cell, 0, (i + 1) * cell - 1, h - 1], fill=col + (255,))
        d.text((i * cell + 2, h + 1), fmt(edges[i]), fill=(25, 25, 25, 255), font=font)
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode('ascii')


def wind_png(spd, dirn):
    img = Image.fromarray(colourise(spd, WIND_EDGES, WIND_RGB, transparent_below=False), 'RGBA')
    draw_arrows(img, spd, dirn)
    return webp_b64(img)


def precip_png(mm):
    return webp_b64(Image.fromarray(colourise(mm, PRECIP_EDGES, PRECIP_RGB), 'RGBA'))


def legends(n):
    wl = legend_b64(WIND_EDGES, WIND_RGB, lambda v: '%d' % v)
    pl = legend_b64(PRECIP_EDGES, PRECIP_RGB, lambda v: '%g' % v)
    return [wl] + [None] * (n - 1), [pl] + [None] * (n - 1)
