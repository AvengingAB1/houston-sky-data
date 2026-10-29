# houston-sky-data

Forecast feed for the **Sunsetglow** app's Houston region (development / demo / GPS hardware-in-the-loop testing).
Published with GitHub Pages: https://avengingab1.github.io/houston-sky-data/

| File | Source | Workflow |
|---|---|---|
| `hrrr_latest.json` | NOAA HRRR 3 km low/mid/high cloud + visibility, 27-33 N / 99.5-91.5 W, +48 h | `hrrr.yml` hourly probe |
| `elements_latest.json` | HRRR 10 m wind (+ arrows) and hourly precipitation overlays | `hrrr.yml` |
| `stations*.json` | MET Norway Locationforecast point forecasts (global) for the app's places, towns, grid | `hrrr.yml` |
| `ecmwf_latest.json`, `ecmwf_elements_latest.json` | ECMWF IFS 9 km via Open-Meteo open data (AWS), to +15 days | `meteo.yml` 08:30/20:30 UTC |
| `meteo_latest.json` | Open-Meteo best_match grid 1 x 1 deg (63 points, once a day) | `meteo.yml` |

HRRR is read from the NOAA Open Data Dissemination bucket on AWS (`noaa-hrrr-bdp-pds`) with HTTP range requests
for only the needed GRIB2 messages (`extract/extract_hrrr.py`). Credits: NOAA/NCEP HRRR; ECMWF (CC BY 4.0);
MET Norway (NLOD / CC BY 4.0); Open-Meteo (CC BY 4.0).
