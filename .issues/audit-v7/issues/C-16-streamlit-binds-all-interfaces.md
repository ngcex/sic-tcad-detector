---
id: C-16
source: raw-C F16
title: Streamlit app listens on all interfaces by default with no auth, 200 MB uploads and unbounded CPU-bound solves
severity: low
category: security
status: fixed
location:
  - .streamlit/config.toml (only a [theme] section; no [server] section)
  - README.md:82 (`uv run streamlit run app/main.py`)
labels: [security, deployment, streamlit]
---

## Summary

With no `[server]` config, Streamlit binds `*` (all IPv4 interfaces), has no authentication,
and allows 200 MB uploads. Anyone on the same network can start multi-minute devsim solves.
Because of C-04, two overlapping solves segfault the server, so this is also a trivial remote
DoS. They can also upload 200 MB CSVs that `pd.read_csv` loads fully into RAM.

The raw file marked this "plausible" because it checked only config defaults. I verified the
bind live, so it is upgraded to open.

## Reproduction / Evidence

```
uv run python -c "from streamlit import config; print(config.get_option('server.address'), config.get_option('server.maxUploadSize'))"
-> None 200

uv run streamlit run app/main.py --server.headless true --server.port 8599 &   # stopped afterwards
lsof -nP -iTCP:8599 -sTCP:LISTEN
-> python3.1 36326 ngcex 6u IPv4 ... TCP *:8599 (LISTEN)
Streamlit log: Local URL: http://localhost:8599
               Network URL: http://192.168.1.230:8599
               External URL: http://<public-ip>:8599
```

`cat .streamlit/config.toml` shows a [theme] section only. The uploader's `type=["csv"]` is a
client-side filter.

## Impact

Low for the intended single-researcher laptop use. On shared or institutional networks it is an
unauthenticated compute endpoint that can be crashed (C-04).

## Proposed fix

```diff
 # .streamlit/config.toml
+[server]
+address = "localhost"
+maxUploadSize = 20
+headless = true
```

In README, add: "for shared use, run with `--server.address 0.0.0.0` behind an authenticating
reverse proxy".

## Test plan

Config only. Manually start the app and check that `lsof -nP -iTCP:8501 -sTCP:LISTEN` shows
`127.0.0.1:8501`. Optionally add a unit test that parses `.streamlit/config.toml` and asserts
`server.address == "localhost"`.

## Effort

S

## Status 2026-10-18 (final-c): fixed (c6c2e51)

- .streamlit/config.toml gains `[server] address = "localhost"` (theme kept). maxUploadSize is set explicitly to Streamlit default 200 MB on purpose: a per-step MC CSV is ~70 B/step (data/synthetic_mc_events.csv 776 kB / 11k steps), so 20 MB would block real Geant4 exports, and with the loopback bind uploads come only from the local user. headless is not set in the file (passed on the CLI for health checks; locally the browser should open). For shared use: `--server.address 0.0.0.0` behind an authenticating proxy (comment in the file).
- Evidence: `uv run streamlit run app/main.py --server.headless true --server.port 8597` -> `curl -s localhost:8597/_stcore/health` = ok; `lsof -nP -iTCP:8597 -sTCP:LISTEN` -> `TCP 127.0.0.1:8597 (LISTEN)` (was `*:8599` at audit time); server killed, port free. tests/test_app_streamlit_config.py asserts server.address == "localhost".
- Scope (checked): Streamlit reads `.streamlit/config.toml` from the launch directory (plus ~/.streamlit). Launched from /tmp (`uv run --project <repo> --directory /tmp streamlit run <repo>/app/main.py --server.headless true --server.port 8598`), the server listened on `*:8598`. The committed config therefore makes localhost the default only when launched from the repo root (the documented command). From any other directory, keep `--server.address localhost`.
