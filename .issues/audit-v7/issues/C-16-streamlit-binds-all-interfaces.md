---
id: C-16
source: raw-C F16
title: Streamlit app listens on all interfaces by default with no auth, 200 MB uploads and unbounded CPU-bound solves
severity: low
category: security
status: open
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
