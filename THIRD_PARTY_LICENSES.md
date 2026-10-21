# Third-party licences and data terms

ETNA itself is MIT-licensed, copyright Cesare Scalia (see [LICENSE](LICENSE)). This file
inventories:

1. the licences of the Python packages ETNA needs at run time;
2. the native solver libraries those packages bundle;
3. the terms of the third-party data that ETNA redistributes in this repository and in
   the wheel (`etna/_data/`).

It is an engineering inventory, **not legal advice**. Licences are given as each project
declares them. Items marked **unverified** have not been checked against a primary
source, or need a human decision (issue #3, H4).

## 1. Python runtime dependencies

ETNA does not vendor or redistribute any of these packages. They are installed by the
user's installer (uv) from their own distributions, under their own licences. The table
covers the locked runtime closure (`uv export --frozen --no-dev --no-emit-project`; the
`dev` and `viz` extras are not included). Each licence is read from the installed package
metadata by `scripts/license_inventory.py`. Regenerate the table after a lock change:

```bash
uv sync --extra dev
uv run python scripts/license_inventory.py --write   # rewrites the table below
uv run python scripts/license_inventory.py --check   # stale table or unreviewed flag -> exit 1
```

Flags: **copyleft** marks a GPL-family declaration (GPL, LGPL, AGPL, MPL or EPL), and
**unknown** marks missing metadata. Every flagged package needs a reviewed note in the
script. At the last generation (2026-10-21, macOS arm64, Python 3.13):

- no GPL, LGPL or AGPL package is declared in the runtime closure;
- the only flag is **certifi** (MPL-2.0, a file-level weak copyleft covering the CA
  bundle, used unmodified);
- devsim needs the separate note in section 2.

<!-- BEGIN GENERATED: scripts/license_inventory.py -->
| Package | Locked version | Declared licence (package metadata) | Flag | Note |
| --- | --- | --- | --- | --- |
| altair | 6.2.2 | BSD License | - | - |
| anyio | 4.15.1 | MIT | - | - |
| attrs | 26.1.0 | MIT | - | - |
| blinker | 1.9.0 | MIT License | - | - |
| cachetools | 7.1.4 | MIT | - | - |
| certifi | 2026.6.17 | MPL-2.0 | copyleft | MPL-2.0 (CA bundle data, unmodified; weak copyleft, file-level). |
| charset-normalizer | 3.4.7 | MIT | - | - |
| click | 8.4.2 | BSD-3-Clause | - | - |
| colorama | 0.4.6 | BSD License | - | Only when `sys_platform == 'win32'`; licence as declared on PyPI. |
| contourpy | 1.3.3 | BSD License | - | - |
| cycler | 0.12.1 | BSD License | - | - |
| devsim | 2.10.0 | Apache License, Version 2.0 | - | Apache-2.0 for devsim itself; its wheel bundles libumfpack_lgpl (UMFPACK 5.1, LGPL-2.1+ per the wheel's NOTICE). See the devsim section below. |
| fonttools | 4.63.0 | MIT | - | - |
| gitdb | 4.0.12 | BSD License | - | - |
| gitpython | 3.2.1 | BSD-3-Clause | - | - |
| h11 | 0.16.0 | MIT | - | - |
| httptools | 0.8.0 | MIT | - | - |
| idna | 3.18 | BSD-3-Clause | - | - |
| itsdangerous | 2.2.0 | BSD License | - | - |
| jinja2 | 3.1.6 | BSD License | - | - |
| jsonschema | 4.26.0 | MIT | - | - |
| jsonschema-specifications | 2025.9.1 | MIT | - | - |
| kiwisolver | 1.5.0 | BSD License | - | - |
| llvmlite | 0.50.0 | BSD-2-Clause AND Apache-2.0 WITH LLVM-exception | - | - |
| markupsafe | 3.0.3 | BSD-3-Clause | - | - |
| matplotlib | 3.11.0 | Python Software Foundation License | - | - |
| narwhals | 2.23.0 | MIT | - | - |
| numba | 0.68.0 | BSD | - | - |
| numpy | 2.5.0 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 | - | - |
| packaging | 26.2 | Apache-2.0 OR BSD-2-Clause | - | - |
| pandas | 3.0.3 | BSD License | - | - |
| pillow | 12.3.0 | MIT-CMU | - | - |
| plotly | 6.8.0 | MIT | - | - |
| protobuf | 7.35.1 | 3-Clause BSD License | - | - |
| pyarrow | 24.0.0 | Apache-2.0 | - | - |
| pydeck | 0.9.2 | Apache License 2.0 | - | - |
| pyparsing | 3.3.2 | MIT | - | - |
| python-dateutil | 2.9.0.post0 | Dual License | - | - |
| python-multipart | 0.0.32 | Apache-2.0 | - | - |
| referencing | 0.37.0 | MIT | - | - |
| requests | 2.34.2 | Apache-2.0 | - | - |
| rpds-py | 2026.6.3 | MIT | - | - |
| scipy | 1.18.0 | BSD License | - | - |
| six | 1.17.0 | MIT | - | - |
| smmap | 5.0.3 | BSD-3-Clause | - | - |
| starlette | 1.3.1 | BSD-3-Clause | - | - |
| streamlit | 1.58.0 | Apache-2.0 | - | - |
| tenacity | 9.1.4 | Apache 2.0 | - | - |
| toml | 0.10.2 | MIT | - | - |
| typing-extensions | 4.16.0 | PSF-2.0 | - | - |
| tzdata | 2026.2 | Apache-2.0 | - | Only when `sys_platform == 'emscripten' or sys_platform == 'win32'`; licence as declared on PyPI. |
| urllib3 | 2.8.0 | MIT | - | - |
| uvicorn | 0.49.0 | BSD-3-Clause | - | - |
| watchdog | 6.0.0 | Apache-2.0 | - | Only when `sys_platform != 'darwin'`; licence as declared on PyPI. |
| websockets | 16.0 | BSD-3-Clause | - | - |
<!-- END GENERATED -->

**Optional extras** (not inventoried above): `viz` (gdstk, manim) and `dev` (pytest,
jupyter, ruff and so on). Users who install them should check their licences. manim
pulls in native libraries (cairo, pango) through the system package manager.

## 2. devsim and its bundled solver libraries (UMFPACK)

- **devsim** is declared **Apache-2.0** in its package metadata ("Apache License, Version
  2.0", classifier "Apache Software License"). Its NOTICE says: "Copyright 2013--2023
  DEVSIM LLC".
- **Bundled UMFPACK (LGPL).** The devsim 2.10.0 wheel installed here (macOS arm64)
  contains a separate shared library, `libumfpack_lgpl.dylib`. It ships at
  `devsim/libumfpack_lgpl.dylib`, plus an identical copy at
  `devsim/umfpack/libumfpack_lgpl.dylib`. devsim's Python solver callback
  (`devsim/umfpack/umfpack_loader.py`) loads the library at run time through `ctypes`,
  by the relative name `./libumfpack_lgpl.<so|dylib>`. ETNA's
  `etna/core/devsim_fastsolve.py` only caches that callback's setup. It does not ship
  or link UMFPACK.
- **What the NOTICE says.** The devsim wheel's `NOTICE` file
  (`devsim-2.10.0.dist-info/licenses/NOTICE`) lists "UMFPACK (Version 5.1)" and "AMD
  Version 2.2", both under the **GNU LGPL 2.1 or later**. It also lists:
  - COLAMD (permissive, University of Florida);
  - SuperLU and LAPACK (BSD-style);
  - Boost (BSL-1.0);
  - Eigen (MPL-2.0);
  - Symdiff (Apache-2.0).
- **Assessment** (not legal advice). ETNA only depends on devsim as a separately
  installed package, and devsim loads UMFPACK as a separate, replaceable shared library.
  On that reading, the LGPL obligations rest with the distributor of the devsim wheel,
  and ETNA's MIT licence is not affected.
- **Observed on macOS arm64** (devsim 2.10.0, 2026-10-21; loaded images listed through
  `_dyld_get_image_name` after a devsim solve, plus devsim's own verbose output):
  - devsim prints "loading UMFPACK 5.1 as direct solver";
  - the library loaded is the bundled copy at `devsim/umfpack/libumfpack_lgpl.dylib`;
  - BLAS: devsim searches `DEVSIM_MATH_LIBS="libopenblas.dylib:liblapack.dylib:libblas.dylib"`,
    reports `libopenblas.dylib` as missing, and loads the system `liblapack.dylib`. The
    images in use are Apple's Accelerate `libBLAS` / `libLAPACK` (part of macOS). No
    Homebrew OpenBLAS is needed.
- **Still UNVERIFIED (open question from issue #3):**
  - that the same UMFPACK build and NOTICE ship in the Linux and Windows devsim wheels;
    only the macOS arm64 wheel was inspected;
  - which BLAS is loaded on Linux. The loader asks for `libopenblas` first, and the CI
    installs `libopenblas-dev` (BSD-3-Clause) on Ubuntu;
  - whether devsim's UMFPACK 5.1 build matches the licence text quoted in the NOTICE.

  A redistributor of a bundle that contains devsim (a Docker image, a frozen app) should
  check these points before shipping.

## 3. Redistributed data

The repository and, where noted, the wheel (`etna/_data/`, listed in `pyproject.toml`
`force-include`) contain the data below. Provenance details are in each directory's
README.

| Data                                                                                                                               | Where                                                                        | In the wheel                                 | Source                                                                                                         | Terms                                                                                                                                                  |
| ---------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| NIST PSTAR / ASTAR stopping powers (Si, C, water, ICRU muscle): CSVs derived from raw responses                                    | `data/nist_star/`                                                            | CSVs and README; raw responses sdist only    | NIST STAR programs (physics.nist.gov), fetched as described in `data/nist_star/README.md`                      | **Unverified.** US-government NIST Standard Reference Data. NIST SRD terms of use apply. They were not opened for this inventory (H4)                  |
| NIST ESTAR electron tables (raw HTML, parsed at run time)                                                                          | `data/nist_estar/raw/`                                                       | yes (6 HTML files)                           | NIST ESTAR (physics.nist.gov), see `data/nist_estar/README.md`                                                 | **Unverified**, as for NIST STAR above                                                                                                                 |
| SR-NIEL proton NIEL in Si and SiC: CSVs and the scraped result pages                                                               | `data/niel/sr_niel_*`                                                        | CSVs only; `.html` pages sdist/repo only     | SR-NIEL web calculator (sr-niel.org), SR-NIEL-7 ver. 11.1, fetched 2026-10-20 by `data/niel/fetch_sr_niel.py`  | **Unverified.** The site's terms of use were not checked. Shipping only the derived CSVs with attribution, or dropping the raw HTML, may be safer (H4) |
| Neutron NIEL / damage energy for Si and SiC                                                                                        | `data/niel/neutron_*`                                                        | CSVs yes; `.csv.gz` pointwise data repo only | **Computed here** with NJOY2016 (HEATR) from ENDF/B-VIII.0 evaluations (NNDC), `scripts/build_neutron_niel.py` | ETNA-generated (MIT). Underlying ENDF/B evaluations: terms **unverified** (generally distributed openly by NNDC)                                       |
| Si neutron damage function (Griffin et al., SAND92-0094, as compiled by Vasilescu & Lindström for RD48)                            | `data/niel/si_neutron_damage_griffin_rd48.csv`                               | yes                                          | archived RD48 `neutrons.xls`, verbatim (see `data/niel/README.md`)                                             | **Unverified**                                                                                                                                         |
| Points digitised from published figures: Petringa et al. 2025, _Physica Medica_ 140, 105675, Fig. 4b; the "Photons" paper, Fig. 4b | `data/validation/`                                                           | yes                                          | digitised in this repository (`lab/microdosimeter/`)                                                           | **Unverified (H4).** Numbers read off figures are facts, not the figure itself, but publisher policy and co-author consent have not been checked       |
| Legacy SiC/water stopping-power CSVs (`stopping_power_{sic,water}.csv`) used by the opt-in `legacy_unsafe` kappa path              | `data/`                                                                      | yes                                          | provenance not documented. The code labels this path as a fabricated artefact of v3.0                          | **Unknown provenance.** Kept only for an explicit legacy opt-in. Candidate for removal                                                                 |
| CCE(LET) tables, layout defaults, measurement templates                                                                            | `data/cce_let_table_*.json`, `data/foundry/`, `data/measurements/templates/` | yes                                          | generated or written in this repository                                                                        | MIT (ETNA)                                                                                                                                             |
| Synthetic MC events                                                                                                                | `data/synthetic_mc_events.csv`                                               | no (sdist/repo only)                         | generated in this repository                                                                                   | MIT (ETNA)                                                                                                                                             |

Open items for a human (issue #3 H4):

- the NIST SRD terms;
- the SR-NIEL terms, and whether to keep the scraped HTML;
- publisher and co-author terms for the digitised points;
- the provenance, or removal, of the legacy stopping-power CSVs.
