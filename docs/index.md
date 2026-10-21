# ETNA: Epitaxial Tcad Numerical Analysis

ETNA is a Python TCAD library, with a Streamlit UI, for **4H-SiC p+/n-/n+ radiation
detectors** used for proton dosimetry and microdosimetry. It solves Poisson plus
electron/hole drift-diffusion in 1D and 2D with [devsim](https://devsim.org/) and two
numba engines that reproduce devsim inside checked parameter boxes. It is a research
code, not a commercial TCAD package and not clinical software.

This site renders the documents in the repository's `docs/` folder. The repository
[README](https://github.com/ngcex/sic-tcad-detector#readme) is the short overview.

## Where to start

| You want to                                                   | Read                                                                                                                       |
| ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| install ETNA and run a first C-V / CCE simulation             | [Usage guide, §3-4](USAGE.md#3-installation)                                                                               |
| know which outputs are verified, calibrated or only predicted | [Usage guide, §2 "What to trust"](USAGE.md#2-what-to-trust) and [§8 "Validity and limits"](USAGE.md#8-validity-and-limits) |
| look up a `run_*` facade and its arguments                    | [Usage guide, §6](USAGE.md#6-the-run_-facades)                                                                             |
| send measurement or foundry data to calibrate the model       | [Data request checklist](DATA_REQUEST_CHECKLIST.md)                                                                        |
| cut a release                                                 | [Releasing](RELEASING.md)                                                                                                  |

## Project

- Source, issues and releases: <https://github.com/ngcex/sic-tcad-detector>. ETNA is
  distributed from GitHub only; there is no PyPI release.
- Changes: [CHANGELOG.md](https://github.com/ngcex/sic-tcad-detector/blob/main/CHANGELOG.md).
- How to cite: [CITATION.cff](https://github.com/ngcex/sic-tcad-detector/blob/main/CITATION.cff)
  (GitHub's "Cite this repository" button).
- Contributing, support and security:
  [CONTRIBUTING.md](https://github.com/ngcex/sic-tcad-detector/blob/main/CONTRIBUTING.md),
  [SUPPORT.md](https://github.com/ngcex/sic-tcad-detector/blob/main/SUPPORT.md),
  [SECURITY.md](https://github.com/ngcex/sic-tcad-detector/blob/main/SECURITY.md).
- Licence: MIT, copyright Cesare Scalia. Third-party licences:
  [THIRD_PARTY_LICENSES.md](https://github.com/ngcex/sic-tcad-detector/blob/main/THIRD_PARTY_LICENSES.md).
