---
id: D-15
title: run_microdosimetry / Streamlit page default sv_width_um=150 vs DESIGN-2's only validated SV width (100 um)
severity: low
category: docs
status: closed 2026-10-18 (no numerical mismatch: mean_chord_length uses 4V/S only with width AND depth, so run_microdosimetry's sv_width_um is inert and l_bar = 2 * thickness = 20 um for 150 and 100 um alike; 91194a5 documents it in the facade docstring and the app help, test_sv_width_does_not_enter_l_bar pins it. Default kept (design decision). The finite-SV l_bar question is physics, see D-INDEX)
location: deliverables/DESIGN-2-microdosimetry-geometry.md:35,97-98; etna/api/simulation.py:843-847 (sv_width_um=150.0, used at ~:888 mean_chord_length); app/workflows/microdosimetry.py:57-58
labels: [docs, microdosimetry, defaults, foundry-deliverable]
source: raw-D F15
---

## Summary

DESIGN-2 gives the SV half-width as 50 µm (a 100 µm SV) and says widths other than 50 µm half-width "are
not quoted" until uniformity sweeps are rerun. The public facade `run_microdosimetry` and the Streamlit
page both default to `sv_width_um=150.0` (full width, passed to `mean_chord_length`). The out-of-the-box
spectrum therefore uses a geometry that the deliverable does not validate.

Nuance that lowers the impact: DESIGN-2 says outright that it "makes no claim about lineal-energy spectra"
(lines 95-96). In the facade, the width only enters the mean chord length l̄, not a CCE-uniformity claim.
The mismatch is real but narrow, so it stays at low.

## Reproduction/Evidence

```
$ sed -n 33,36p deliverables/DESIGN-2-microdosimetry-geometry.md   -> SV half-width **50 µm** (→ 100 µm SV)
$ sed -n 95,99p ...                                                  -> other widths not quoted; no claim about spectra
$ sed -n 843,848p etna/api/simulation.py                             -> sv_width_um: float = 150.0
$ grep mean_chord_length etna/api/simulation.py                      -> l_bar = mean_chord_length(sv_thickness_um, sv_width_um=sv_width_um)
$ sed -n 55,59p app/workflows/microdosimetry.py                      -> value=150.0
```

## Impact

A user comparing app spectra with the DESIGN-2 geometry gets a different l̄ and so a shifted y-scale.
This is a minor consistency problem.

## Proposed fix

```diff
-    sv_width_um: float = 150.0,
+    sv_width_um: float = 100.0,  # DESIGN-2 validated SV (half-width 50 um)
```

Make the same change in `app/workflows/microdosimetry.py:58` (`value=100.0`). Alternatively, keep 150 and
add a line to DESIGN-2 §caveats noting that the spectrum tooling default differs from the geometry spec.

## Test plan

`tests/test_api_microdosimetry.py`: assert that the default `l_bar_um` equals
`mean_chord_length(10.0, sv_width_um=100.0)`. Update any app test that pins 150.

## Effort

S
