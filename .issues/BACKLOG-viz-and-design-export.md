# Backlog: visualisation and engineering export (start ONLY after the n+ default flip is finished and validated by TCAD)

Authorised by the user (2026-10-17): uv env may be changed to add the optional deps.

1. Spec sheet of the stack (JSON/CSV/PDF) from DeviceConfig + run_all: layers, doping, area,
   bias, V_fd and reach-through (with definition), CCE, dark current. + dimensioned
   cross-section (SVG/PNG).
2. Manim animations (optional dep, `uv add --optional viz manim`; macOS may need
   `brew install pkg-config cairo pango`; ffmpeg is present): scene 1 stack in scale,
   scene 2 depletion growth + E(x) vs bias with real run_field/run_cv data. App button
   "Export animation" + st.video; Plotly stays the live view.
3. GDSII layout (gdstk) of active area / p+ edge / guard ring / pads. Label "reference
   layout, not TCAD-validated"; needs real dimensions from the owner (guard ring and pads
   are hand-scored in DESIGN-1/2). Do not invent dimensions.
4. Mesh/profile export (CSV/VTK) for cross-checking with other TCAD tools.
5. Validation variants (advisor order): microdosimeter preset with tilted-alpha CCE vs bias
   (Micro Fig. 4b, shape only), 1D graded epi/n+ buffer test for the W creep, finite 2D p+
   implant with oxide/field plate (edge C-V, guard ring).
