"""Measurement inputs: load owner-supplied C-V, SIMS/SRP and I(T) files, validate them,
compare them with the model and turn them into model inputs.

Nothing here changes the physics: the model side is the public ``etna.run_cv`` /
``run_cv_batch`` / ``run_dark_current`` facades on a ``DeviceConfig``. The repository
ships header-only templates (``data/measurements/templates/``, columns and units in
``data/measurements/README.md``) and no measurement data.

- C-V (``etna.measurements.cv``): ``load_cv_csv`` -> ``CVData`` (per-area C, W = eps*A/C),
  ``compare_cv`` (rms of W and C vs V), ``depletion_profile`` (differential N(W)),
  ``invert_cv`` (the lab/cv_gap knot inversion -> ``donor_table``).
- SIMS/SRP (``etna.measurements.sims``): ``load_sims_csv`` -> ``DepthProfile``,
  ``resample_profile`` -> ``ProfileTable`` (donor_table + per-knot uncertainty, p+ N_A
  and junction depth), ``ProfileTable.to_device_config``.
- I(T) (``etna.measurements.iv_t``): ``load_it_csv`` -> ``ITData``, ``fit_arrhenius``,
  ``analyze_it`` (data E_a per bias vs the model's E_a by the same fit).
- Unit / plausibility checks (``etna.measurements.plausibility``): obvious unit
  mix-ups are refused with a hint, unusual values warned
  (``MeasurementPlausibilityWarning``, ``plausibility_report``).
- Outputs: ``to_json()`` on every result, ``load_donor_table_json``,
  ``report_markdown``, figures in ``etna.measurements.plots``.
"""

from etna.measurements._csv import MeasurementFileError
from etna.measurements.cv import (
    CVComparison,
    CVData,
    CVInversion,
    DepletionProfile,
    compare_cv,
    default_knots_um,
    depletion_profile,
    invert_cv,
    load_cv_csv,
)
from etna.measurements.io import (
    load_donor_table_json,
    read_donor_table_json,
    report_markdown,
    results_json,
)
from etna.measurements.iv_t import (
    ArrheniusFit,
    EaComparison,
    ITAnalysis,
    ITData,
    analyze_it,
    fit_arrhenius,
    load_it_csv,
    model_current_vs_T,
)
from etna.measurements.plausibility import (
    MeasurementPlausibilityWarning,
    plausibility_report,
)
from etna.measurements.sims import (
    DepthProfile,
    ProfileTable,
    load_sims_csv,
    model_donor_profile,
    profile_vs_model,
    resample_profile,
)

__all__ = [
    "MeasurementFileError",
    "MeasurementPlausibilityWarning",
    "plausibility_report",
    "CVData",
    "CVComparison",
    "CVInversion",
    "DepletionProfile",
    "load_cv_csv",
    "compare_cv",
    "depletion_profile",
    "invert_cv",
    "default_knots_um",
    "DepthProfile",
    "ProfileTable",
    "load_sims_csv",
    "resample_profile",
    "model_donor_profile",
    "profile_vs_model",
    "ITData",
    "ArrheniusFit",
    "EaComparison",
    "ITAnalysis",
    "load_it_csv",
    "fit_arrhenius",
    "analyze_it",
    "model_current_vs_T",
    "load_donor_table_json",
    "read_donor_table_json",
    "report_markdown",
    "results_json",
]
