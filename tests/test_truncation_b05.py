"""audit-v7 B-05: truncated DD / C-V sweeps keep the converged points, flagged.

- cce_vs_bias_at_fluence returns the converged prefix with truncated=True instead of
  losing every point to the RuntimeError of the first unreachable bias.
- sweep_cce_vs_temperature pairs each CCE with the voltage the sweep returned (never
  by position in the requested voltages) and carries a truncated column.
"""

import numpy as np
import pytest


def test_cce_vs_bias_at_fluence_keeps_the_converged_prefix(monkeypatch):
    import etna.core.drift_diffusion as dd
    from etna.core.charge_collection import cce_vs_bias_at_fluence

    real_ramp = dd.ramp_bias
    calls = []

    def ramp(device_info, V_target, **kw):
        calls.append(V_target)
        if len(calls) > 1:
            raise RuntimeError("ramp_bias: failed to converge (test)")
        return real_ramp(device_info, V_target, **kw)

    monkeypatch.setattr(dd, "ramp_bias", ramp)
    result = cce_vs_bias_at_fluence(np.array([-20.0, -10.0]), fluence=1e12)
    assert result["truncated"] is True
    assert result["requested_v_stop"] == -20.0
    assert result["voltages"].tolist() == [-10.0]  # ramp order: -10 V first
    assert len(result["cce_values"]) == len(result["valid"]) == len(result["reason"])
    assert 0.0 < result["cce_values"][0] <= 1.0


def test_dd_temperature_sweep_pairs_cce_with_the_returned_voltages(monkeypatch):
    import etna.core.charge_collection as cc
    from etna.core.temperature_sweep import sweep_cce_vs_temperature

    def fake_cce_vs_bias(V_range, device_kwargs=None):
        # requested [-30, -10, -20]; the sweep reached -10 and -20 only
        return {
            "voltages": np.array([-10.0, -20.0]),
            "cce_values": np.array([0.5, 0.7]),
            "truncated": True,
        }

    monkeypatch.setattr(cc, "cce_vs_bias", fake_cce_vs_bias)
    df = sweep_cce_vs_temperature([300.0], [-30.0, -10.0, -20.0], method="dd")
    assert list(df.columns) == ["T", "V", "CCE", "truncated"]
    assert dict(zip(df["V"], df["CCE"])) == {-10.0: 0.5, -20.0: 0.7}
    assert df["truncated"].all()


def test_hecht_temperature_sweep_uses_the_reached_voltages(monkeypatch):
    import etna.core.cv_analysis as cva
    import etna.core.drift_diffusion as dd
    from etna.core.charge_collection import hecht_cce
    from etna.core.temperature_sweep import sweep_cce_vs_temperature

    monkeypatch.setattr(dd, "create_dd_device", lambda **kw: {"device_name": "stub"})
    monkeypatch.setattr(
        cva,
        "cv_sweep",
        lambda device_info, V: {
            "voltages": np.array([-10.0]),
            "depletion_widths": np.array([5e-4]),
            "truncated": True,
        },
    )
    df = sweep_cce_vs_temperature([300.0], [-10.0, -20.0], method="hecht")
    assert df["V"].tolist() == [-10.0]
    assert df["CCE"].iloc[0] == pytest.approx(float(hecht_cce(-10.0, d=5e-4, T=300.0)))
    assert df["truncated"].tolist() == [True]
