"""etna.design on the batched engine: evaluator parity, bounds, robust mode, facade."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

import etna
from etna import DeviceConfig, PhysicsModels
from etna.api.simulation import _ENGINE_BOX, validated_domain_reason
from etna.design import (
    DesignSettings,
    RobustSpec,
    TOLERANCES,
    area_bounds,
    design_space,
    evaluate_configs,
    make_problem,
    run_design,
)
from etna.design.dosimeter import perturb

CFGS = [
    DeviceConfig(),
    DeviceConfig(N_D_bulk=1e14, epi_thickness_um=12.0),
    DeviceConfig(N_D_bulk=5e13, N_D_junction=3e15, L_transition_um=1.0, area_cm2=0.1),
]
CHEAP = ["sensitivity_nC_per_Gy", "min:dark_current_A"]


class TestEvaluator:
    def test_parity_with_the_facades(self):
        """The vectorised evaluator reproduces run_cv / run_cce / run_dark_current /
        run_transit_time at the working bias (same engine kernels, one batch)."""
        m = evaluate_configs(CFGS, DesignSettings(V_work=-10.0))
        assert m["ok"].all()
        assert m["group"].tolist() == [0, 1, 0]  # 10 um configs share one batch
        hf = PhysicsModels(high_field_mobility="lades2000")
        for i, c in enumerate(CFGS):
            cv = etna.run_cv(c, 0.0, -10.0, 2)
            assert m["capacitance_F"][i] == pytest.approx(cv.y[-1], rel=1e-12)
            assert m["W_dep_um"][i] == pytest.approx(
                cv.metadata["depletion_widths"][-1] * 1e4, rel=1e-12
            )
            cce = etna.run_cce(c, -10.0, -10.0, 1)
            assert m["cce"][i] == pytest.approx(cce.y[-1], abs=1e-9)
            dk = etna.run_dark_current(c, 0.0, -10.0, 2)
            assert m["dark_current_A"][i] == pytest.approx(dk.y[-1], rel=1e-9)
            tt = etna.run_transit_time(c, hf, -10.0, -10.0, 1)
            t = max(tt.metadata["t_e"][0], tt.metadata["t_h"][0])
            assert m["transit_time_s"][i] == pytest.approx(t, rel=1e-12)
        # sensitivity is linear in the area at fixed epi physics
        assert m["sensitivity_nC_per_Gy"][2] / m["sensitivity_nC_per_Gy"][
            0
        ] == pytest.approx(2.5 * m["cce_uniform"][2] / m["cce_uniform"][0], rel=1e-12)

    def test_expensive_metrics_match_facades(self):
        c = CFGS[0]
        s = DesignSettings(V_work=-10.0, fluence=1e13)
        m = evaluate_configs([c], s, metrics=("cce_fluence",))
        r = etna.run_radiation_damage(
            c, fluences=np.array([1e13]), V_bias=-10.0, proton_energy_MeV=62.0
        )
        assert m["cce_fluence"][0] == r.y[0]
        assert m["solver_damage"][0] == "batched+traps"

    def test_outside_box_is_never_sent_to_devsim(self, monkeypatch):
        import etna.api.simulation as sim

        def boom(*a, **k):
            raise AssertionError("devsim path called")

        monkeypatch.setattr(sim, "_run_cv_devsim", boom)
        monkeypatch.setattr(sim, "_run_cce_devsim", boom)
        bad = [DeviceConfig(T=500.0), DeviceConfig(epi_thickness_um=40.0)]
        m = evaluate_configs(bad + CFGS[:1])
        assert m["ok"].tolist() == [False, False, True]
        assert all("no devsim" in r for r in m["reason"][:2])
        assert np.isnan(m["cce"][:2]).all()


class TestBounds:
    def test_bounds_come_from_the_engine_box(self):
        sp = design_space()
        v = {x.name: x for x in sp.variables}
        assert (v["N_D_bulk"].lo, v["N_D_bulk"].hi) == _ENGINE_BOX["N_D_bulk"]
        assert (v["epi_thickness_um"].lo, v["epi_thickness_um"].hi) == _ENGINE_BOX[
            "epi_thickness_um"
        ]
        assert v["epi_thickness_um"].step == 0.25 and v["N_D_bulk"].log
        sp = design_space(bounds={"epi_thickness_um": (8.0, 12.0)})
        assert sp.variables[0].lo == 8.0 and sp.variables[0].hi == 12.0

    def test_robust_bounds_keep_every_perturbation_inside_the_box(self):
        sp = design_space(robust=True)
        v = {x.name: x for x in sp.variables}
        assert v["N_D_bulk"].lo == pytest.approx(2e13 / 0.8)
        assert v["L_transition_um"].lo == pytest.approx(0.6)
        assert v["epi_thickness_um"].lo == 6.0 and v["epi_thickness_um"].hi == 29.0
        names, D, _, _ = RobustSpec().deltas()
        corners = sp.decode(
            np.array([[0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [0, 1, 0, 1, 0]], float)
        )
        for row in corners:
            c = dataclasses.replace(DeviceConfig(), **dict(zip(sp.names, row)))
            for d in D:
                assert validated_domain_reason(perturb(c, names, d, TOLERANCES)) is None

    def test_area_bound_is_flagged_placeholder(self):
        lo, hi, note = area_bounds()
        assert lo < 0.04 < hi  # the documented 4 mm2 dosimeter fits
        assert "PLACEHOLDER" in note


class TestRunDesign:
    def test_facade_pareto_and_figure(self):
        r = run_design(
            objectives=CHEAP, n_init=6, n_iter=1, batch_size=4, seed=0,
            breakdown_check=False, bounds={"epi_thickness_um": (6.0, 14.0)},
        )  # fmt: skip
        md = r.metadata
        assert r.sim_type == "design" and md["solver"] == "batched"
        assert md["n_evaluations"] == 10 and md["n_feasible"] >= 2
        assert len(r.x) == len(md["pareto"]) == len(md["pareto_configs"]) >= 1
        sens = [p["metrics"]["sensitivity_nC_per_Gy"] for p in md["pareto"]]
        assert np.allclose(r.x, sens)
        assert set(md["alt_physics"]) == set(etna.design.ALT_PHYSICS)
        assert md["notes"] and md["throughput_designs_per_s"] > 0
        from etna.viz import figure_for

        assert figure_for(r).etna_kind == "design"

    def test_single_objective_with_constraint(self):
        r = run_design(
            objectives=["min:V_fd"], constraints=["sensitivity_nC_per_Gy>=14"],
            n_init=6, n_iter=2, batch_size=3, seed=1, alt_physics=False,
            breakdown_check=False,
        )  # fmt: skip
        b = r.metadata["best"]
        assert b is not None and b["metrics"]["sensitivity_nC_per_Gy"] >= 14
        assert r.y[-1] == pytest.approx(b["metrics"]["V_fd"])
        assert np.all(np.diff(r.y[np.isfinite(r.y)]) <= 0)  # best-so-far

    def test_robust_worst_case_is_no_better_than_nominal(self):
        rs = RobustSpec(mode="worst", n_samples=2)
        pr = make_problem(objectives=[etna.design.Objective("sensitivity_nC_per_Gy", "max")],
                          constraints=[], robust=rs)  # fmt: skip
        P = pr.space.decode(
            np.array([[0.3, 0.4, 0.5, 0.5, 0.5], [0.6, 0.2, 0.3, 0.7, 0.2]])
        )
        m = pr.evaluate(P)
        assert m["ok"].all()
        assert np.all(m["sensitivity_nC_per_Gy"] <= m["sensitivity_nC_per_Gy_nominal"])
        assert np.all(
            m["dark_current_A"] == m["dark_current_A_nominal"]
        )  # not an objective

    def test_microdosimeter_problem_with_mocked_cell(self, monkeypatch):
        import etna.design.micro as micro

        def fake(params, base, V_work, n_scan=8, cce_threshold=0.95):
            ro, rc = params["outer_radius_um"], params["column_radius_um"]
            return {
                "efficient_fraction": 1.0 - (rc + 3.0) ** 2 / ro**2,
                "cce_min": 0.8, "cce_mean_area": 0.9,
                "cell_area_um2": np.pi * ro**2, "scan_spacing_um": ro / (n_scan - 1),
            }  # fmt: skip

        monkeypatch.setattr(micro, "evaluate_cell", fake)
        r = run_design(
            problem="microdosimeter",
            n_init=6,
            n_iter=1,
            batch_size=2,
            seed=0,
            cache=False,  # the mocked cell must not enter the process-wide cache
        )
        md = r.metadata
        assert md["problem"] == "microdosimeter" and md["solver"] == "fv2d"
        assert md["pareto"] and md["pareto_configs"] == []
        assert "annular-trench" in md["notes"][0]

    def test_damage_objective_without_cache(self, tmp_path):
        """String per-point metrics (solver_damage) survive the uncached path and the
        checkpoint (regression: they broke the Pareto rows)."""
        kw = dict(
            objectives=["sensitivity_nC_per_Gy", "cce_fluence"], constraints=[],
            n_init=3, n_iter=0,
            seed=2, cache=False, alt_physics=False, breakdown_check=False,
            bounds={"epi_thickness_um": (8.0, 10.0)}, checkpoint=tmp_path / "d.npz",
        )  # fmt: skip
        r = run_design(method="random", **kw)
        md = r.metadata
        assert md["n_evaluations"] == 3 and md["pareto"]
        assert set(md["history"]["metrics"]["solver_damage"]) <= {"batched+traps", ""}
        assert "solver_damage" not in md["pareto"][0]["metrics"]
        r2 = run_design(method="random", **kw)  # resumed from the checkpoint
        assert np.array_equal(r2.metadata["history"]["P"], md["history"]["P"])
