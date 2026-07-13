# Phase 40 Plan Verification — Geometry Viewer

**Verified:** 2026-07-13
**Plans checked:** 40-01-PLAN.md, 40-02-PLAN.md
**Verdict:** ISSUES FOUND (1 blocker — mechanical/process gate; substantive plan content passes)

---

## Summary

Goal-backward review of the two Phase 40 plans against the phase goal ("users can see
the electric field overlaid on a 2D device cross-section... adapting automatically to
1D/2D") and the three requirements VIZ-01/VIZ-02/VIZ-03. The substantive plan content —
requirement coverage, the three integration seams flagged in RESEARCH.md, dependency
graph, scope, and locked-decision compliance — all check out against source. One
process-gate blocker (missing `40-VALIDATION.md`) and one warning (unresolved Open
Questions heading) were found and must be addressed before/alongside execution.

---

## Dimension 1: Requirement Coverage — PASS

| Requirement | Plans | Tasks | Status |
|---|---|---|---|
| VIZ-01 (2D heatmap) | 40-01, 40-02 | 40-01 T1/T2 (builder+test), 40-02 T1/T2 (wiring+2D-route test) | Covered |
| VIZ-02 (1D bar, same interface) | 40-01, 40-02 | 40-01 T1/T2, 40-02 T1 (branch keeps line charts + adds bar) | Covered |
| VIZ-03 (quantity dropdown, no re-solve) | 40-01, 40-02 | 40-01 T1 (QUANTITIES), 40-02 T1 (selectbox) + T2 (call-counter test) | Covered |

Both plans declare `requirements: [VIZ-01, VIZ-02, VIZ-03]` in frontmatter. All three
requirement IDs from ROADMAP.md Phase 40 appear in both plans' frontmatter. No gaps.

## Dimension 2: Task Completeness — PASS

All 4 tasks (2 per plan) have `<files>`, `<action>`, `<verify automated>`, and `<done>`.
Actions are concrete (exact dict contents, exact function bodies, exact grep/ast gates),
not vague. `<verify>` commands are runnable shell one-liners tied to concrete assertions
(import success, grep counts, `ast.parse`, `pytest -x`).

## Dimension 3: Dependency Correctness — PASS

- 40-01: `depends_on: []` → Wave 1. Correct — it only creates new files, no dependency on the page.
- 40-02: `depends_on: ["40-01"]` → Wave 2. Correct — 40-02 imports `build_geometry_figure`/`QUANTITIES` from the module 40-01 creates. No cycle, no forward reference, wave numbers consistent with dependency direction.

## Dimension 4: Key Links Planned — PASS

- `geometry_viewer.py` → `petringa.MeshData` (type-hint import): confirmed `MeshData` is
  re-exported at top-level (`petringa/__init__.py:7,24` — `from petringa.api.results import MeshData, SimResult`; `__all__` includes `"MeshData"`). Plan's `from petringa import MeshData` will work.
- `field_map.py` → `build_geometry_figure(result.mesh, quantity)`: 40-02 Task 1 step 4-6 explicitly imports and calls it; verify gate greps for the literal call string. Wired, not just declared.
- `field_map.py` → `session_state['field_result']`: verified present in current file (line 42, 50) and plan preserves the caching mechanic; 40-02 Task 2 adds a call-counter test that empirically proves no re-solve, which is the strongest possible check.

## Dimension 5: Scope Sanity — PASS

- 40-01: 2 tasks, 2 files (new module + new test file). Within budget.
- 40-02: 2 tasks, 2 files (one page edit, one test file edit). Within budget.
No plan exceeds the 2-3 task target or the file-count thresholds.

## Dimension 6: Verification Derivation — PASS

`must_haves.truths` in both plans are behavior-observable at the builder/page level
("A pure builder converts a 2D MeshData into a Plotly heatmap trace", "A 2D device
config no longer hits st.stop()...", "quantity dropdown... re-renders... without
re-running run_field"). Artifacts map to truths; `key_links` cover the two load-bearing
wiring points (MeshData import, build_geometry_figure call site).

## Dimension 7: Context Compliance — N/A

No CONTEXT.md exists for Phase 40 (confirmed by RESEARCH.md's own `<user_constraints>`
section: "No CONTEXT.md exists... discuss-phase was not run"). RESEARCH.md's Assumptions
Log (A1–A7) substitutes for locked decisions and both plans explicitly cite and adopt
A1/A2/A3/A4/A5/A6/A7 by ID. No contradiction found; no deferred-idea scope creep (lateral
slice picking and region/contact overlays are explicitly excluded from both plans, matching RESEARCH's "Deferred Ideas").

## Dimension 7b: Scope Reduction Detection — PASS (no reduction found)

Scanned both plans for reduction language ("v1", "static for now", "will be wired
later", "stub", etc.). None found describing the three VIZ requirements themselves.
The only "skip"/"out of scope" items (2D CSV download, region/contact overlay,
lateral-slice picking) are RESEARCH-flagged discretionary omissions with no
corresponding requirement (VIZ-01/02/03 do not mention CSV export or overlays) —
these are legitimate scope boundaries, not silent reductions of an in-scope requirement.

## Dimension 7c: Architectural Tier Compliance — PASS

RESEARCH.md's Architectural Responsibility Map is present and both plans respect it:
physics/mesh-extraction stays in the API tier (`run_field`, untouched by this phase);
interpolation + figure construction is correctly placed in the app/UI tier
(`geometry_viewer.py`); no task pushes solve or interpolation logic into the wrong tier.
The "never call devsim directly" constraint (security/architecture-relevant given it's
a locked decision) is enforced by an explicit grep gate in 40-01's acceptance criteria.

## Dimension 8: Nyquist Compliance

RESEARCH.md contains a `## Validation Architecture` section (confirmed, lines 351-390),
so this dimension is **NOT skipped**.

### Check 8e — VALIDATION.md Existence (Gate) — **BLOCKING FAIL**

```
$ ls .planning/phases/40-geometry-viewer/*-VALIDATION.md
No such file
```

`40-VALIDATION.md` does not exist in the phase directory (only `40-01-PLAN.md`,
`40-02-PLAN.md`, `40-PATTERNS.md`, `40-RESEARCH.md` are present). Per the plan-checker
contract, this is a **BLOCKING FAIL** regardless of the quality of RESEARCH.md's inline
"Validation Architecture" section — the dedicated `VALIDATION.md` artifact is missing.

**Remediation:** Re-run `/gsd:plan-phase 40 --research` (or the phase's validation-doc
generation step) to materialize `40-VALIDATION.md` from the RESEARCH.md content that
already exists (the substance — test framework, requirements→test map, AppTest accessor
facts — is already written in RESEARCH.md §Validation Architecture and only needs to be
extracted/promoted to the dedicated file).

### Checks 8a–8d — informational (evaluated against RESEARCH.md's inline content, pending 8e remediation)

- **8a (automated verify presence):** Every `auto` task in both plans has an
  `<automated>` command. No task relies on a Wave-0 "MISSING" placeholder needing a
  cross-reference. PASS.
- **8b (feedback latency):** All automated commands are `uv run python -c "..."` (fast,
  <1s) or `uv run pytest tests/test_app_*.py -x` (single-file, fast, no full-suite/watch
  flags, no E2E browser tooling). PASS.
- **8c (sampling continuity):** 4 tasks total across 2 waves, each with an automated
  verify. No window of 3 consecutive unverified tasks. PASS.
- **8d (Wave 0 completeness):** No `<automated>MISSING</automated>` references exist in
  either plan requiring a Wave 0 test-first task. N/A/PASS.

**Dimension 8 overall: FAIL** (blocked solely by 8e — the missing dedicated file).

## Dimension 9: Cross-Plan Data Contracts — PASS

Only one shared data path: `MeshData` produced by `run_field` (unmodified, outside this
phase) and consumed by both `geometry_viewer.py` (40-01) and `field_map.py` (40-02). No
transformation conflicts — 40-01 defines the pure read-only transform, 40-02 only calls
it; neither plan mutates `MeshData` in place. No preservation-mechanism concern applies.

## Dimension 10: CLAUDE.md Compliance

No `./CLAUDE.md` exists in the working directory (not found). **Dimension 10: SKIPPED
(no CLAUDE.md found).** (Project conventions in this repo are instead carried via
STATE.md decisions, which both plans cite correctly — e.g. `uv run pytest`, module-attribute
facade referencing for monkeypatch mockability, DeviceConfig unhashable → no `@st.cache_data`.)

## Dimension 11: Research Resolution — **WARNING**

RESEARCH.md line 415 has a `## Open Questions` heading **without** the required
`(RESOLVED)` suffix. Per contract this must be flagged even though the two questions
listed have clear, actionable dispositions written inline:

1. "Does 2D `run_field` converge for a realistic 2D DeviceConfig?" — disposition given:
   build+test against synthetic MeshData; treat browser 2D verification as best-effort;
   a non-convergence is upstream physics, not a Phase 40 defect. Both plans correctly
   adopt this disposition (test against synthetic `MeshData`, never call `run_field`/devsim
   in tests — enforced by explicit no-devsim/no-run_field grep gates in both plans' acceptance criteria).
2. "Should regions/contacts be overlaid on the heatmap?" — disposition given: out of
   scope, correctly excluded from both plans.

Since both questions have de facto resolutions that the plans correctly implement, this
is scored **WARNING**, not blocker — but the heading should be updated to
`## Open Questions (RESOLVED)` for process hygiene and to prevent future re-litigation.

## Dimension 12: Pattern Compliance — PASS

PATTERNS.md exists and maps all 4 files (2 new, 2 modified) to analogs. Both plans
reference the correct analogs verbatim:
- `geometry_viewer.py` → structural analog `results.py` (purity discipline, signature
  shape) correctly cited, with the correct caveat that the griddata/heatmap body itself
  has "No Analog Found" and must come from RESEARCH §3's verified sketch instead — both
  plans respect this distinction explicitly in their `<read_first>` blocks.
- `field_map.py` → self-edit analog, correctly surgical (delete guard, branch render, add
  selectbox) rather than a rewrite.
- Shared pattern (persistent-key selectbox from `device_sidebar.py:69-71`) is referenced
  and correctly mirrored (`index=0`, `key=`).

---

## Verification of the Three Flagged Integration Seams (explicit focus per orchestrator prompt)

**1. Removing the `st.stop()` 2D guard.** Confirmed present in current `field_map.py`
(lines 32-38, verified by direct read). 40-02 Task 1 step 1 explicitly deletes it; the
task's own `<verify automated>` asserts `'1D-only' not in src`, and acceptance criteria
add a redundant `grep -c '1D-only' ... returns 0` check. This is a hard gate, not just an
instruction — a plan that failed to delete the string would fail its own automated verify. Adequately covered.

**2. Skipping the empty-array line-chart builder for 2D.** Confirmed by reading
`petringa/api/simulation.py:301-311`: 2D `run_field` returns `x_out=np.array([])`,
`y_out=np.array([])`. `build_field_figures` (results.py:70-85) and `to_csv_bytes`
(results.py:130-141, field branch) both index `result.x`/`result.metadata["potential"]`/`result.metadata["net_doping"]` directly — for 2D, `metadata["potential"]`/`metadata["net_doping"]` are actually full node arrays (not empty, see simulation.py:295-299) while `result.x`/`result.y` are empty, so calling these builders on a 2D result would either silently plot nothing (empty x with non-empty y → shape mismatch/empty chart) or raise inside `pd.DataFrame` construction (`to_csv_bytes`, column-length mismatch: `depth_um=result.x` empty vs `ElectricField_V_per_cm=result.y` empty vs `NetDoping_cm-3=result.metadata["net_doping"]` full-length — actually column lengths for the field branch would be inconsistent between x/y and the metadata-derived columns). 40-02 Task 1 step 5 correctly branches on `result.mesh.y_coords is None` and explicitly excludes both builders from the 2D branch. This is the correct and sufficient fix — verified against source, not just plan prose.

**3. The "no re-solve on dropdown change" mechanic verified with something stronger than
"no exception".** 40-02 Task 2 step 5 (`test_selectbox_change_does_not_resolve`)
mandates a **call counter** that must remain at exactly 1 after `at.selectbox[0].select("Net doping").run()` — this is a behavioral proof of no-re-solve, not merely an absence-of-crash check. The task's acceptance criteria make this explicit: "A selectbox-no-resolve test exists using a call counter that stays at 1". This satisfies the orchestrator's specific concern; it is a materially stronger check than `at.exception == []` alone (which the plan also keeps, but not as the sole evidence).

All three seams are correctly identified and their fixes are traceable to concrete,
gated task steps — not aspirational prose.

---

## Issues

```yaml
issues:
  - dimension: nyquist_compliance
    severity: blocker
    check: "8e — VALIDATION.md existence"
    description: "40-VALIDATION.md does not exist in .planning/phases/40-geometry-viewer/ even though RESEARCH.md contains a ## Validation Architecture section, which is the trigger condition for this dimension's mandatory gate."
    plan: null
    fix_hint: "Re-run /gsd:plan-phase 40 --research (or the phase's validation-doc generation step) to materialize 40-VALIDATION.md, promoting the existing RESEARCH.md §Validation Architecture content (test framework, requirements->test map, AppTest accessor facts) into the dedicated file."

  - dimension: research_resolution
    severity: warning
    description: "RESEARCH.md's '## Open Questions' heading (line 415) lacks the required (RESOLVED) suffix, even though both listed questions have inline dispositions that both plans correctly implement."
    file: "40-RESEARCH.md"
    unresolved_questions:
      - "Does 2D run_field converge for a realistic 2D DeviceConfig? (disposition given: build/test against synthetic MeshData; browser 2D verification best-effort)"
      - "Should regions/contacts be overlaid on the heatmap? (disposition given: out of scope)"
    fix_hint: "Update the heading to '## Open Questions (RESOLVED)' to reflect that both items have actionable, plan-adopted dispositions."
```

## Additional Note (non-blocking, must be stated per advisor review)

VIZ-01's success criterion ("a 2D heatmap appears... after running a field map
simulation") is verified in these plans only via a **synthetic `MeshData` + mocked
`run_field`** (both the 40-01 pure unit test and 40-02's AppTest fakes never exercise a
real 2D devsim solve — by design, per RESEARCH's documented open risk that 2D `ramp_bias`
convergence is unverified this session). This is a reasonable and correctly-scoped
posture for a UI-wiring phase (consistent with Phase 39's precedent for the 1D solver
non-convergence issue), and the plans are honest about it rather than silently
overclaiming. It is called out here so execution/verification downstream does not treat
"tests pass" as proof that a real 2D field-map run renders a heatmap in the browser —
that remains an open, upstream-owned risk, not a Phase 40 plan defect.

---

## Recommendation

Fix the Dimension 8 Check 8e blocker (create `40-VALIDATION.md`) before proceeding to
execution — this is a process/gate requirement, not a rewrite of plan content. The
substantive plan content (task actions, requirement coverage, the three integration
seams, dependency graph, scope) requires no changes. Recommend also updating the
RESEARCH.md Open Questions heading to `(RESOLVED)` in the same pass, since the
underlying work is already done.
