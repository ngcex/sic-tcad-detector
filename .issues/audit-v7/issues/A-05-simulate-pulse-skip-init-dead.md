---
id: A-05
title: simulate_pulse(skip_init=...) is a dead parameter; per-pulse I_dark in simulate_pulse_train includes residual carriers from the previous pulse
severity: medium
category: bug
status: fixed 2026-10-18 (WP-A 03c7c17: initialize() stores the dark baseline, pulse train reuses it)
location: etna/core/transient.py:232 (param), 255-258 (docstring), ~279-281 (unconditional I_dark extraction), ~506-516 (simulate_pulse_train passes skip_init=(i > 0))
labels: [audit-v7, slice-A, transient, flash, dead-parameter]
source: raw-A-core-solvers.md F5
related: A-11 (lists skip_init as dead code)
---

## Summary

The docstring says that `skip_init=True` skips the dark-current measurement. The function
body never reads `skip_init`. `I_dark = extract_contact_current(...)` always runs at the
start of each pulse.

For pulse i>0 in a train, the device still carries the tail of pulse i-1. Each
`pulse_times[i]["I_dark"]` is therefore the decaying signal current, not the dark current.
The train-level `I_dark` (measured once before the loop) is correct.

## Reproduction / Evidence

- `grep -rnw skip_init etna app tests scripts notebooks` returns only
  transient.py:232, :255 and :515. That is the definition, the docstring and the
  pass-through, with no read in the body.
- I read the code at transient.py ~279: `I_dark = extract_contact_current(self.device_info, contact=self.contact)`
  runs unconditionally.
- simulate_pulse_train (~493-516) measures the train-level I_dark once, then calls
  `simulate_pulse(..., skip_init=(i > 0))` while keeping device state.

## Impact

`compute_transient_cce(result_i, ...)` on pulse i>0 subtracts an inflated baseline and
under-reports CCE for later pulses. That defeats the inter-pulse memory analysis the train
function is designed for. Per-pulse CCE is used only by notebooks 08 and 08b, so the
impact is currently latent (from the raw audit, not re-run).

## Proposed fix

```python
# TransientSolver
def initialize(self):
    ...
    self._I_dark = extract_contact_current(self.device_info, contact=self.contact)

def simulate_pulse(..., skip_init=False):
    if skip_init and getattr(self, "_I_dark", None) is not None:
        I_dark = self._I_dark
    else:
        I_dark = extract_contact_current(self.device_info, contact=self.contact)
        self._I_dark = I_dark
```

## Test plan

- Monkeypatch `extract_contact_current` with a call counter. Calling
  `simulate_pulse(skip_init=True)` after `initialize()` must not perform the baseline
  measurement (count = number of time steps only).
- In a 2-pulse train, `pulse_times[1]["I_dark"] == pulse_times[0]["I_dark"] == result["I_dark"]`.

## Effort

S
