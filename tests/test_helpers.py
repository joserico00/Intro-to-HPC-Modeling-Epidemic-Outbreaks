"""Offline checks for epidemic_helpers.py and the reference numbers in expected_outputs/.

No downloads and no notebooks: the committed outbreak record and the helper functions are
enough, so a change in the modelling code fails here rather than in a participant's notebook.

Run:  python tests/test_helpers.py
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from epidemic_helpers import (  # noqa: E402
    align_to_threshold,
    analytic_final_size,
    growth_rate,
    rmse,
    run_many,
    simulate_sir,
    simulate_stochastic,
    sir_day,
)

N = 10_000
I0 = 5
TRUE_R0 = 2.5
ok = 0


def check(label, condition):
    global ok
    assert condition, f"FAILED: {label}"
    ok += 1
    print(f"  ok  {label}")


observed = np.loadtxt(REPO_ROOT / "data" / "observed_outbreak.csv", delimiter=",", skiprows=1)[:, 1]
reference = json.loads((REPO_ROOT / "expected_outputs" / "reference_metrics.json").read_text())

print("the recorded outbreak")
check("150 days of daily case counts", len(observed) == 150)
check("its final size is the reference value", int(observed.sum()) + I0 == reference["observed_final_size"])
regenerated, final_size = simulate_stochastic(TRUE_R0, N=N, I0=I0, days=150, gamma=0.25,
                                              rng=np.random.default_rng(11))
check("data/generate_data.py reproduces the committed record exactly",
      np.array_equal(regenerated, observed.astype(int)) and final_size == reference["observed_final_size"])

print("the equation model (Group A's tools)")
S, I, R, new_cases = simulate_sir(TRUE_R0)
check("nobody is created or lost", np.allclose(S + I + R, N))
check("one day of the model moves people from S to I to R",
      sir_day(9_000.0, 100.0, 900.0, TRUE_R0 * 0.25, 0.25, N)[0] < 9_000)
check("an epidemic at R0 = 2.5 takes off and burns out", new_cases.sum() > 5_000 and I[-1] < 1)
check("an epidemic at R0 below 1 never takes off", simulate_sir(0.8)[3].sum() < 100)
check("vaccinating above the herd-immunity threshold prevents it, below it does not",
      simulate_sir(TRUE_R0, vaccinated_frac=0.7)[3].sum() < 100
      and simulate_sir(TRUE_R0, vaccinated_frac=0.3)[3].sum() > 1_000)

print("measuring R0 the way Group A does")
obs_aligned = align_to_threshold(observed, 20)
candidates = [1.5, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5]
errors = [rmse(obs_aligned, align_to_threshold(simulate_sir(R0)[3], 20)) for R0 in candidates]
best_R0 = candidates[int(np.argmin(errors))]
check(f"the best-fitting candidate is the reference value (got {best_R0})", best_R0 == reference["groupA_R0"])
check("the herd-immunity formula gives the reference threshold",
      round(1 - 1 / reference["groupB_R0"], 2) == reference["threshold_formula"])

print("the coin-flip model (Group B's tools)")
curves, finals = run_many(TRUE_R0, n_runs=50, seed=0)
check("every run is recorded", curves.shape == (50, 150) and len(finals) == 50)
check("the same seed gives the same runs", np.array_equal(run_many(TRUE_R0, n_runs=5, seed=0)[1],
                                                          finals[:5]))
check("most runs take off and some fizzle", (finals > 1_000).sum() > 25 and (finals < 500).sum() > 0)
rng = np.random.default_rng(4)
fizzled = sum(1 for _ in range(500) if simulate_stochastic(TRUE_R0, I0=1, rng=rng)[1] < 500)
check(f"a single introduction fizzles about as often as the reference ({fizzled / 500:.3f})",
      abs(fizzled / 500 - reference["p_fizzle_one_case"]) < 0.1)

print("curve tools")
check("alignment starts the curve at the 20th case", np.cumsum(obs_aligned)[0] >= 1
      and obs_aligned.sum() < observed.sum())
check("a curve that never reaches the threshold aligns to an empty array",
      align_to_threshold([1, 1, 1], 20).size == 0)
check("a curve compared with itself has no error", rmse(obs_aligned, obs_aligned) == 0)
check("growth is positive while the outbreak climbs", growth_rate(obs_aligned) > 0)
check("the final-size equation agrees with the simulated epidemic",
      abs(analytic_final_size(TRUE_R0) - new_cases.sum() / N) < 0.05)

print(f"\nAll {ok} checks passed.")
