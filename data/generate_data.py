"""
generate_data.py — create the 'observed outbreak' dataset locally (no downloads).

Produces data/observed_outbreak.csv: 150 days of daily new cases from ONE
stochastic epidemic in a town of N = 10,000 people, generated with a hidden
'true' R0 = 2.5 (mean infectious period 4 days, gamma = 0.25) and a fixed seed.

Participants play disease detectives: they are NOT told R0 = 2.5 — estimating
it from this curve is the first scientific task in both group notebooks.
The file is tiny and committed; this script exists for reproducibility.

Run from the repository root:  python data/generate_data.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from epidemic_helpers import simulate_stochastic  # noqa: E402

TRUE_R0 = 2.5
SEED = 11          # chosen so the seeded run "takes off" (some seeds fizzle!)


def main():
    rng = np.random.default_rng(SEED)
    new_cases, final_size = simulate_stochastic(TRUE_R0, N=10_000, I0=5,
                                                days=150, gamma=0.25, rng=rng)
    out = os.path.join(os.path.dirname(__file__), "observed_outbreak.csv")
    with open(out, "w") as f:
        f.write("day,new_cases\n")
        for d, c in enumerate(new_cases):
            f.write(f"{d},{c}\n")
    print(f"Wrote {out}")
    print(f"  final size: {final_size} of 10000 people "
          f"({100 * final_size / 10_000:.1f}%)")
    print(f"  peak day:   {int(np.argmax(new_cases))}, "
          f"peak cases: {int(new_cases.max())}/day")


if __name__ == "__main__":
    main()
