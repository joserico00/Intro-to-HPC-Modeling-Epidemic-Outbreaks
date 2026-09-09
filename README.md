# Intro to HPC Bootcamp: Modeling Epidemic Outbreaks

> **Research question:** When does an outbreak explode — and how much
> vaccination stops it?

A disease swept through a town of 10,000 people. Participants get only the
public-health record — 150 days of daily case counts — and play disease
detectives: **measure the outbreak's R₀**, then **find the vaccination level
that would have prevented it**, and check it against the century-old
herd-immunity formula **v\* = 1 − 1/R₀**.

Two groups answer the same question with opposite modeling philosophies, and
the joint presentation's payoff is watching them agree:

| | Group A — the equation model | Group B — the coin-flip model |
|---|---|---|
| Open | `Group_A/01_epidemics_group_A.ipynb` | `Group_B/01_epidemics_group_B.ipynb` |
| Model | SIR difference equations (smooth averages) | stochastic simulation (whole people, real luck) |
| Builds from scratch | `sir_day` + `simulate_sir` (~20 lines) | `simulate_outbreak` (~18 lines) |
| Special insight | a *sharp* vaccination cliff at 1 − 1/R₀ | a *soft* cliff + stochastic extinction (32% of single-case introductions fizzle on luck alone) |

Project leads and curious students: `01_epidemics_complete_master.ipynb` is the
thorough deep dive — both methods, a 🧭 Python survival guide, 🧩 function maps,
and master-only mathematics (the final-size equation, the growth-rate shortcut
and its bias).

**No Python experience required.** Every code cell is preceded by a
plain-language "What the next cell does" note; the only cell participants edit
is clearly marked ✏️; helper machinery is hidden 🔒 and documented in the
master's function map.

## Reading without Jupyter

`notebook-preview.html` renders all three notebooks — with real executed
outputs and figures — in any browser: tabs for Group A, Group B, and the Master
deep dive. `master-vs-groups.html` shows how the master maps onto the
simplified group paths.

## Quick start (runs anywhere)

```bash
cd <participant_work_directory>
git clone <repository_url>
cd Intro-to-HPC-Bootcamp-Modeling-Epidemic-Outbreaks

# conda/mamba:
conda env create -f environment.yml && conda activate epidemic-modeling
# or venv/pip:
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt

python -m ipykernel install --user --name epidemic-modeling \
    --display-name "Python (epidemic-modeling)"
jupyter lab      # open your group's notebook; pick the epidemic-modeling kernel
```

On an HPC site, follow the current site's documentation for the module command,
approved environment directory, and Jupyter portal — confirm names with your
mentor, never assume last year's. No batch jobs needed: everything runs
in-kernel in seconds. (Group B's many-simulations design *is* the HPC array-job
pattern, and the master's §10 says how it scales.)

## Expected results (reference values)

| Measurement | Group A (equations) | Group B (coin flips) |
|---|---|---|
| Fitted R₀ (hidden truth: 2.5) | ≈ 2.4 | ≈ 2.5 |
| Vaccination threshold | ≈ 0.60–0.65 | ≈ 0.60 |
| Herd-immunity formula 1 − 1/R₀ | 0.58 | 0.60 |
| P(fizzle from one case) | — (equations can't) | ≈ 0.32 (theory 1/R₀ = 0.40) |

Reference figures and metrics live in `expected_outputs/`.

## Data provenance

The outbreak is synthetic: [`data/generate_data.py`](data/generate_data.py)
runs one seeded stochastic epidemic (hidden truth R₀ = 2.5, mean infectious
period 4 days). Fully reproducible, license-free, no real patients.
Real-data extensions: WHO measles line lists, Our World in Data COVID curves.

## Time budget

The whole pathway — notebook, experiments, and presentation preparation — fits
the program's **10–12 hour** participant total; the headline figure (the
vaccination cliff) is produced about two-thirds through the notebook.

## Design lineage

Structure and pedagogy follow the Intro to HPC Bootcamp project-lead guide and
its sibling projects (regulatory genomics, chaotic-weather forecasting, urban
air quality): thorough master at the root, simplified Group A/B copies, one ✏️
edit area, ✅ checkpoints, 📊 variable reference, and a shared final experiment.
