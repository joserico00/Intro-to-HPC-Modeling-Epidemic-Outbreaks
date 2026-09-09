"""
epidemic_helpers.py — tested helper code for the "Modeling Epidemic Outbreaks"
bootcamp project.

This module hides routine infrastructure so the notebooks can focus on the
scientific choices. Everything is plain NumPy — no external libraries — so the
project runs on any machine with Python + NumPy.

Hidden should mean low-attention, not mysterious: every function has a short
docstring, and the thorough master notebook's 🧩 function map explains each
tool's job in plain English.
"""

from __future__ import annotations

import numpy as np

# ----------------------------------------------------------------------------
# 1. The equation-based model (Group A): SIR
# ----------------------------------------------------------------------------

def sir_day(S, I, R, beta, gamma, N):
    """One day of the SIR model: tomorrow's counts from today's counts.

    S: susceptible (could catch it), I: infectious (spreading it),
    R: recovered (immune). beta = R0 * gamma is the transmission rate;
    gamma = recovery rate (1 / average days infectious).
    The infection probability 1 - exp(-beta*I/N) is exactly the average of
    the coin-flip rule used by the stochastic (Group B) model — the two
    models are deterministic and random twins.
    """
    new_infections = S * (1.0 - np.exp(-beta * I / N))
    recoveries = gamma * I
    return (S - new_infections,
            I + new_infections - recoveries,
            R + recoveries,
            new_infections)


def simulate_sir(R0, N=10_000, I0=5, vaccinated_frac=0.0, days=150,
                 gamma=0.25):
    """Run the daily SIR equations; returns arrays (S, I, R, new_cases).

    R0 : average people one sick person infects in a fully susceptible crowd.
    vaccinated_frac : fraction of the population immune from day 0.
    """
    beta = R0 * gamma
    S = np.empty(days + 1); I = np.empty(days + 1); R = np.empty(days + 1)
    new_cases = np.empty(days)
    S[0] = N * (1 - vaccinated_frac) - I0
    I[0] = I0
    R[0] = N * vaccinated_frac
    for d in range(days):
        S[d + 1], I[d + 1], R[d + 1], new_cases[d] = sir_day(
            S[d], I[d], R[d], beta, gamma, N)
    return S, I, R, new_cases


# ----------------------------------------------------------------------------
# 2. The individual-based model (Group B): stochastic chain of infections
# ----------------------------------------------------------------------------

def simulate_stochastic(R0, N=10_000, I0=5, vaccinated_frac=0.0, days=150,
                        gamma=0.25, rng=None):
    """One random epidemic: whole people, coin flips, different every run.

    Each day: every susceptible person independently escapes (or not) the
    infection pressure from the I infectious people; each infectious person
    recovers with probability gamma. Returns (new_cases, final_size).
    """
    if rng is None:
        rng = np.random.default_rng()
    beta = R0 * gamma
    S = int(round(N * (1 - vaccinated_frac))) - I0
    I = I0
    new_cases = np.zeros(days, dtype=int)
    for d in range(days):
        p_infected_today = 1.0 - np.exp(-beta * I / N)
        infections = rng.binomial(S, p_infected_today)
        recoveries = rng.binomial(I, gamma)
        S -= infections
        I += infections - recoveries
        new_cases[d] = infections
        if I == 0:
            break
    final_size = int(new_cases.sum()) + I0
    return new_cases, final_size


def run_many(R0, n_runs=200, seed=0, **kwargs):
    """Repeat the random epidemic n_runs times.

    Returns (curves, final_sizes): curves is an array (n_runs, days) of daily
    new cases; final_sizes is an array (n_runs,) of total people ever infected.
    """
    rng = np.random.default_rng(seed)
    days = kwargs.get("days", 150)
    curves = np.zeros((n_runs, days), dtype=int)
    finals = np.zeros(n_runs, dtype=int)
    for i in range(n_runs):
        nc, fs = simulate_stochastic(R0, rng=rng, **kwargs)
        curves[i, :len(nc)] = nc
        finals[i] = fs
    return curves, finals


# ----------------------------------------------------------------------------
# 3. Working with outbreak curves
# ----------------------------------------------------------------------------

def align_to_threshold(new_cases, threshold=20):
    """Shift a daily-cases curve so day 0 = the day cumulative cases pass
    `threshold`. Random epidemics take off at random times; aligning them at
    the same milestone makes curves comparable."""
    cum = np.cumsum(new_cases)
    idx = np.argmax(cum >= threshold)
    if cum[idx] < threshold:                 # never reached the threshold
        return np.array([])
    return np.asarray(new_cases)[idx:]


def growth_rate(new_cases, start=5, end=25):
    """Early exponential growth rate r: the slope of log(daily cases) during
    the growth phase (days `start`..`end` after alignment).
    R0 can then be estimated as 1 + r / gamma."""
    days = np.arange(start, end)
    y = np.asarray(new_cases, dtype=float)[start:end]
    good = y > 0
    slope, _ = np.polyfit(days[good], np.log(y[good]), 1)
    return slope


def rmse(a, b):
    """Typical difference between two curves (compared over the shorter one)."""
    n = min(len(a), len(b))
    return float(np.sqrt(np.mean((np.asarray(a[:n], float)
                                  - np.asarray(b[:n], float)) ** 2)))


def analytic_final_size(R0, tol=1e-10):
    """Solve the classic final-size equation  Z = 1 - exp(-R0 * Z)  for the
    fraction of the population ultimately infected (master notebook only)."""
    Z = 0.9
    for _ in range(200):
        Z_new = 1.0 - np.exp(-R0 * Z)
        if abs(Z_new - Z) < tol:
            break
        Z = Z_new
    return Z
