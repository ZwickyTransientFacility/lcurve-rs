#!/usr/bin/env python3
"""
100-parameter-set regression test: Rust vs C++ lroche.

Generates 100 model parameter files with varied physical configurations,
runs both C++ and Rust lroche on each, and compares flux output.
Fails if any model exceeds 1e-5 max relative error.

Usage:
    cd /fred/oz480/mcoughli/lcurve/lcurve-rs
    module load gcccore/12.3.0 pgplot/5.2.2
    python3 test_data/run_regression100.py
"""

import os
import sys
import subprocess
import tempfile
import shutil
import random
import math

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
CPP_LROCHE = "/fred/oz480/mcoughli/lcurve/local/bin/lcurve/lroche"
RUST_LROCHE = os.path.join(PROJECT_DIR, "target", "release", "lroche")
CARGO_HOME = os.path.join(PROJECT_DIR, ".cargo_home")

# Thresholds
MAX_REL_ERROR_THRESHOLD = 1e-5
NTIME = 50
TIME1 = -0.2
TIME2 = 1.2

# Seed for reproducibility
random.seed(42)


def make_param_line(name, value, rng=0, dstep=0, vary=0):
    """Format a parameter line: name = value range dstep vary"""
    return f"{name:20s} = {value} {rng} {dstep} {vary}"


def make_flag_line(name, value):
    """Format a computational/flag parameter line: name = value"""
    return f"{name:20s} = {value}"


def base_model():
    """Return a dict of default parameter values matching test_model.dat."""
    return {
        # Physical parameters (value, range, dstep, vary)
        "q":              (0.5, 0.1, 0.01, 0),
        "iangle":         (82.0, 2.0, 0.1, 0),
        "r1":             (0.015, 0.005, 0.001, 0),
        "r2":             (-1, 0.01, 0.001, 0),
        "cphi3":          (0.015, 0.001, 0.001, 0),
        "cphi4":          (0.017, 0.001, 0.001, 0),
        "t1":             (15000, 500, 100, 0),
        "t2":             (3500, 200, 50, 0),
        "spin1":          (1, 0.001, 0.001, 0),
        "spin2":          (1, 0.001, 0.001, 0),
        "ldc1_1":         (0.4, 0.01, 0.01, 0),
        "ldc1_2":         (0.0, 0.01, 0.01, 0),
        "ldc1_3":         (0.0, 0.01, 0.01, 0),
        "ldc1_4":         (0.0, 0.01, 0.01, 0),
        "ldc2_1":         (0.6, 0.01, 0.01, 0),
        "ldc2_2":         (0.0, 0.01, 0.01, 0),
        "ldc2_3":         (0.0, 0.01, 0.01, 0),
        "ldc2_4":         (0.0, 0.01, 0.01, 0),
        "velocity_scale": (0, 1, 1, 0),
        "beam_factor1":   (0, 0.1, 0.02, 0),
        "beam_factor2":   (0, 0.1, 0.002, 0),
        "deltat":         (0, 0.001, 0.001, 0),
        "t0":             (0.0, 0.0001, 1e-05, 0),
        "period":         (1.0, 1e-06, 1e-06, 0),
        "gravity_dark1":  (0.25, 0.0001, 0.0001, 0),
        "gravity_dark2":  (0.08, 0.0001, 0.0001, 0),
        "absorb":         (0.5, 0.001, 0.001, 0),
        "slope":          (0, 0.01, 1e-05, 0),
        "quad":           (0, 0.01, 1e-05, 0),
        "cube":           (0, 0.01, 1e-05, 0),
        "third":          (0, 0.01, 1e-05, 0),
        "rdisc1":         (0.05, 0.01, 0.0001, 0),
        "rdisc2":         (0.35, 0.01, 0.02, 0),
        "height_disc":    (0.1, 0.01, 1e-05, 0),
        "beta_disc":      (1.5, 0.01, 1e-05, 0),
        "temp_disc":      (3000, 50, 40, 0),
        "texp_disc":      (-1.8, 0.2, 0.001, 0),
        "lin_limb_disc":  (0.3, 0.02, 0.0001, 0),
        "quad_limb_disc": (0, 0.02, 0.0001, 0),
        "radius_spot":    (0.35, 0.005, 0.01, 0),
        "length_spot":    (0.02, 0.002, 0.005, 0),
        "height_spot":    (0.05, 0.005, 1e-05, 0),
        "expon_spot":     (0, 0.2, 0.1, 0),
        "epow_spot":      (1, 0.1, 0.1, 0),
        "angle_spot":     (140, 2, 2, 0),
        "yaw_spot":       (0, 2, 2, 0),
        "temp_spot":      (15000, 400, 200, 0),
        "tilt_spot":      (90, 5, 2, 0),
        "cfrac_spot":     (0.2, 0.02, 0.008, 0),
    }


def base_flags():
    """Return a dict of default flag/computational parameter values."""
    return {
        "delta_phase":  1e-07,
        "nlat1f":       20,
        "nlat2f":       40,
        "nlat1c":       10,
        "nlat2c":       20,
        "npole":        0,
        "nlatfill":     0,
        "nlngfill":     0,
        "lfudge":       0.05,
        "llo":          0,
        "lhi":          -50,
        "phase1":       0.05,
        "phase2":       0.45,
        "wavelength":   550,
        "roche1":       0,
        "roche2":       1,
        "eclipse1":     1,
        "eclipse2":     1,
        "glens1":       0,
        "tperiod":      0.15,
        "gdark_bolom1": 1,
        "gdark_bolom2": 1,
        "mucrit1":      0,
        "mucrit2":      0,
        "limb1":        "Poly",
        "limb2":        "Poly",
        "use_radii":    1,
        "mirror":       0,
        "add_disc":     0,
        "nrad":         40,
        "opaque":       0,
        "add_spot":     0,
        "nspot":        100,
        "iscale":       0,
    }


# Parameter ordering for the model file
PARAM_ORDER = [
    "q", "iangle", "r1", "r2", "cphi3", "cphi4",
    "t1", "t2", "spin1", "spin2",
    "ldc1_1", "ldc1_2", "ldc1_3", "ldc1_4",
    "ldc2_1", "ldc2_2", "ldc2_3", "ldc2_4",
    "velocity_scale", "beam_factor1", "beam_factor2",
    "deltat", "t0", "period",
    "gravity_dark1", "gravity_dark2", "absorb",
    "slope", "quad", "cube", "third",
    "rdisc1", "rdisc2", "height_disc", "beta_disc",
    "temp_disc", "texp_disc", "lin_limb_disc", "quad_limb_disc",
    "radius_spot", "length_spot", "height_spot",
    "expon_spot", "epow_spot", "angle_spot", "yaw_spot",
    "temp_spot", "tilt_spot", "cfrac_spot",
]

FLAG_ORDER = [
    "delta_phase", "nlat1f", "nlat2f", "nlat1c", "nlat2c",
    "npole", "nlatfill", "nlngfill", "lfudge",
    "llo", "lhi", "phase1", "phase2", "wavelength",
    "roche1", "roche2", "eclipse1", "eclipse2", "glens1",
    "tperiod", "gdark_bolom1", "gdark_bolom2",
    "mucrit1", "mucrit2", "limb1", "limb2",
    "use_radii", "mirror", "add_disc", "nrad", "opaque",
    "add_spot", "nspot", "iscale",
]


def write_model_file(path, params, flags, comment=""):
    """Write a model parameter file."""
    with open(path, "w") as f:
        if comment:
            f.write(f"# {comment}\n")
        for name in PARAM_ORDER:
            val, rng, dstep, vary = params[name]
            f.write(make_param_line(name, val, rng, dstep, vary) + "\n")
        for name in FLAG_ORDER:
            f.write(make_flag_line(name, flags[name]) + "\n")


def uniform(lo, hi):
    return random.uniform(lo, hi)


def choice(options):
    return random.choice(options)


def eggleton_rl(q_donor):
    """Eggleton (1983) Roche lobe radius in units of binary separation.

    q_donor = mass of star whose Roche lobe is computed / mass of companion.
    """
    q23 = q_donor ** (2.0 / 3.0)
    q13 = q_donor ** (1.0 / 3.0)
    return 0.49 * q23 / (0.6 * q23 + math.log(1.0 + q13))


# ---------------------------------------------------------------------------
# Model generators — each returns (params, flags, description)
# ---------------------------------------------------------------------------

def gen_basic_stellar(idx):
    """Models 1-20: Basic stellar parameter variations."""
    p = base_model()
    f = base_flags()
    q = uniform(0.05, 1.0)
    iangle = uniform(60, 89)
    r1 = uniform(0.005, 0.05)
    t1 = uniform(8000, 30000)
    t2 = uniform(2500, 6000)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["r2"] = (-1, 0.01, 0.001, 0)  # Roche-lobe filling
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    # Vary absorb too
    absorb = uniform(0.0, 1.0)
    p["absorb"] = (round(absorb, 3), 0.001, 0.001, 0)
    desc = f"basic_stellar q={q:.3f} i={iangle:.1f} r1={r1:.4f}"
    return p, f, desc


def gen_spherical_secondary(idx):
    """Models 21-30: Spherical (non-Roche) secondary."""
    p = base_model()
    f = base_flags()
    q = uniform(0.1, 0.8)
    iangle = uniform(65, 88)
    r1 = uniform(0.01, 0.04)
    # Cap r2 at 90% of secondary Roche lobe to stay physical
    rl2 = eggleton_rl(q)
    r2 = uniform(0.1, min(0.4, 0.9 * rl2))
    t1 = uniform(10000, 25000)
    t2 = uniform(3000, 5500)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["r2"] = (round(r2, 4), 0.01, 0.001, 0)  # Positive = spherical
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    f["roche2"] = 0  # Spherical star 2
    desc = f"spherical_secondary r2={r2:.3f} roche2=0"
    return p, f, desc


def gen_eclipse_combos(idx):
    """Models 31-40: Eclipse on/off combinations."""
    p = base_model()
    f = base_flags()
    q = uniform(0.2, 0.9)
    iangle = uniform(70, 88)
    r1 = uniform(0.01, 0.03)
    t1 = uniform(10000, 20000)
    t2 = uniform(3000, 5000)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    # Cycle through eclipse combos
    combo = (idx - 31) % 4
    e1 = combo // 2
    e2 = combo % 2
    f["eclipse1"] = e1
    f["eclipse2"] = e2
    desc = f"eclipse_combo e1={e1} e2={e2}"
    return p, f, desc


def gen_gravity_darkening(idx):
    """Models 41-50: Gravity darkening variations."""
    p = base_model()
    f = base_flags()
    q = uniform(0.2, 0.8)
    iangle = uniform(65, 87)
    r1 = uniform(0.01, 0.04)
    t1 = uniform(10000, 25000)
    t2 = uniform(3000, 5500)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    gd1 = uniform(0.08, 0.32)
    gd2 = uniform(0.08, 0.32)
    absorb = uniform(0.0, 1.0)
    p["gravity_dark1"] = (round(gd1, 4), 0.0001, 0.0001, 0)
    p["gravity_dark2"] = (round(gd2, 4), 0.0001, 0.0001, 0)
    p["absorb"] = (round(absorb, 3), 0.001, 0.001, 0)
    # Vary bolometric flag
    gb1 = choice([0, 1])
    gb2 = choice([0, 1])
    f["gdark_bolom1"] = gb1
    f["gdark_bolom2"] = gb2
    desc = f"gravity_dark gd1={gd1:.3f} gd2={gd2:.3f} bolom={gb1},{gb2}"
    return p, f, desc


def gen_limb_darkening(idx):
    """Models 51-60: Limb darkening coefficient variations."""
    p = base_model()
    f = base_flags()
    q = uniform(0.2, 0.8)
    iangle = uniform(70, 87)
    r1 = uniform(0.01, 0.035)
    t1 = uniform(10000, 22000)
    t2 = uniform(3000, 5000)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    # Vary LDC coefficients
    p["ldc1_1"] = (round(uniform(0.1, 0.8), 3), 0.01, 0.01, 0)
    p["ldc1_2"] = (round(uniform(-0.5, 0.5), 3), 0.01, 0.01, 0)
    p["ldc1_3"] = (round(uniform(-0.3, 0.3), 3), 0.01, 0.01, 0)
    p["ldc1_4"] = (round(uniform(-0.2, 0.2), 3), 0.01, 0.01, 0)
    p["ldc2_1"] = (round(uniform(0.2, 0.9), 3), 0.01, 0.01, 0)
    p["ldc2_2"] = (round(uniform(-0.5, 0.5), 3), 0.01, 0.01, 0)
    p["ldc2_3"] = (round(uniform(-0.3, 0.3), 3), 0.01, 0.01, 0)
    p["ldc2_4"] = (round(uniform(-0.2, 0.2), 3), 0.01, 0.01, 0)
    # Choose limb darkening law
    limb_choice = choice(["Poly", "Claret"])
    f["limb1"] = limb_choice
    f["limb2"] = choice(["Poly", "Claret"])
    desc = f"limb_dark limb1={f['limb1']} limb2={f['limb2']}"
    return p, f, desc


def gen_with_disc(idx):
    """Models 61-75: With accretion disc."""
    p = base_model()
    f = base_flags()
    q = uniform(0.1, 0.8)
    iangle = uniform(65, 88)
    r1 = uniform(0.008, 0.03)
    t1 = uniform(12000, 28000)
    t2 = uniform(2800, 5000)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    f["add_disc"] = 1
    # Cap rdisc2 within primary Roche lobe
    rl1 = eggleton_rl(1.0 / q)
    rdisc1 = r1 + uniform(0.01, 0.03)  # inner > r1
    rdisc2 = uniform(0.2, min(0.45, 0.9 * rl1))
    p["rdisc1"] = (round(rdisc1, 5), 0.01, 0.0001, 0)
    p["rdisc2"] = (round(rdisc2, 4), 0.01, 0.02, 0)
    p["height_disc"] = (round(uniform(0.02, 0.15), 4), 0.01, 1e-05, 0)
    p["beta_disc"] = (round(uniform(1.0, 2.0), 3), 0.01, 1e-05, 0)
    p["temp_disc"] = (round(uniform(2000, 8000), 0), 50, 40, 0)
    p["texp_disc"] = (round(uniform(-2.5, -0.5), 3), 0.2, 0.001, 0)
    p["lin_limb_disc"] = (round(uniform(0.1, 0.6), 3), 0.02, 0.0001, 0)
    p["quad_limb_disc"] = (round(uniform(-0.2, 0.3), 3), 0.02, 0.0001, 0)
    # Vary opaque for some
    f["opaque"] = choice([0, 1])
    desc = f"disc rdisc2={rdisc2:.3f} opaque={f['opaque']}"
    return p, f, desc


def gen_with_spot(idx):
    """Models 76-85: With bright spot (no disc)."""
    p = base_model()
    f = base_flags()
    q = uniform(0.1, 0.8)
    iangle = uniform(65, 88)
    r1 = uniform(0.01, 0.03)
    t1 = uniform(12000, 25000)
    t2 = uniform(3000, 5000)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    f["add_disc"] = 1  # Spot requires disc to be enabled
    f["add_spot"] = 1
    # Disc params (need valid disc for spot)
    rl1 = eggleton_rl(1.0 / q)
    rdisc1 = r1 + uniform(0.01, 0.02)
    rdisc2 = uniform(0.25, min(0.45, 0.9 * rl1))
    p["rdisc1"] = (round(rdisc1, 5), 0.01, 0.0001, 0)
    p["rdisc2"] = (round(rdisc2, 4), 0.01, 0.02, 0)
    p["height_disc"] = (round(uniform(0.03, 0.12), 4), 0.01, 1e-05, 0)
    p["beta_disc"] = (round(uniform(1.0, 2.0), 3), 0.01, 1e-05, 0)
    p["temp_disc"] = (round(uniform(2500, 6000), 0), 50, 40, 0)
    p["texp_disc"] = (round(uniform(-2.0, -0.8), 3), 0.2, 0.001, 0)
    # Spot params
    p["radius_spot"] = (round(rdisc2 * uniform(0.8, 1.0), 4), 0.005, 0.01, 0)
    p["length_spot"] = (round(uniform(0.01, 0.05), 4), 0.002, 0.005, 0)
    p["height_spot"] = (round(uniform(0.02, 0.1), 4), 0.005, 1e-05, 0)
    p["expon_spot"] = (round(uniform(0, 2), 3), 0.2, 0.1, 0)
    p["epow_spot"] = (round(uniform(0.5, 2.0), 3), 0.1, 0.1, 0)
    p["angle_spot"] = (round(uniform(100, 180), 1), 2, 2, 0)
    p["yaw_spot"] = (round(uniform(-20, 20), 1), 2, 2, 0)
    p["temp_spot"] = (round(uniform(8000, 25000), 0), 400, 200, 0)
    p["tilt_spot"] = (round(uniform(60, 120), 1), 5, 2, 0)
    p["cfrac_spot"] = (round(uniform(0.05, 0.5), 3), 0.02, 0.008, 0)
    desc = f"spot angle={p['angle_spot'][0]:.0f} temp={p['temp_spot'][0]:.0f}"
    return p, f, desc


def gen_disc_and_spot(idx):
    """Models 86-95: Both disc and bright spot with varied params."""
    p = base_model()
    f = base_flags()
    q = uniform(0.15, 0.9)
    iangle = uniform(65, 88)
    r1 = uniform(0.008, 0.03)
    t1 = uniform(10000, 28000)
    t2 = uniform(2800, 5500)
    p["q"] = (round(q, 4), 0.1, 0.01, 0)
    p["iangle"] = (round(iangle, 2), 2.0, 0.1, 0)
    p["r1"] = (round(r1, 5), 0.005, 0.001, 0)
    p["t1"] = (round(t1, 0), 500, 100, 0)
    p["t2"] = (round(t2, 0), 200, 50, 0)
    f["add_disc"] = 1
    f["add_spot"] = 1
    # Disc params — cap within primary Roche lobe
    rl1 = eggleton_rl(1.0 / q)
    rdisc1 = r1 + uniform(0.01, 0.025)
    rdisc2 = uniform(0.2, min(0.45, 0.9 * rl1))
    p["rdisc1"] = (round(rdisc1, 5), 0.01, 0.0001, 0)
    p["rdisc2"] = (round(rdisc2, 4), 0.01, 0.02, 0)
    p["height_disc"] = (round(uniform(0.02, 0.15), 4), 0.01, 1e-05, 0)
    p["beta_disc"] = (round(uniform(1.0, 2.0), 3), 0.01, 1e-05, 0)
    p["temp_disc"] = (round(uniform(2000, 7000), 0), 50, 40, 0)
    p["texp_disc"] = (round(uniform(-2.5, -0.5), 3), 0.2, 0.001, 0)
    p["lin_limb_disc"] = (round(uniform(0.1, 0.5), 3), 0.02, 0.0001, 0)
    p["quad_limb_disc"] = (round(uniform(-0.1, 0.3), 3), 0.02, 0.0001, 0)
    f["opaque"] = choice([0, 1])
    # Spot params
    p["radius_spot"] = (round(rdisc2 * uniform(0.85, 1.0), 4), 0.005, 0.01, 0)
    p["length_spot"] = (round(uniform(0.01, 0.04), 4), 0.002, 0.005, 0)
    p["height_spot"] = (round(uniform(0.02, 0.08), 4), 0.005, 1e-05, 0)
    p["expon_spot"] = (round(uniform(0, 1.5), 3), 0.2, 0.1, 0)
    p["epow_spot"] = (round(uniform(0.5, 2.0), 3), 0.1, 0.1, 0)
    p["angle_spot"] = (round(uniform(110, 170), 1), 2, 2, 0)
    p["yaw_spot"] = (round(uniform(-15, 15), 1), 2, 2, 0)
    p["temp_spot"] = (round(uniform(8000, 22000), 0), 400, 200, 0)
    p["tilt_spot"] = (round(uniform(70, 110), 1), 5, 2, 0)
    p["cfrac_spot"] = (round(uniform(0.05, 0.4), 3), 0.02, 0.008, 0)
    # Also vary gravity/limb for coverage
    gd1 = uniform(0.08, 0.32)
    gd2 = uniform(0.08, 0.32)
    p["gravity_dark1"] = (round(gd1, 4), 0.0001, 0.0001, 0)
    p["gravity_dark2"] = (round(gd2, 4), 0.0001, 0.0001, 0)
    desc = f"disc+spot opaque={f['opaque']}"
    return p, f, desc


def gen_edge_cases(idx):
    """Models 96-100: Edge cases — low inclination, high q, t2<0, mirror, roche1."""
    p = base_model()
    f = base_flags()
    case = idx - 96
    if case == 0:
        # Low inclination — no eclipses expected
        p["iangle"] = (30.0, 2.0, 0.1, 0)
        p["q"] = (0.5, 0.1, 0.01, 0)
        desc = "edge_low_incl i=30"
    elif case == 1:
        # High mass ratio
        p["q"] = (1.0, 0.1, 0.01, 0)
        p["iangle"] = (80.0, 2.0, 0.1, 0)
        p["r1"] = (0.01, 0.005, 0.001, 0)
        desc = "edge_high_q q=1.0"
    elif case == 2:
        # Negative t2 (means t2 = |t2| * t1)
        p["t2"] = (-0.3, 200, 50, 0)
        p["iangle"] = (78.0, 2.0, 0.1, 0)
        desc = "edge_neg_t2 t2=-0.3"
    elif case == 3:
        # Mirror symmetry mode
        p["q"] = (0.4, 0.1, 0.01, 0)
        p["iangle"] = (82.0, 2.0, 0.1, 0)
        f["mirror"] = 1
        desc = "edge_mirror mirror=1"
    elif case == 4:
        # Roche star 1
        p["q"] = (0.6, 0.1, 0.01, 0)
        p["iangle"] = (80.0, 2.0, 0.1, 0)
        p["r1"] = (0.02, 0.005, 0.001, 0)
        f["roche1"] = 1
        desc = "edge_roche1 roche1=1"
    else:
        desc = "edge_default"
    return p, f, desc


# ---------------------------------------------------------------------------
# Model generation dispatcher
# ---------------------------------------------------------------------------

def generate_model(idx):
    """Generate model parameters for model number idx (1-100)."""
    if 1 <= idx <= 20:
        return gen_basic_stellar(idx)
    elif 21 <= idx <= 30:
        return gen_spherical_secondary(idx)
    elif 31 <= idx <= 40:
        return gen_eclipse_combos(idx)
    elif 41 <= idx <= 50:
        return gen_gravity_darkening(idx)
    elif 51 <= idx <= 60:
        return gen_limb_darkening(idx)
    elif 61 <= idx <= 75:
        return gen_with_disc(idx)
    elif 76 <= idx <= 85:
        return gen_with_spot(idx)
    elif 86 <= idx <= 95:
        return gen_disc_and_spot(idx)
    elif 96 <= idx <= 100:
        return gen_edge_cases(idx)
    else:
        raise ValueError(f"Model index {idx} out of range")


# ---------------------------------------------------------------------------
# Running lroche and parsing output
# ---------------------------------------------------------------------------

def run_cpp_lroche(model_path, output_path):
    """Run C++ lroche and return (success, error_msg)."""
    cmd = (
        f'echo "1" | {CPP_LROCHE} {model_path} none '
        f'{TIME1} {TIME2} {NTIME} 0.001 1 0.0 12345 1 {output_path} none'
    )
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0:
        return False, f"C++ exit code {result.returncode}: {result.stderr[:200]}"
    if not os.path.exists(output_path):
        return False, "C++ produced no output file"
    return True, ""


def run_rust_lroche(model_path, output_path):
    """Run Rust lroche and return (success, error_msg)."""
    env = os.environ.copy()
    env["CARGO_HOME"] = CARGO_HOME
    cmd = [
        RUST_LROCHE, model_path, "none",
        f"--time1={TIME1}", f"--time2={TIME2}", f"--ntime={NTIME}",
        f"--output={output_path}",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=120, env=env,
    )
    if result.returncode != 0:
        return False, f"Rust exit code {result.returncode}: {result.stderr[:200]}"
    if not os.path.exists(output_path):
        return False, "Rust produced no output file"
    return True, ""


def parse_lightcurve(path):
    """Parse lroche output file. Returns list of (time, flux) tuples."""
    points = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                time = float(parts[0])
                flux = float(parts[2])
                points.append((time, flux))
    return points


def compare_outputs(cpp_path, rust_path):
    """Compare two lightcurve files. Returns (max_rel_err, mean_rel_err, n_points)."""
    cpp_data = parse_lightcurve(cpp_path)
    rust_data = parse_lightcurve(rust_path)

    if len(cpp_data) != len(rust_data):
        return None, None, 0

    max_rel = 0.0
    sum_rel = 0.0
    n = len(cpp_data)

    for (tc, fc), (tr, fr) in zip(cpp_data, rust_data):
        # Check times match
        if abs(tc - tr) > 1e-12:
            return None, None, 0
        # Relative error on flux
        denom = max(abs(fc), abs(fr), 1e-30)
        rel = abs(fc - fr) / denom
        max_rel = max(max_rel, rel)
        sum_rel += rel

    mean_rel = sum_rel / n if n > 0 else 0.0
    return max_rel, mean_rel, n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Create temp directory for all model/output files
    tmpdir = tempfile.mkdtemp(prefix="lroche_regtest_")
    print(f"Working directory: {tmpdir}")
    print(f"C++ binary: {CPP_LROCHE}")
    print(f"Rust binary: {RUST_LROCHE}")
    print()

    results = []  # (idx, desc, max_rel, mean_rel, n_points, status, msg)
    n_pass = 0
    n_fail = 0
    n_error = 0
    n_skip = 0

    for idx in range(1, 101):
        params, flags, desc = generate_model(idx)

        model_path = os.path.join(tmpdir, f"model_{idx:03d}.dat")
        cpp_out = os.path.join(tmpdir, f"cpp_{idx:03d}.dat")
        rust_out = os.path.join(tmpdir, f"rust_{idx:03d}.dat")

        write_model_file(model_path, params, flags, comment=f"Model {idx}: {desc}")

        # Run C++
        ok_cpp, err_cpp = run_cpp_lroche(model_path, cpp_out)
        if not ok_cpp:
            # Check if Rust also fails — if so it's a bad parameter combo, skip.
            # If Rust succeeds but C++ doesn't, also skip (C++ limitation).
            ok_rust, err_rust = run_rust_lroche(model_path, rust_out)
            if not ok_rust:
                print(f"  [{idx:3d}] SKIP  both fail: {err_cpp[:80]}")
                results.append((idx, desc, None, None, 0, "SKIP", f"Both: {err_cpp}"))
            else:
                print(f"  [{idx:3d}] SKIP  C++ fails (Rust ok): {err_cpp[:80]}")
                results.append((idx, desc, None, None, 0, "SKIP", f"C++ only: {err_cpp}"))
            n_skip += 1
            continue

        # Run Rust
        ok_rust, err_rust = run_rust_lroche(model_path, rust_out)
        if not ok_rust:
            print(f"  [{idx:3d}] ERROR (Rust): {err_rust}")
            results.append((idx, desc, None, None, 0, "ERROR", f"Rust: {err_rust}"))
            n_error += 1
            continue

        # Compare
        max_rel, mean_rel, n_points = compare_outputs(cpp_out, rust_out)
        if max_rel is None:
            msg = "Output mismatch (different lengths or times)"
            print(f"  [{idx:3d}] ERROR: {msg}")
            results.append((idx, desc, None, None, n_points, "ERROR", msg))
            n_error += 1
            continue

        if max_rel > MAX_REL_ERROR_THRESHOLD:
            status = "FAIL"
            n_fail += 1
            marker = "FAIL"
        else:
            status = "PASS"
            n_pass += 1
            marker = "pass"

        print(f"  [{idx:3d}] {marker}  max_rel={max_rel:.2e}  mean_rel={mean_rel:.2e}  {desc}")
        results.append((idx, desc, max_rel, mean_rel, n_points, status, ""))

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    print()
    print("=" * 90)
    print(f"{'SUMMARY':^90s}")
    print("=" * 90)
    print(f"{'#':>4s}  {'Status':>6s}  {'Max Rel':>10s}  {'Mean Rel':>10s}  {'Pts':>4s}  Description")
    print("-" * 90)
    for idx, desc, max_rel, mean_rel, n_points, status, msg in results:
        if max_rel is not None:
            print(f"{idx:4d}  {status:>6s}  {max_rel:10.2e}  {mean_rel:10.2e}  {n_points:4d}  {desc}")
        else:
            print(f"{idx:4d}  {status:>6s}  {'N/A':>10s}  {'N/A':>10s}  {n_points:4d}  {desc} -- {msg}")
    print("-" * 90)
    print(f"  PASS: {n_pass}   FAIL: {n_fail}   SKIP: {n_skip}   ERROR: {n_error}   TOTAL: {len(results)}")
    print("=" * 90)

    # Find worst-case
    valid = [(r[2], r[0], r[1]) for r in results if r[2] is not None]
    if valid:
        worst_rel, worst_idx, worst_desc = max(valid)
        print(f"  Worst case: model {worst_idx} ({worst_desc}), max_rel_error = {worst_rel:.2e}")

    # Clean up
    shutil.rmtree(tmpdir)
    print(f"\nCleaned up {tmpdir}")

    if n_fail > 0 or n_error > 0:
        print(f"\nFAILED: {n_fail} models exceeded threshold, {n_error} Rust errors")
        sys.exit(1)
    else:
        msg = f"\nALL {n_pass} MODELS PASSED (max relative error < {MAX_REL_ERROR_THRESHOLD})"
        if n_skip > 0:
            msg += f", {n_skip} skipped (C++ failures)"
        print(msg)
        sys.exit(0)


if __name__ == "__main__":
    main()
