# Installation

## Prerequisites

- [Rust toolchain](https://rustup.rs/) (stable)

For Python bindings:

- Python >= 3.9
- [maturin](https://github.com/PyO3/maturin)
- numpy

## CLI binary

```bash
cd lcurve-rs
cargo build --release
```

The binary is at `target/release/lroche`. Add it to your `PATH` or copy it to a local `bin/` directory.

## Python bindings

```bash
pip install maturin numpy
cd lcurve-rs/python
maturin develop --release
```

This builds and installs the `lcurve_rs` package into the active virtualenv or conda environment.

!!! tip "Virtual environment"
    Maturin requires a virtualenv or conda environment. Create one with:
    ```bash
    python -m venv .venv
    source .venv/bin/activate
    ```

## Building a wheel

To build a distributable wheel:

```bash
cd lcurve-rs/python
maturin build --release
# Wheel at target/wheels/lcurve_rs-*.whl
pip install target/wheels/lcurve_rs-*.whl
```

## Controlling parallelism

lcurve-rs uses [Rayon](https://github.com/rayon-rs/rayon) for multi-core parallelism. By default, it uses all available cores. To limit thread count:

```bash
export RAYON_NUM_THREADS=4
```
