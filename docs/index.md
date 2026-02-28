# lcurve-rs

Fast Rust light curve engine for eclipsing white-dwarf binaries, with Python bindings.

A port of Tom Marsh's C++ [lcurve](http://github.com/trmrsh/cpp-lcurve/tree/master) and associated packages ([cpp-roche](https://github.com/trmrsh/cpp-roche), [cpp-subs](https://github.com/trmrsh/cpp-subs)), rewritten in Rust for memory safety, multi-core parallelism, and easy installation.

## Key Features

- **Roche geometry** — Roche-lobe-filling stars with gravity darkening and limb darkening (polynomial + Claret 4-coefficient)
- **Accretion disc** — opaque/transparent disc, power-law temperature profile, bright spot, disc edge
- **Star spots** — up to 3 Gaussian spots on the primary, 2 on the secondary, plus uniform equatorial spots
- **Irradiation** — mutual irradiation with configurable absorption efficiency
- **Parallel** — Rayon multithreading across time points and grid elements
- **Python bindings** — `lcurve_rs` package via PyO3/maturin with numpy arrays

## Quick Example

```python
import lcurve_rs
import numpy as np

model = lcurve_rs.Model("model.dat")
result = model.light_curve(time1=-0.2, time2=1.2, ntime=1000)

import matplotlib.pyplot as plt
plt.plot(result.times, result.flux)
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.show()
```

## Getting Started

See the [Installation](getting-started/installation.md) and [Quick Start](getting-started/quickstart.md) guides.
