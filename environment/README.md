# Phase A environment reconstruction

The published software labels do not uniquely identify an installable
environment. This directory therefore distinguishes the target reconstruction
from tested platform-specific candidates.

## Target reconstruction

`phase-a-paper-reconstruction.yml` pins the versions recoverable from the paper,
notebook metadata, saved output, and the closest contemporary environment still
present on the development workstation. Create it with:

```bash
conda env create -f environment/phase-a-paper-reconstruction.yml
conda activate cistar-paper-reproduction
python -m pip install --no-deps \
  "git+https://github.com/IDAES/idaes-pse.git@66935c80a5aafc3ffc9ab3d387e488cddd4f233b"
```

IDAES commit `66935c80a5aafc3ffc9ab3d387e488cddd4f233b` is the first official
IDAES commit that pins the 6.4.2 release of the IDAES Pyomo fork. It reports
`2.0.0.dev3`, contains every IDAES import used by this repository, and predates
the subsequent Pyomo 6.4.3 transition. It is the strongest initial candidate,
not proof of Kanishka's exact private environment.

The environment file uses the public PyPI Pyomo 6.4.2 release for portability.
If model behavior differs, test the IDAES Pyomo fork from the URL embedded in
the historical IDAES `setup.py` and record that as a separate candidate.

## Candidate A1: Apple Silicon compatibility smoke test

Candidate A1 was created on 2026-09-11 by cloning a surviving 2022-era local
environment and replacing IDAES 1.13.0 with commit `66935c80` using
`pip install --no-deps`. Its relevant versions are:

- macOS 26.6.2 on arm64;
- Python 3.10.8 (platform deviation from notebook Python 3.8.5);
- IDAES 2.0.0.dev3 at commit `66935c80`;
- Pyomo 6.4.2;
- NumPy 1.23.5, pandas 1.5.1, SciPy 1.9.3, and matplotlib 3.6.2;
- bundled IPOPT 3.14.10 with MUMPS; and
- external IPOPT 3.14.19 with HSL MA27 and MA57.

All repository Python modules import under A1. The M5/Bakken flowsheet builds
with 24,610 component data objects and four degrees of freedom before
initialization. Small nonlinear solves terminate optimally with both MA27 and
MA57. A1 is suitable for compatibility diagnostics, but it does not satisfy the
paper stack gate because Python, platform, and IPOPT differ.

The first full M5/Bakken optimization rerun at zero carbon tax also terminated
optimally. Starting from the exact unit, constrained, and costed checkpoint
chain, the eight-degree-of-freedom problem converged in 20 IPOPT iterations and
reproduced every migrated-CSV quantity at or beyond its stored precision. See
`reproducibility/runs/2026-09-11-candidate-a1-m5-bakken-tax0.json`. This proves
one case is computationally reproducible under A1, while the paper-stack and
full-matrix gates remain open.

Run the repeatable smoke checks with:

```bash
conda run -n cistar-paper-reproduction-a1 \
  python reproducibility/smoke_environment.py \
  --ipopt /absolute/path/to/hsl-enabled/ipopt
```

The paper used IPOPT 3.13.2 with MA27. The locally retained IPOPT 3.13.2 binary
has MUMPS only, while the locally built IPOPT 3.14.19 has both MA27 and MA57.
Solver-version and linear-solver effects must therefore be reported separately.
