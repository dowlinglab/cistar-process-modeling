# Historical environment evidence

## Published stack

The main article reports:

- IDAES-PSE `2.0.0.dev3`;
- Pyomo `6.4.2`;
- IPOPT `3.13.2` with HSL MA27;
- Lenovo ThinkPad P14s Gen 1, 1.8 GHz quad-core Intel Core i7, 16 GB RAM; and
- 5,515 continuous variables and 5,507 equality constraints.

The paper reports approximate runtimes of one hour for sequential unit
initialization, less than five minutes for flowsheet convergence with costing,
and less than ten minutes for each optimization. It also reports up to a fourfold
initialization speedup with Pyomo 6.5.1, while retaining Pyomo 6.4.2 for the
published results.

## Evidence embedded in the repository

- All original flowsheet/convergence/optimization notebooks declare Python
  `3.8.5` in `metadata.language_info.version`.
- The three analysis notebooks were later executed with Python `3.8.16`.
- Saved warnings in all three analysis notebooks expose the original environment
  name and installation path:
  `idaes-200dev3-pyomo-642` under Python 3.8.
- Saved solver output repeatedly identifies IPOPT `3.13.2`, MA27, and
  `tol=1e-06`.
- The notebooks instantiate `SolverFactory('ipopt')` without recording an
  explicit executable path. Therefore, the exact IPOPT build and linked HSL
  library are not identified by committed source.

## IDAES Git-history constraint

The official IDAES source reported `2.0.0.dev3` across a long development
interval, from commit `6bdf3dbf` (2022-05-20) until the version advanced near
commit `c2fbb63c` (2023-05-31). The version string does not identify a unique
commit.

Dependency history makes a Pyomo 6.4.2-compatible source checkout most plausible
around August-November 2022. A May 2023 IDAES checkout could still print
`2.0.0.dev3` while expecting a newer/custom Pyomo around 6.5.1. The saved output
dates and paper text therefore do not uniquely resolve the source revision.

## Reconstruction policy

Phase A will call the selected stack a **defensible reconstruction**, not the
exact historical environment, unless additional evidence identifies the source
commit and solver build. Candidate evaluation will:

1. bound compatible IDAES commits using imports and API calls in the notebooks;
2. test the strongest Pyomo 6.4.2 candidate first;
3. test a slightly later stable/release candidate if a required development
   contribution is missing;
4. record exact source commits, package hashes, solver output, and platform; and
5. compare MA27 and MA57 from identical saved initial points.

The retired historical IDAES package channel tested during archaeology returned
404, so the reconstruction cannot rely on that channel remaining available.
