# Reproducibility, IDAES Modernization, and Private Reactor Development Plan

## Objective

Reproduce every reported result from Ghosh et al. (2024), modernize the
flowsheet from IDAES 2.0.0.dev3 to the current stable IDAES 2.12.0, and prepare
an isolated private development path for a null-space Gaussian-process reactor.
All numerical comparisons will distinguish exact reproduction, agreement within
a stated tolerance, and explained non-reproduction.

## Repository and branch topology

### Public work

1. `main`
   - Published repository snapshot.
2. `codex/reproduce-published-results`
   - Branched from `main`.
   - Recovers and pins the historical software environment.
   - Adds scripted, testable workflows for all paper cases.
   - Opens PR A into `main`.
3. `codex/idaes-2.12-modernization`
   - Branched from the head of `codex/reproduce-published-results` so that it
     directly builds on the reproducibility infrastructure.
   - Migrates and validates the complete workflow under IDAES 2.12.0.
   - Opens PR B into `main`; it remains unmerged until reproduction is judged
     reliable.

### Private work

Use a separate private clone/repository rather than a private branch in the
public clone:

- `origin`: private repository, with push access.
- `upstream`: `https://github.com/dowlinglab/cistar-process-modeling.git`, used
  for fetching public history.
- Private development branch: `codex/null-space-reactor-development`, starting
  from the accepted Phase B commit.
- The private surrogate-training repository remains a separately versioned
  dependency. Do not use a public-facing private Git submodule.
- Never configure the public repository as a push URL in the private working
  clone.

When the reactor integration is ready for release, create a fresh branch from
the then-current public `main` and cherry-pick or squash only audited,
publication-ready commits. Verify that no private data, model artifacts,
credentials, private URLs, or unwanted development history are present before
opening the public PR.

## Phase A - Historical reproduction

### Environment recovery

1. Inventory version evidence from the paper, Supporting Information, notebook
   metadata, serialized outputs, and IDAES release history.
2. Attempt the exact published stack first:
   - IDAES 2.0.0.dev3
   - Pyomo 6.4.2
   - Python 3.8.5
   - IPOPT 3.13.2 with MA27
3. If the development build cannot be recovered, identify its source commit and
   test the closest stable release after checking whether project-specific
   contributions were required.
4. Add a portable environment specification, a fully resolved lock or explicit
   package manifest, solver provenance, and one-command setup instructions.
5. Also test MA57, recording solver options, status, iterations, wall time, and
   numerical differences.

### Reproducibility harness

1. Extract notebook orchestration into importable modules and thin CLI entry
   points without changing equations.
2. Preserve the notebooks as documented analyses while eliminating hidden
   execution-order dependencies.
3. Record model statistics, degrees of freedom, scaling state, initialization
   status, solver termination, constraint residuals, and objective components.
4. Add smoke tests for model construction and each unit initialization stage.
5. Add regression tables that compare regenerated values with the checked-in
   reference results.

### Complete result matrix

Reproduce, at minimum:

- ROK models M2, M3, M4, and M5 for the Bakken base case.
- All reported carbon-tax cases for M5.
- Bakken, the Eagle Ford aggregate, and EF-1 through EF-12.
- Initialized and optimized composite curves.
- Stream and heat-integration tables.
- Product flow and composition plots.
- MSP, TAC, hydrogen rebate, LHV, upstream/downstream emissions, carbon
  efficiency, heating utility, cooling utility, recycle behavior, and all
  reported operating decisions.

### Phase A gate

- Every published table/figure datum has a machine-readable comparison result.
- Tolerances are declared before judging the result.
- Every failure is classified as environment, solver, nondeterminism,
  implementation defect, missing provenance, or irreproducible result.
- Setup and execution succeed from a clean environment using documented
  commands.
- PR A is pushed and opened against public `main`.

## Phase B - IDAES 2.12.0 migration

1. Branch from the Phase A head.
2. Build a new isolated environment pinned to IDAES 2.12.0 and compatible
   Python, Pyomo, and numerical dependencies.
3. Migrate public APIs before changing scientific equations.
4. Validate in increasing scope:
   - property and reaction packages;
   - each unit operation independently;
   - translators and phase-boundary behavior;
   - open-loop unit sequence;
   - recycle closure;
   - costing and emissions;
   - embedded heat integration;
   - optimization.
5. Use current IDAES initializers where appropriate, while retaining a legacy
   initialization path until numerical parity is demonstrated.
6. Use the IDAES Diagnostics Toolbox and scaling tools to record structural
   singularities, poorly scaled variables/constraints, near-parallel
   constraints, extreme Jacobian entries, and residuals.
7. Compare MA27 and MA57 using identical starting points and tolerances.
8. Run the complete Phase A regression matrix under IDAES 2.12.0.

### Phase B gate

- All cases converge reliably from documented initialization states, or each
  exception has a minimized reproducer and evidence-backed explanation.
- Scientific outputs satisfy the declared parity tolerances or have a traced
  cause for the difference.
- No new result silently overwrites the historical reference data.
- PR B is pushed and opened against public `main` and held for review.

## Phase C - Private null-space GP reactor development

Phase C begins only in the isolated private clone.

### Version 1 scope

- Legacy10 data.
- Independent GPyTorch Matérn-3/2 GP using the accepted rotated basis as the
  scientific reference.
- scikit-learn parity implementation.
- PySMO prototype/embedding benchmark, without requiring it to reproduce a
  different kernel or training algorithm exactly.
- GP predictive mean only.
- Process-scale formulation based on intensive conversion/product-carbon
  quantities and space velocity, not direct multiplication of laboratory outlet
  flows.
- Propylene-rich reactive feed with inert/pass-through handling.
- Explicit C10-C18 pseudo-components.
- Optional convex-hull/trust-region domain constraint and extrapolation
  diagnostics.
- Modular nonnegativity handling, beginning with a constrained correction
  budget and leaving clean extension points for lexicographic and objective-
  penalty formulations.
- Catalyst mass may differ substantially from the laboratory value but remains
  coupled to throughput by the space-velocity definition and training-domain
  limits.

### Phase C validation gates

1. Standalone numerical parity with GPyTorch predictions.
2. Exact elemental conservation to numerical tolerance.
3. Nonnegative outlet flows within the configured correction budget.
4. Derivative checks for every algebraic GP expression.
5. Standalone IDAES unit initialization over nominal and boundary cases.
6. Integration into the modernized flowsheet only after standalone tests pass.
7. Mean-only process optimization before any uncertainty propagation.
8. LOO refit/reoptimization as the first uncertainty extension after this
   96-hour goal.

## 96-hour execution priorities

1. Phase A environment recovery and one end-to-end reference case.
2. Phase A full published result matrix and PR A.
3. Phase B model construction and unit-by-unit migration.
4. Phase B full matrix, discrepancy analysis, and PR B.
5. Establish the private Phase C repository topology and write the reactor
   interface/specification; implement only after the private destination exists.

The manuscript repositories and manuscript text remain out of scope until the
computational results stabilize.
