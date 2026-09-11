# Published-result reproducibility audit

This directory keeps the published targets separate from checked-in result
snapshots and future reruns.

## Files

- `published_reference.json` transcribes the quantitative values printed in the
  main article and Supporting Information. It also declares transcription
  tolerances equal to one half-unit in the final displayed digit.
- `result_manifest.json` maps every main-text and supporting figure/table to its
  notebook, archived data, and required comparison method.
- `audit_published_tables.py` compares Tables S4-S6 with the two checked-in CSV
  families.
- `runs/` contains immutable structured records for executed reproduction
  attempts, including failed or interrupted attempts that affect provenance.

The transcription tolerances are not solver regression tolerances. Phase A will
add tighter, metric-specific solver tolerances after repeat-run variability has
been measured with the reconstructed environment.

## Run the current audit

From the repository root:

```bash
python reproducibility/audit_published_tables.py
python reproducibility/audit_published_tables.py --snapshot migrated
python -m unittest discover -s tests -v
```

After creating a historical candidate environment, rerun the M5/Bakken tax
sequence without modifying archived checkpoints or CSVs:

```bash
python reproducibility/run_m5_bakken_tax_series.py \
  --ipopt /absolute/path/to/hsl-enabled/ipopt \
  --linear-solver ma27 \
  --tee \
  --output /path/outside/the/repository/m5-bakken-tax-series.json
```

Tax rates are expressed in USD/kg CO2e, matching the optimization notebook;
the corresponding reported table values are 1,000 times larger in USD/tonne.
To isolate a difficult sequential transition, `--initial-optimal-tax 0.19
--tax-rates 0.41` loads the archived USD 190/tonne optimum before solving the
USD 410/tonne case.

The default `postprocessed` snapshot uses the root-level CSV files created by
commit `957e363`, which recalculated TAC and MSP using cooling water above 303 K
and refrigerated water from 288-303 K. The `migrated` snapshot uses the CSV files
under `results/` that arrived in the initial public migration commit `2295e03`.
Neither snapshot is silently treated as the paper.

## Known provenance findings

The published process/downstream emissions in Tables S4-S6 and Figures 3 and 6
do not match the CSV values. For example, Table S4 reports 9.17 g CO2e/MJ for M2,
while both checked-in CSV snapshots report 8.876.

The public implementation constrains both upstream and downstream emissions
using this denominator in `src/emissions_calculations.py`:

```text
fuel_energy + gas_energy + h2_energy
```

The paper describes emissions per unit of liquid fuel, and the printed values
appear to use an earlier normalization convention. The public repository starts
with a bulk migration from a private research repository, so its Git history
does not contain the equation version that generated the paper numbers. The
audit labels this unresolved issue `EMISSIONS-NORMALIZATION-001` and continues
to fail it under `--strict`.

The audit also labels 22 postprocessed decision/cost values
`PUBLISHED-SNAPSHOT-DRIFT-001`. They exceed literal print-rounding tolerance,
which establishes that the SI tables and the checked-in higher-precision CSVs
were not generated from precisely the same saved optimization snapshot. The
differences are reported rather than erased by widening tolerances.

Finally, the migrated snapshot has 29 additional MSP/TAC differences labeled
`LEGACY-ECONOMICS-POSTPROCESSING-001`. Commit `957e363` explains these: it
recalculated economics using separate cooling-water and refrigerated-water
grades, producing the root-level `postprocessed` CSV family.

## Evidence sources

- Ghosh et al. (2024), DOI `10.1021/acssuschemeng.4c00933`.
- Main article Figures 1-8, especially the numeric labels in Figures 3, 6, and 7.
- Supporting Information Figures S1-S8 and Tables S1-S6.
- `results/solution_data.xlsx`, containing 23 stream sheets and 46 initialized/
  optimal heat-integration sheets plus a README sheet.
- Saved notebook metadata and IPOPT output.

The source PDFs are not redistributed in this repository.
