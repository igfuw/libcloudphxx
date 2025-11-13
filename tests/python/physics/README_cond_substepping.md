# Condensation Substepping Tests

This directory contains tests for condensation with various substepping configurations.

## Files

- `lgrngn_cond_substepping.py` - Main test script that runs condensation simulations
- `lgrngn_cond_substepping_test.py` - Comparison script that validates results against reference data
- `lgrngn_cond_substepping_plot.py` - Plotting script to visualize results

## Workflow

### 1. Run the test and generate results

```bash
python lgrngn_cond_substepping.py
```

This will:
- Run condensation tests with various configurations (mixing, const_p, exact substepping, adaptive, etc.)
- Test different RH formulas (pv_cc, rv_cc, pv_tet, rv_tet)
- Test different numbers of substeps
- Save results to `test_results/lgrngn_cond_substepping_results.csv`

### 2. Create reference data (first time only)

After verifying that the results are correct, save them as reference data:

```bash
python lgrngn_cond_substepping.py --save-ref
```

This creates `test_results/lgrngn_cond_substepping_refdata.csv` which will be used for future comparisons.

### 3. Compare results with reference data

On subsequent runs, compare new results against the reference:

```bash
python lgrngn_cond_substepping_test.py
```

This will:
- Load both the results and reference data
- Compare each variable with appropriate tolerances
- Report any differences that exceed tolerances
- Exit with code 0 if all tests pass, 1 if any fail

### 4. Visualize results

Generate plots of the results:

```bash
python lgrngn_cond_substepping_plot.py
```

This creates plots in `test_results/plots/` showing:
- Supersaturation vs. number of substeps
- Theta and rv leakage vs. number of substeps
- Activated concentration vs. number of substeps
- Radius moments vs. number of substeps
- Execution time vs. number of substeps

Each variable is plotted for all RH formulas with logarithmic x-axis.

## Test Tolerances

The comparison script uses the following default tolerances:

- `ss` (supersaturation): 1% relative tolerance
- `th_diff`, `rv_diff` (leakage): Absolute tolerances (1e-5, 1e-6)
- `act`, `mr`, `sr`, `tr` (concentration and moments): 2% relative tolerance
- `exectime`: Not compared (machine-dependent)

## Configuration Options Tested

- **mixing**: Communication of rv/th changes between SDs after each substep
- **constp**: Constant pressure vs. variable pressure
- **exact_sstp**: Per-particle vs. per-cell substepping
- **adaptive**: Adaptive determination of substep count
- **sstp_cond**: Number of condensation substeps (1, 2, 3, 4, 6, 8, 32)
- **sstp_cond_act**: Number of substeps for activating droplets (1, 8)
- **RH_formula**: Relative humidity formula (pv_cc, rv_cc, pv_tet, rv_tet)

## Example Complete Workflow

```bash
# Initial run - establish reference data
python lgrngn_cond_substepping.py --save-ref

# Generate initial plots
python lgrngn_cond_substepping_plot.py

# After making code changes, run test again
python lgrngn_cond_substepping.py

# Compare with reference
python lgrngn_cond_substepping_test.py

# Generate new plots to compare visually
python lgrngn_cond_substepping_plot.py
```

## Updating Reference Data

If intentional changes are made that affect results:

1. Verify the new results are correct
2. Update reference data: `python lgrngn_cond_substepping.py --save-ref`
3. Document the reason for the change in git commit message
