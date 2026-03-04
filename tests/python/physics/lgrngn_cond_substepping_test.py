#!/usr/bin/env python3
"""
Test script to compare condensation substepping results against reference data.
This script runs the lgrngn_cond_substepping.py test and compares the results
with previously stored reference data.
"""

import sys
import os
import csv
import math


def _read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _to_bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("true", "1", "t", "yes", "y"):
        return True
    if s in ("false", "0", "f", "no", "n"):
        return False
    raise ValueError(f"Cannot parse boolean from: {v!r}")


def _to_number_or_str(v):
    """Best-effort parse to float; keep original string if not numeric."""
    if v is None:
        return ""
    s = str(v).strip()
    if s == "":
        return ""
    try:
        return float(s)
    except ValueError:
        return s


def _config_key(row, config_cols):
    key = []
    for c in config_cols:
        if c in ("mixing", "constp", "exact_sstp", "adaptive"):
            key.append(_to_bool(row[c]))
        elif c in ("RH_formula", "sstp_cond", "sstp_cond_act"):
            # Those are integers in practice, but allow float-in-CSV (e.g. "1.0")
            key.append(int(float(row[c])))
        else:
            key.append(row[c])
    return tuple(key)


def _get_float(row, col):
    try:
        return float(row[col])
    except Exception as e:
        raise ValueError(f"Column {col!r} is missing or non-numeric in row: {row}") from e


def _safe_max(values):
    # max() but returns NaN for empty/all-NaN lists
    finite = [v for v in values if not (isinstance(v, float) and math.isnan(v))]
    if not finite:
        return float("nan")
    return max(finite)


def compare_results(results_file, refdata_file, tolerances=None):
    """
    Compare results CSV with reference data CSV.
    
    Parameters:
    -----------
    results_file : str
        Path to the results CSV file
    refdata_file : str
        Path to the reference data CSV file
    tolerances : dict, optional
        Dictionary mapping column names to absolute or relative tolerances.
        Example: {'ss': {'rtol': 0.01}, 'exectime': {'atol': 0.1}}
    
    Returns:
    --------
    bool : True if all comparisons pass, False otherwise
    """

    # Default tolerances
    if tolerances is None:
        tolerances = {
            'ss': {'rtol': 1.5e-2},  # 1% relative tolerance for supersaturation
            'th_diff': {'atol': 1e-5},  # absolute tolerance for theta leak
            'rv_diff': {'atol': 1e-6},  # absolute tolerance for rv leak
            'act': {'rtol': 1.5e-2},  # relative tolerance for activated concentration
            'mr': {'rtol': 1.5e-2},  #  relative tolerance for mean radius
            'sr': {'rtol': 1.5e-2},  #  relative tolerance for second moment
            'tr': {'rtol': 1.5e-2},  #  relative tolerance for third moment
            'act_post_evap': {'rtol': 1.5e-2},
            'gccn_post_evap': {'rtol': 1.5e-2},
            'th_post_cond': {'rtol': 1e-4},
            'rv_post_cond': {'rtol': 1e-3},
        }

    # Load data
    try:
        results_rows = _read_csv_rows(results_file)
        ref_rows = _read_csv_rows(refdata_file)
    except FileNotFoundError as e:
        print(f"Error: Could not find file: {e}")
        return False

    # Check if both CSVs have the same number of rows
    if len(results_rows) != len(ref_rows):
        print(f"Error: Row count mismatch. Results: {len(results_rows)}, Reference: {len(ref_rows)}")
        return False

    # Check columns
    results_cols = set(results_rows[0].keys()) if results_rows else set()
    ref_cols = set(ref_rows[0].keys()) if ref_rows else set()
    if results_cols != ref_cols:
        print("Error: Column mismatch.")
        print(f"Results columns: {sorted(results_cols)}")
        print(f"Reference columns: {sorted(ref_cols)}")
        return False

    config_cols = ['mixing', 'constp', 'exact_sstp', 'adaptive', 'RH_formula', 'sstp_cond', 'sstp_cond_act']

    results_map = {}
    for r in results_rows:
        k = _config_key(r, config_cols)
        if k in results_map:
            print(f"Error: Duplicate configuration in results: {dict(zip(config_cols, k))}")
            return False
        results_map[k] = r

    ref_map = {}
    for r in ref_rows:
        k = _config_key(r, config_cols)
        if k in ref_map:
            print(f"Error: Duplicate configuration in reference: {dict(zip(config_cols, k))}")
            return False
        ref_map[k] = r

    if set(results_map.keys()) != set(ref_map.keys()):
        missing_in_ref = set(results_map.keys()) - set(ref_map.keys())
        missing_in_results = set(ref_map.keys()) - set(results_map.keys())
        if missing_in_ref:
            print(f"Error: {len(missing_in_ref)} configurations missing in reference.")
        if missing_in_results:
            print(f"Error: {len(missing_in_results)} configurations missing in results.")
        return False

    # Compare each variable
    all_pass = True
    comparison_cols = [col for col in tolerances.keys() if col in results_cols]

    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)

    for col in comparison_cols:
        tol = tolerances[col]

        # Calculate differences
        if 'rtol' in tol:
            # Relative tolerance
            rtol = tol['rtol']
            rel_diffs = []
            worst = None  # (rel_diff, key, result, ref)
            for k in results_map.keys():
                r_res = _get_float(results_map[k], col)
                r_ref = _get_float(ref_map[k], col)
                if r_ref == 0:
                    rd = float("nan")
                else:
                    rd = abs((r_res - r_ref) / r_ref)
                rel_diffs.append(rd)
                if not (isinstance(rd, float) and math.isnan(rd)):
                    if worst is None or rd > worst[0]:
                        worst = (rd, k, r_res, r_ref)
            max_rel_diff = _safe_max(rel_diffs)
            passed = all((d <= rtol) for d in rel_diffs if not (isinstance(d, float) and math.isnan(d)))
            if isinstance(max_rel_diff, float) and math.isnan(max_rel_diff):
                passed = True

            print(f"\n{col}:")
            print(f"  Max relative difference: {max_rel_diff:.6e} (tolerance: {rtol:.6e})")

            if not passed:
                print(f"  FAILED: Relative difference exceeds tolerance")
                failing_count = sum(1 for d in rel_diffs if not (isinstance(d, float) and math.isnan(d)) and d > rtol)
                print(f"  Number of failing configurations: {failing_count}")
                print(f"  Example failing case:")
                if worst is not None:
                    rd, k, r_res, r_ref = worst
                    print(f"    Config: {dict(zip(config_cols, k))}")
                    print(f"    Result: {r_res:.6e}")
                    print(f"    Reference: {r_ref:.6e}")
                    print(f"    Relative diff: {rd:.6e}")
                all_pass = False
            else:
                print(f"  PASSED")

        elif 'atol' in tol:
            # Special handling for rv_diff and th_diff: check if |result| <= |ref| + tolerance
            if col in ['rv_diff', 'th_diff']:
                # For leakage variables, we want to verify that the absolute value doesn't increase
                # It's acceptable if |result| <= |ref| * (1 + 0.015), i.e., within 1.5% increase
                tolerance_margin = 0.015  # 1.5%
                passed = True
                max_result = -float("inf")
                max_ref = -float("inf")
                worst = None  # (excess, key, abs_res, abs_ref, max_allowed)
                for k in results_map.keys():
                    res_abs = abs(_get_float(results_map[k], col))
                    ref_abs = abs(_get_float(ref_map[k], col))
                    max_allowed = ref_abs * (1 + tolerance_margin)
                    max_result = max(max_result, res_abs)
                    max_ref = max(max_ref, ref_abs)
                    if res_abs > max_allowed:
                        passed = False
                        excess = res_abs - max_allowed
                        if worst is None or excess > worst[0]:
                            worst = (excess, k, res_abs, ref_abs, max_allowed)

                print(f"\n{col}:")
                print(f"  Max |result|: {max_result:.6e}")
                print(f"  Max |reference|: {max_ref:.6e}")
                print(f"  Tolerance: |ref| * (1 + {tolerance_margin})")

                if not passed:
                    print(f"  FAILED: |result| exceeds |reference| * (1 + {tolerance_margin}) for some configurations")
                    failing_count = 0
                    for k in results_map.keys():
                        res_abs = abs(_get_float(results_map[k], col))
                        ref_abs = abs(_get_float(ref_map[k], col))
                        if res_abs > ref_abs * (1 + tolerance_margin):
                            failing_count += 1
                    print(f"  Number of failing configurations: {failing_count}")
                    print(f"  Example failing case:")
                    if worst is not None:
                        excess, k, res_abs, ref_abs, max_allowed = worst
                        print(f"    Config: {dict(zip(config_cols, k))}")
                        print(f"    |result|: {res_abs:.6e}")
                        print(f"    |reference|: {ref_abs:.6e}")
                        print(f"    Max allowed: {max_allowed:.6e}")
                        print(f"    Excess: {excess:.6e}")
                    all_pass = False
                else:
                    print(f"  PASSED: |result| <= |reference| * (1 + {tolerance_margin}) for all configurations")
            else:
                # Standard absolute tolerance for other variables
                atol = tol['atol']
                diffs = []
                worst = None  # (abs_diff, key, result, ref)
                for k in results_map.keys():
                    r_res = _get_float(results_map[k], col)
                    r_ref = _get_float(ref_map[k], col)
                    d = abs(r_res - r_ref)
                    diffs.append(d)
                    if worst is None or d > worst[0]:
                        worst = (d, k, r_res, r_ref)
                max_abs_diff = max(diffs) if diffs else float("nan")
                passed = all(d <= atol for d in diffs)

                print(f"\n{col}:")
                print(f"  Max absolute difference: {max_abs_diff:.6e} (tolerance: {atol:.6e})")

                if not passed:
                    print(f"  FAILED: Absolute difference exceeds tolerance")
                    failing_count = sum(1 for d in diffs if d > atol)
                    print(f"  Number of failing configurations: {failing_count}")
                    print(f"  Example failing case:")
                    if worst is not None:
                        d, k, r_res, r_ref = worst
                        print(f"    Config: {dict(zip(config_cols, k))}")
                        print(f"    Result: {r_res:.6e}")
                        print(f"    Reference: {r_ref:.6e}")
                        print(f"    Absolute diff: {d:.6e}")
                    all_pass = False
                else:
                    print(f"  PASSED")

    # Additional check: act_post_evap should equal gccn_post_evap
    # After evaporation, only larger mode particles should have r > 0.5 microns
    if 'act_post_evap' in results_cols and 'gccn_post_evap' in results_cols:
        print("\n" + "-"*80)
        print("ADDITIONAL CHECK: act_post_evap == gccn_post_evap")
        print("-"*80)

        # Allow small tolerance for floating point comparison
        equality_tol = 1e-10
        diffs = []
        worst = None  # (diff, key, act, gccn)
        for k in results_map.keys():
            act = _get_float(results_map[k], 'act_post_evap')
            gccn = _get_float(results_map[k], 'gccn_post_evap')
            d = abs(act - gccn)
            diffs.append(d)
            if worst is None or d > worst[0]:
                worst = (d, k, act, gccn)
        max_diff = max(diffs) if diffs else float("nan")
        equality_check = all(d <= equality_tol for d in diffs)

        print(f"Max absolute difference: {max_diff:.6e} (tolerance: {equality_tol:.6e})")

        if not equality_check:
            print(f"FAILED: act_post_evap != gccn_post_evap for some configurations")
            failing_count = sum(1 for d in diffs if d > equality_tol)
            print(f"Number of failing configurations: {failing_count}")
            print(f"Example failing case:")
            if worst is not None:
                d, k, act, gccn = worst
                print(f"  Config: {dict(zip(config_cols, k))}")
                print(f"  act_post_evap: {act:.6e}")
                print(f"  gccn_post_evap: {gccn:.6e}")
                print(f"  Difference: {d:.6e}")
            all_pass = False
        else:
            print(f"PASSED: act_post_evap == gccn_post_evap for all configurations")

    # Additional check: supersaturation should be in expected range
    # GCCNs condense even at ss<0
    if 'ss' in results_cols:
        print("\n" + "-"*80)
        print("ADDITIONAL CHECK: ss within expected range")
        print("-"*80)

        ss_min = -0.96
        ss_max = -0.71
        ss_values = []
        worst_low = None  # (ss, key)
        worst_high = None  # (ss, key)
        for k in results_map.keys():
            v = _get_float(results_map[k], 'ss')
            ss_values.append(v)
            if worst_low is None or v < worst_low[0]:
                worst_low = (v, k)
            if worst_high is None or v > worst_high[0]:
                worst_high = (v, k)
        min_ss = min(ss_values) if ss_values else float("nan")
        max_ss = max(ss_values) if ss_values else float("nan")
        range_check = all((v >= ss_min and v <= ss_max) for v in ss_values)

        print(f"Expected range: [{ss_min}, {ss_max}]")
        print(f"Actual range:   [{min_ss:.3f}, {max_ss:.3f}]")

        if not range_check:
            print(f"FAILED: ss outside expected range for some configurations")
            failing_count = sum(1 for v in ss_values if v < ss_min or v > ss_max)
            print(f"Number of failing configurations: {failing_count}")

            if worst_low is not None and worst_low[0] < ss_min:
                v, k = worst_low
                print(f"Example case with ss too low:")
                print(f"  Config: {dict(zip(config_cols, k))}")
                print(f"  ss: {v:.6f} (min allowed: {ss_min})")

            if worst_high is not None and worst_high[0] > ss_max:
                v, k = worst_high
                print(f"Example case with ss too high:")
                print(f"  Config: {dict(zip(config_cols, k))}")
                print(f"  ss: {v:.6f} (max allowed: {ss_max})")

            all_pass = False
        else:
            print(f"PASSED: ss within expected range for all configurations")

    print("\n" + "="*80)
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*80 + "\n")

    return all_pass


def main():
    """Main test function"""
    
    # Define file paths
    results_file = "test_results/lgrngn_cond_substepping_results.csv"

    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))    
    # Define file paths relative to script location
    refdata_file = os.path.join(script_dir, "refdata", "lgrngn_cond_substepping_refdata.csv")
    
    
    # Check if reference data exists
    if not os.path.exists(refdata_file):
        print(f"Warning: Reference data file not found: {refdata_file}")
        print("To create reference data, run:")
        print(f"  cp {results_file} {refdata_file}")
        print("\nThis will use the current results as the reference for future tests.")
        return 1
    
    # Check if results file exists
    if not os.path.exists(results_file):
        print(f"Error: Results file not found: {results_file}")
        print("Please run lgrngn_cond_substepping.py first to generate results.")
        return 1
    
    # Run comparison
    success = compare_results(results_file, refdata_file)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
