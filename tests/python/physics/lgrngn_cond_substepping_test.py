#!/usr/bin/env python3
"""
Test script to compare condensation substepping results against reference data.
This script runs the lgrngn_cond_substepping.py test and compares the results
with previously stored reference data.
"""

import sys
import os
import pandas as pd
import numpy as np

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
        results = pd.read_csv(results_file)
        refdata = pd.read_csv(refdata_file)
    except FileNotFoundError as e:
        print(f"Error: Could not find file: {e}")
        return False
    
    # Check if both dataframes have the same shape
    if results.shape != refdata.shape:
        print(f"Error: Shape mismatch. Results: {results.shape}, Reference: {refdata.shape}")
        return False
    
    # Check if both dataframes have the same columns
    if not set(results.columns) == set(refdata.columns):
        print(f"Error: Column mismatch.")
        print(f"Results columns: {sorted(results.columns)}")
        print(f"Reference columns: {sorted(refdata.columns)}")
        return False
    
    # Merge on configuration columns to align rows
    config_cols = ['mixing', 'constp', 'exact_sstp', 'adaptive', 'RH_formula', 'sstp_cond', 'sstp_cond_act']
    merged = results.merge(refdata, on=config_cols, suffixes=('_result', '_ref'))
    
    if len(merged) != len(results):
        print(f"Error: Not all configurations matched. {len(merged)} out of {len(results)} matched.")
        return False
    
    # Compare each variable
    all_pass = True
    comparison_cols = [col for col in tolerances.keys() if col in results.columns]
    
    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)
    
    for col in comparison_cols:
        result_col = f"{col}_result"
        ref_col = f"{col}_ref"
        
        if result_col not in merged.columns or ref_col not in merged.columns:
            continue
        
        tol = tolerances[col]
        
        # Calculate differences
        if 'rtol' in tol:
            # Relative tolerance
            rtol = tol['rtol']
            # Avoid division by zero
            ref_values = merged[ref_col].replace(0, np.nan)
            rel_diff = np.abs((merged[result_col] - merged[ref_col]) / ref_values)
            max_rel_diff = rel_diff.max()
            passed = (rel_diff <= rtol).all() or np.isnan(max_rel_diff)
            
            print(f"\n{col}:")
            print(f"  Max relative difference: {max_rel_diff:.6e} (tolerance: {rtol:.6e})")
            
            if not passed:
                print(f"  FAILED: Relative difference exceeds tolerance")
                failing_rows = merged[rel_diff > rtol]
                print(f"  Number of failing configurations: {len(failing_rows)}")
                print(f"  Example failing case:")
                if len(failing_rows) > 0:
                    idx = rel_diff.idxmax()
                    print(f"    Config: {merged.loc[idx, config_cols].to_dict()}")
                    print(f"    Result: {merged.loc[idx, result_col]:.6e}")
                    print(f"    Reference: {merged.loc[idx, ref_col]:.6e}")
                    print(f"    Relative diff: {rel_diff.loc[idx]:.6e}")
                all_pass = False
            else:
                print(f"  PASSED")
        
        elif 'atol' in tol:
            # Special handling for rv_diff and th_diff: check if |result| <= |ref| + tolerance
            if col in ['rv_diff', 'th_diff']:
                # For leakage variables, we want to verify that the absolute value doesn't increase
                # It's acceptable if |result| <= |ref| * (1 + 0.015), i.e., within 1.5% increase
                ref_abs = np.abs(merged[ref_col])
                result_abs = np.abs(merged[result_col])
                tolerance_margin = 0.015  # 1.5%
                max_allowed = ref_abs * (1 + tolerance_margin)
                passed = (result_abs <= max_allowed).all()
                max_result = result_abs.max()
                max_ref = ref_abs.max()
                
                print(f"\n{col}:")
                print(f"  Max |result|: {max_result:.6e}")
                print(f"  Max |reference|: {max_ref:.6e}")
                print(f"  Tolerance: |ref| * (1 + {tolerance_margin})")
                
                if not passed:
                    print(f"  FAILED: |result| exceeds |reference| * (1 + {tolerance_margin}) for some configurations")
                    failing_mask = result_abs > max_allowed
                    failing_rows = merged[failing_mask]
                    print(f"  Number of failing configurations: {len(failing_rows)}")
                    print(f"  Example failing case:")
                    if len(failing_rows) > 0:
                        idx = (result_abs - max_allowed).idxmax()
                        print(f"    Config: {merged.loc[idx, config_cols].to_dict()}")
                        print(f"    |result|: {result_abs.loc[idx]:.6e}")
                        print(f"    |reference|: {ref_abs.loc[idx]:.6e}")
                        print(f"    Max allowed: {max_allowed.loc[idx]:.6e}")
                        print(f"    Excess: {(result_abs.loc[idx] - max_allowed.loc[idx]):.6e}")
                    all_pass = False
                else:
                    print(f"  PASSED: |result| <= |reference| * (1 + {tolerance_margin}) for all configurations")
            else:
                # Standard absolute tolerance for other variables
                atol = tol['atol']
                abs_diff = np.abs(merged[result_col] - merged[ref_col])
                max_abs_diff = abs_diff.max()
                passed = (abs_diff <= atol).all()
                
                print(f"\n{col}:")
                print(f"  Max absolute difference: {max_abs_diff:.6e} (tolerance: {atol:.6e})")
                
                if not passed:
                    print(f"  FAILED: Absolute difference exceeds tolerance")
                    failing_rows = merged[abs_diff > atol]
                    print(f"  Number of failing configurations: {len(failing_rows)}")
                    print(f"  Example failing case:")
                    if len(failing_rows) > 0:
                        idx = abs_diff.idxmax()
                        print(f"    Config: {merged.loc[idx, config_cols].to_dict()}")
                        print(f"    Result: {merged.loc[idx, result_col]:.6e}")
                        print(f"    Reference: {merged.loc[idx, ref_col]:.6e}")
                        print(f"    Absolute diff: {abs_diff.loc[idx]:.6e}")
                    all_pass = False
                else:
                    print(f"  PASSED")
    
    # Additional check: act_post_evap should equal gccn_post_evap
    # After evaporation, only larger mode particles should have r > 0.5 microns
    if 'act_post_evap' in results.columns and 'gccn_post_evap' in results.columns:
        print("\n" + "-"*80)
        print("ADDITIONAL CHECK: act_post_evap == gccn_post_evap")
        print("-"*80)
        
        act_post_evap = merged['act_post_evap_result']
        gccn_post_evap = merged['gccn_post_evap_result']
        
        # Allow small tolerance for floating point comparison
        equality_tol = 1e-10
        abs_diff = np.abs(act_post_evap - gccn_post_evap)
        max_diff = abs_diff.max()
        equality_check = (abs_diff <= equality_tol).all()
        
        print(f"Max absolute difference: {max_diff:.6e} (tolerance: {equality_tol:.6e})")
        
        if not equality_check:
            print(f"FAILED: act_post_evap != gccn_post_evap for some configurations")
            failing_rows = merged[abs_diff > equality_tol]
            print(f"Number of failing configurations: {len(failing_rows)}")
            print(f"Example failing case:")
            if len(failing_rows) > 0:
                idx = abs_diff.idxmax()
                print(f"  Config: {merged.loc[idx, config_cols].to_dict()}")
                print(f"  act_post_evap: {act_post_evap.loc[idx]:.6e}")
                print(f"  gccn_post_evap: {gccn_post_evap.loc[idx]:.6e}")
                print(f"  Difference: {abs_diff.loc[idx]:.6e}")
            all_pass = False
        else:
            print(f"PASSED: act_post_evap == gccn_post_evap for all configurations")
    
    # Additional check: supersaturation should be in expected range
    # GCCNs condense even at ss<0
    if 'ss' in results.columns:
        print("\n" + "-"*80)
        print("ADDITIONAL CHECK: ss within expected range")
        print("-"*80)
        
        ss_min = -0.96
        ss_max = -0.71
        ss_values = merged['ss_result']
        
        min_ss = ss_values.min()
        max_ss = ss_values.max()
        range_check = (ss_values >= ss_min).all() and (ss_values <= ss_max).all()
        
        print(f"Expected range: [{ss_min}, {ss_max}]")
        print(f"Actual range:   [{min_ss:.3f}, {max_ss:.3f}]")
        
        if not range_check:
            print(f"FAILED: ss outside expected range for some configurations")
            failing_rows_low = merged[ss_values < ss_min]
            failing_rows_high = merged[ss_values > ss_max]
            failing_count = len(failing_rows_low) + len(failing_rows_high)
            print(f"Number of failing configurations: {failing_count}")
            
            if len(failing_rows_low) > 0:
                idx = ss_values.idxmin()
                print(f"Example case with ss too low:")
                print(f"  Config: {merged.loc[idx, config_cols].to_dict()}")
                print(f"  ss: {ss_values.loc[idx]:.6f} (min allowed: {ss_min})")
            
            if len(failing_rows_high) > 0:
                idx = ss_values.idxmax()
                print(f"Example case with ss too high:")
                print(f"  Config: {merged.loc[idx, config_cols].to_dict()}")
                print(f"  ss: {ss_values.loc[idx]:.6f} (max allowed: {ss_max})")
            
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
    refdata_file = "refdata/lgrngn_cond_substepping_refdata.csv"
    
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
