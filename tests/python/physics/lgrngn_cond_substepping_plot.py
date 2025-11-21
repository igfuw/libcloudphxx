import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the results
df = pd.read_csv("test_results/lgrngn_cond_substepping_results.csv")

variables = [
    ('ss', 'Supersaturation [%]'),
    ('th_diff', 'theta leaked [K]'),
    ('rv_diff', 'rv leaked [1]'),
    ('act', 'Activated concentration [1/g]'),
    ('mr', 'Mean radius [um]'),
    ('sr', 'Second moment of radius [m^2]'),
    ('tr', 'Third moment of radius [m^3]'),
    ('exectime', 'Execution time [s]'),
]

RH_formula_names = {
    0 : "pv_cc", 1: "rv_cc", 2 :"pv_tet", 3 : "rv_tet"
}

group_cols = ['mixing', 'constp', 'exact_sstp', 'adaptive', 'sstp_cond_act']
RH_formulas = sorted(df['RH_formula'].unique())

os.makedirs("test_results/plots", exist_ok=True)

for var, ylabel in variables:
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    axs = axs.flatten()
    for i, RH_formula in enumerate(RH_formulas):
        ax = axs[i]
        subdf = df[df['RH_formula'] == RH_formula]
        for name, group in subdf.groupby(group_cols):
            label = ', '.join(f"{col}={val}" for col, val in zip(group_cols, name))
            ax.plot(group['sstp_cond'], group[var], marker='o', label=label, linestyle='--')
        ax.set_title(f'RH_formula={RH_formula_names[RH_formula]}')
        ax.set_xlabel('Number of substeps')
        ax.set_ylabel(ylabel)
        ax.set_xscale('log')
        if var == 'exectime':
            ax.set_yscale('log')
        ax.legend(fontsize='small', loc='best')
    plt.tight_layout()
    plt.savefig(f"test_results/plots/{var}_vs_substeps.png")
    plt.close()