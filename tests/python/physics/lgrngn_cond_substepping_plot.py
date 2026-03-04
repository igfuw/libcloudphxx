import csv
import matplotlib.pyplot as plt
import os


def _read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _to_bool(v):
    s = str(v).strip().lower()
    if s in ("true", "1", "t", "yes", "y"):
        return True
    if s in ("false", "0", "f", "no", "n"):
        return False
    raise ValueError(f"Cannot parse boolean from: {v!r}")


def _to_float(v):
    return float(str(v).strip())


def _to_int(v):
    return int(float(str(v).strip()))


# Load the results
rows = _read_csv_rows("test_results/lgrngn_cond_substepping_results.csv")

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
RH_formulas = sorted({ _to_int(r['RH_formula']) for r in rows })

os.makedirs("test_results/plots", exist_ok=True)

def _group_key(row):
    return (
        _to_bool(row['mixing']),
        _to_bool(row['constp']),
        _to_bool(row['exact_sstp']),
        _to_bool(row['adaptive']),
        _to_int(row['sstp_cond_act']),
    )


for var, ylabel in variables:
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    axs = axs.flatten()
    for i, RH_formula in enumerate(RH_formulas):
        ax = axs[i]
        subrows = [r for r in rows if _to_int(r['RH_formula']) == RH_formula]

        grouped = {}
        for r in subrows:
            k = _group_key(r)
            grouped.setdefault(k, []).append(r)

        for key, group in grouped.items():
            label = ', '.join(
                f"{col}={val}" for col, val in zip(group_cols, key)
            )
            group_sorted = sorted(group, key=lambda r: _to_int(r['sstp_cond']))
            x = [_to_int(r['sstp_cond']) for r in group_sorted]
            y = [_to_float(r[var]) for r in group_sorted]
            ax.plot(x, y, marker='o', label=label, linestyle='--')
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