"""
Validation of the Explicit Formula: Comparative Analysis

Corrections applied:
  - Coefficient +4 in the explicit formula.
  - Constant -16 from the pole at s = 0.
  - Common height filter between zeta(s) and L(s, chi_-4).

The script also prints the maximum heights of each list to
verify that the L(s, chi_-4) data covers the full required range.
"""

import mpmath as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------------------------------------
# 1. Configuration and precision
# ------------------------------------------------------------
mp.dps = 50

def load_data(filename):
    """Load a three-column file: gamma  Re(value)  Im(value) with full precision."""
    gammas = []
    values = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Sanitize PARI/GP space-padded scientific notation
            line = line.replace(" E-", "E-").replace(" E+", "E+")
            line = line.replace(" e-", "e-").replace(" e+", "e+")

            parts = line.split()
            if len(parts) >= 3:
                g = mp.mpf(parts[0])
                re = mp.mpf(parts[1])
                im = mp.mpf(parts[2])
                gammas.append(g)
                values.append(mp.mpc(re, im))
    return gammas, values

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------
print("Loading high-precision data...")
try:
    gammas_L, zeta_L = load_data("zeta_zeros_L4_prec.txt")
    _, Lprime_L = load_data("Lprime_zeros_L4_prec.txt")
    gammas_Z, zetaprime_Z = load_data("zetaprime_zeros_zeta_prec.txt")
    _, L_Z = load_data("L_zeros_zeta_prec.txt")
except FileNotFoundError as e:
    print(f"Error: file not found: {e.filename}")
    exit(1)

print(f"Total zeros of L(s, chi_-4): {len(gammas_L)}")
print(f"Total zeros of zeta(s):      {len(gammas_Z)}")

# ------------------------------------------------------------
# 3. Maximum height verification
# ------------------------------------------------------------
max_L = max(gammas_L)
max_Z = max(gammas_Z)
print(f"\nMaximum height of L(s, chi_-4): {max_L}")
print(f"Maximum height of zeta(s):       {max_Z}")

# Common height: use the smaller one for full coverage
T_common = min(max_L, max_Z)
print(f"Common height chosen: T_common = {T_common}")

# ------------------------------------------------------------
# 4. Filter both lists to g <= T_common
# ------------------------------------------------------------
mask_L = [g <= T_common for g in gammas_L]
gammas_L   = [g for g, m in zip(gammas_L, mask_L) if m]
zeta_L     = [z for z, m in zip(zeta_L, mask_L) if m]
Lprime_L   = [lp for lp, m in zip(Lprime_L, mask_L) if m]

mask_Z = [g <= T_common for g in gammas_Z]
gammas_Z   = [g for g, m in zip(gammas_Z, mask_Z) if m]
zetaprime_Z = [zp for zp, m in zip(zetaprime_Z, mask_Z) if m]
L_Z        = [lv for lv, m in zip(L_Z, mask_Z) if m]

print(f"\nAfter filtering:")
print(f"Zeros of L(s, chi_-4): {len(gammas_L)}")
print(f"Zeros of zeta(s):      {len(gammas_Z)}")

if max_L < 10000:
    print("\nWARNING: The L(s, chi_-4) list does not reach height 10000.")
    print("For validating X = 10^8, you must generate more zeros (Option A).")
else:
    print("\nThe L(s, chi_-4) list covers up to height 10000; OK for X = 10^8.")

# ------------------------------------------------------------
# 5. Precompute coefficients
# ------------------------------------------------------------
coeffs = []   # (gamma, coefficient) with coeff = 4 / (rho * zeta_K'(rho))

# Zeros of L(s, chi_-4): zeta_K'(rho) = zeta(rho) * L'(rho)
for g, z, Lp in zip(gammas_L, zeta_L, Lprime_L):
    rho = mp.mpc(0.5, g)
    zetaK_prime = z * Lp
    coeff = 4 / (rho * zetaK_prime)
    coeffs.append((g, coeff))

# Zeros of zeta(s): zeta_K'(rho) = zeta'(rho) * L(rho)
for g, zp, Lv in zip(gammas_Z, zetaprime_Z, L_Z):
    rho = mp.mpc(0.5, g)
    zetaK_prime = zp * Lv
    coeff = 4 / (rho * zetaK_prime)
    coeffs.append((g, coeff))

coeffs.sort(key=lambda x: x[0])
print(f"\nPrecomputed {len(coeffs)} coefficients.")

# ------------------------------------------------------------
# 6. Explicit sum function
# ------------------------------------------------------------
def explicit_sum(X, use_moving_truncation=False):
    """
    Computes the explicit formula sum.
    If use_moving_truncation is True, cuts off frequencies above T = sqrt(X).
    """
    X_mp = mp.mpf(X)
    total = mp.mpf(-16)   # contribution from the pole at s = 0

    T = mp.sqrt(X_mp) if use_moving_truncation else None

    for g, coeff in coeffs:
        if T is not None and g > T:
            break
        rho = mp.mpc(0.5, g)
        term = coeff * (X_mp ** rho)
        total += 2 * mp.re(term)

    return total

# ------------------------------------------------------------
# 7. Evaluation at checkpoints
# ------------------------------------------------------------
X_vals = [1e5, 1e6, 1e7, 1e8]
S_abs_sieve = [172, 2020, 5884, 1184]

print("\nEvaluating explicit formula...")
results = []
for X, S in zip(X_vals, S_abs_sieve):
    S_exp_trunc = explicit_sum(X, use_moving_truncation=True)
    abs_trunc = float(abs(S_exp_trunc))
    err_trunc = abs(abs_trunc - S) / S

    S_exp_full = explicit_sum(X, use_moving_truncation=False)
    abs_full = float(abs(S_exp_full))
    err_full = abs(abs_full - S) / S

    results.append((X, S, abs_trunc, err_trunc, abs_full, err_full))

# ------------------------------------------------------------
# 8. Save and print table
# ------------------------------------------------------------
columns = ["X", "|S(X)| Sieve", "|S_exp| T=sqrt(X)", "Err Trunc", "|S_exp| Full", "Err Full"]
df = pd.DataFrame(results, columns=columns)
df.to_csv("explicit_validation_comparative.csv", index=False)
print("\nResults saved to 'explicit_validation_comparative.csv'.")

print("\n" + "="*95)
print("                           COMPARATIVE VALIDATION TABLE IV")
print("===============================================================================================")
print(df.to_string(index=False, formatters={
    "X": lambda x: f"{x:.0e}",
    "|S(X)| Sieve": lambda x: f"{int(x):d}",
    "|S_exp| T=sqrt(X)": lambda x: f"{x:.2f}",
    "Err Trunc": lambda x: f"{x:.4f}",
    "|S_exp| Full": lambda x: f"{x:.2f}",
    "Err Full": lambda x: f"{x:.4f}"
}))
print("===============================================================================================\n")

# ------------------------------------------------------------
# 9. Generate plot
# ------------------------------------------------------------
print("Generating dense plot using dynamic truncation...")
X_plot = np.logspace(5, 8, 30)
S_plot_trunc = []

for X in X_plot:
    S_exp = explicit_sum(X, use_moving_truncation=True)
    S_plot_trunc.append(float(abs(S_exp)))

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)
fig, ax = plt.subplots(figsize=(9, 6))

ax.loglog(X_plot, S_plot_trunc, 'b-', linewidth=2.5, label=r'$|S_{\text{explicit}}(X)|$ (with $T=\sqrt{X}$)')
ax.loglog(X_plot, np.sqrt(X_plot), 'k--', linewidth=1.8, label=r'$X^{1/2}$ reference scale')

x_vals_float = [float(x) for x in X_vals]
ax.scatter(x_vals_float, [r[2] for r in results], color='red', s=90, zorder=5, label='Sieve Checkpoints (Truncated)')

ax.set_xlabel(r'$X$', fontsize=16)
ax.set_ylabel('Amplitude', fontsize=16)
ax.legend(loc="upper left")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
plt.tight_layout()

plt.savefig("explicit_validation_precise.pdf", dpi=300)
print("Plot saved as 'explicit_validation_precise.pdf'")