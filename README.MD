# 🧪 Möbius–Gauss Lab Suite

[![License](https://img.shields.io/badge/License-Research%20Use-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue?logo=c%2B%2B)](https://isocpp.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://python.org)
[![PARI/GP](https://img.shields.io/badge/PARI%2FGP-2.17.3-red?logo=pari)](https://pari.math.u-bordeaux.fr/)
[![OpenMP](https://img.shields.io/badge/OpenMP-4.5-green?logo=openmp)](https://www.openmp.org/)
[![Paper](https://img.shields.io/badge/Paper-View%20PDF-brightgreen)](https://arxiv.org/abs/xxxx.xxxxx)

---

This repository contains the **code, data, and scripts** used in the paper

> ****Spectral Resonances and the Möbius Function over ℤ[i]****  
> *by Vitor M. Pozza*  
> *to appear in Mathematical Proceedings of the Cambridge Philosophical Society*

The project focuses on the **computational and numerical verification** of analytic number theory results concerning the Möbius function \(\mu(\alpha)\) on the Gaussian integers \(\mathbb{Z}[i]\).

---

## 🎯 Overview

The main objectives of this project are:

- ✅ **Exact enumeration of μ(α) for all non-zero Gaussian integers with norm N(α) ≤ 10⁸.
- ✅ **Validation of the truncated explicit formula for the summatory function S(X) = Σ μ(α), with N(α) ≤ X.
- ✅ **Numerical analysis** of angular equidistribution.
- ✅ **Computation of shifted autocorrelations** (Chowla-type statistics).
- ✅ **Distribution of the values −1, 0 and 1 of μ.
- ✅ **High-precision data for the zeros of ζ(s) and L(s, χ₋₄) and their derivatives.

---

## 📁 Repository Structure

```text
mobius-gauss-project/
├── mobius_tables.cpp                  # C++ sieve implementation
├── mobius_tables.exe                  # compiled executable (Windows)
├── Makefile.mak                       # build script
├── generate_data.gp                   # PARI/GP script for zeros/derivatives
├── validate_explicit_formula.py       # explicit formula validation
├── analyze_mobius.py                  # main statistical analysis
├── angular_analysis.py                # angular equidistribution analysis
├── requirements.txt                   # Python dependencies
├── README.MD                          # this file
├── figures/                           # generated plots
│   ├── figure1_summatory.pdf
│   ├── figure2_autocorrelation.pdf
│   ├── figure3_distribution.pdf
│   ├── figure4_explicit_formula_amplitude.pdf
│   ├── figure5_angular_equidistribution.pdf
├── mobius_data_1e8.csv                # full sieve data (≈2.96 GB)
├── table1_summatory.csv               # Table 1 data
├── table2_autocorr.csv                # Table 2 data
├── table3_distribution.csv            # Table 3 proportions
├── table3_distribution_full.csv       # Table 3 full counts
├── angular_data.csv                   # angular deviation data
├── explicit_validation_comparative.csv # Table IV comparisons
├── zeros_zeta.txt                     # zeta zeros (partial)
├── zeros_L4.txt                       # L(s,chi_-4) zeros (partial)
├── zeros_zeta_raw.txt                 # zeta zeros up to height 10050
├── zeros_L4_raw.txt                   # L(s,chi_-4) zeros up to height 10050
├── zetaprime_zeros_zeta_prec.txt      # zeta'(rho) at zeta zeros
├── L_zeros_zeta_prec.txt              # L(rho,chi_-4) at zeta zeros
├── zeta_zeros_L4_prec.txt             # zeta(rho) at L-zeros
├── Lprime_zeros_L4_prec.txt           # L'(rho,chi_-4) at L-zeros
└── zeta_zeros_zeta_prec.txt           # auxiliary zeta(rho) at zeta zeros

🔍 Detailed File Descriptions
⚙️ C++ Implementation
File	Description
mobius_tables.cpp	High-performance sieve in C++17 with OpenMP parallelisation. Generates all Gaussian primes with norm ≤ 104104, initialises a lattice grid, marks prime multiples, tracks prime factor parity, detects non-square-free elements, and writes μ(α)μ(α) for all 0<N(α)≤1080<N(α)≤108.
mobius_tables.exe	Compiled executable for Windows (produced from the C++ source).
Makefile.mak	Build script for compiling with a C++17 compiler (requires OpenMP support).
🧮 PARI/GP Scripts
File	Description
generate_data.gp	Computes the non-trivial zeros of ζ(s)ζ(s) and L(s,χ−4)L(s,χ−4​) up to height 10,05010,050, evaluates derivatives and values, and writes high-precision data files for the Python validation script.
🐍 Python Scripts
File	Description
validate_explicit_formula.py	Loads zero data, constructs the coefficients 4/(ρζQ(i)′(ρ))4/(ρζQ(i)′​(ρ)), evaluates the truncated explicit formula with T=XT=X
​, compares with sieve values, and generates Table IV and Figure 4.
analyze_mobius.py	Reads the full sieve data, computes S(X)S(X), C(X)C(X), and the distribution of μμ values; produces Tables 1–3 and Figures 1–3.
angular_analysis.py	Computes angular deviations for fixed sectors (width 1∘1∘) and produces the data for Figure 5.
requirements.txt	Python package dependencies (see below).
📊 Data Files
Full Sieve Data
File	Description
mobius_data_1e8.csv	Complete enumeration of μ(α) for all non-zero Gaussian integers with N(α) ≤ 10⁸. Large file (approx. 2.96 GB). Format: a,b,mu, where a, b ∈ Z, a² + b² ≤ 10⁸, and mu = μ(a + bi).
Summary Tables
File	Description
table1_summatory.csv	Data for the summatory function S(X) at X = 10⁵, 10⁶, 10⁷, 10⁸.
table2_autocorr.csv	Shifted autocorrelation C(X)=∑μ(α)μ(α+1)C(X)=∑μ(α)μ(α+1).
table3_distribution.csv	Proportions of μ(α)=−1,0,1μ(α)=−1,0,1.
table3_distribution_full.csv	Full counts of μμ values.
Angular Analysis
File	Description
angular_data.csv	Angular deviation Δθ(X)Δθ​(X) for sectors of width 1∘1∘ at X=108X=108. Used to generate Figure 5.
Explicit Formula Validation
File	Description
explicit_validation_comparative.csv	Comparison between sieve values and finite spectral reconstructions (Table IV).
Zeros and Spectral Data
File	Description
zeros_zeta.txt / zeros_L4.txt	Partial lists of zeros (used in early validation).
zeros_zeta_raw.txt / zeros_L4_raw.txt	Full lists of zeros computed by PARI/GP up to height 10,050.
zetaprime_zeros_zeta_prec.txt	ζ′(1/2+iγ)ζ′(1/2+iγ) at zeros of ζ(s)ζ(s). Format: gamma Re(zetaprime) Im(zetaprime).
L_zeros_zeta_prec.txt	L(1/2+iγ,χ−4)L(1/2+iγ,χ−4​) at zeros of ζ(s)ζ(s).
zeta_zeros_L4_prec.txt	ζ(1/2+iγ)ζ(1/2+iγ) at zeros of L(s,χ−4)L(s,χ−4​).
Lprime_zeros_L4_prec.txt	L′(1/2+iγ,χ−4)L′(1/2+iγ,χ−4​) at zeros of L(s,χ−4)L(s,χ−4​).
zeta_zeros_zeta_prec.txt	Auxiliary file with ζ(1/2+iγ)ζ(1/2+iγ) at zeros of ζ(s)ζ(s).
📦 Data Sources

    The zeros of ζ(s)ζ(s) used in the validation were initially taken from Odlyzko's tables:

        A. M. Odlyzko, On the distribution of spacings between zeros of the zeta function, Math. Comp. 48 (1987), 273–308.

    The final high-precision zeros and derivatives were generated with PARI/GP using the script generate_data.gp, up to height 10,05010,050.

    All integer arithmetic in the sieve was performed with strict 64-bit integers to avoid floating-point errors.

🚀 Computational Pipeline

The complete numerical workflow consists of the following steps:

    Sieve computation
    Compile and run the C++ code:
    bash

    make -f Makefile.mak
    ./mobius_tables.exe

    This produces mobius_data_1e8.csv.

    Main analysis
    Run the Python analysis script:
    bash

    python analyze_mobius.py

    This generates Tables 1–3 and Figures 1–3.

    Generation of spectral data
    Run the PARI/GP script:
    bash

    gp -q generate_data.gp

    This produces the high-precision zero and derivative files.

    Explicit formula validation
    Run the Python validation script:
    bash

    python validate_explicit_formula.py

    This produces Table IV and Figure 4.

    Angular equidistribution analysis
    Run:
    bash

    python angular_analysis.py

    This produces angular_data.csv and Figure 5.

🛠 Requirements
C++

    C++17 compiler with OpenMP support (e.g., GCC, MSVC, Clang).

    The compiled Windows executable mobius_tables.exe is also included.

PARI/GP

    PARI/GP version 2.17.3 or later.

Python

The following Python packages are required (see requirements.txt):
text

numpy
pandas
matplotlib
seaborn
mpmath

Install them with:
bash

pip install -r requirements.txt

📝 Reproducibility Notes

    The complete sieve data file `mobius_data_1e8.csv` is approximately 2.96 GB
    and is therefore not hosted directly in this repository. It is available
    for download from Zenodo: [(https://doi.org/10.5281/zenodo.22182052)]

    The high-precision zero files were generated with PARI/GP and contain more digits than necessary for most analyses; they are provided for complete reproducibility.

    The validation scripts assume that all data files are located in the same directory as the scripts.

📄 Citation

If you use this code or data, please cite the associated paper:
bibtex

@article{Pozza2026,
  author  = {Vitor M. Pozza},
  title   = {Spectral Resonances and the M\"obius Function over {\(\mathbb{Z}[i]\)}},
  journal = {Math. Proc. Camb. Phil. Soc.},
  year    = {2026},
  note    = {to appear}
}

📧 Contact

    Author: Vitor M. Pozza

    Email: pozza.data.analyst@outlook.com

    Institution: Sosígenes Costa Institute – Federal University of Southern Bahia

⚖️ License

Copyright (c) 2026 Vitor M. Pozza. All rights reserved.

This repository, including all scripts, theoretical algorithms, and processed datasets, is provided strictly for academic research, peer review, and scientific reproducibility purposes. 

No commercial use, redistribution, or modification of the software and data is permitted without explicit prior written authorization from the author.
