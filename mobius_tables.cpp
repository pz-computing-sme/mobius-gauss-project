/**
 * ============================================================================
 * Gaussian Integer Möbius Sieve - FINAL WITH COMPLETE VALIDATIONS (FIXED)
 * ============================================================================
 * 
 * Computes μ(α) for all Gaussian integers α = a + bi with norm N(α) ≤ 10⁸.
 * Additionally, aggregates μ(α) by angular sector (bins) for anisotropy analysis,
 * and exports individual (a, b, mu) data for the Möbius phase diagram.
 * 
 * VALIDATIONS IMPLEMENTED IN THIS C++ CODE (integer arithmetic):
 *   1. Conjugate symmetry: μ(a,b) == μ(a,-b) [checked, not fixed]
 *   2. Unit symmetry: μ(α) == μ(iα) for all four units
 *   3. Multiplicativity: μ(αβ) == μ(α)μ(β) for 10,000 random coprime pairs
 *   4. Density convergence: compared with 1/ζ_Q(i)(2)
 * 
 * ============================================================================
 * IMPORTANT NOTE: COMPARISON WITH THE EXPLICIT FORMULA (Theorem 5.3)
 * ============================================================================
 * The explicit formula validation (Table IV and Figure 4 of the article) is
 * NOT implemented in this C++ sieve, and FOR GOOD REASON:
 * 
 *   - The explicit formula requires evaluating ζ(s) and L(s, χ₋₄) near the
 *     critical line Re(s) = 1/2 with 50+ decimal digits of precision.
 *   - Standard C++ double precision (≈15 digits) is wholly inadequate and
 *     would produce massive numerical errors due to cancellation.
 *   - Implementing arbitrary-precision complex zeta/L-functions in C++
 *     from scratch is impractical and would obscure the core sieve logic.
 * 
 * This validation is CORRECTLY performed by the external Python scripts:
 *   - generate_data_explicit.py  : loads zeros from PARI/GP, computes L(ρ)
 *   - validate_explicit_formula.py : evaluates the explicit sum and
 *                                     produces Table IV and Figure 4.
 * 
 * These scripts use the mpmath library (arbitrary precision) and are the
 * appropriate tools for this analytic validation. The user MUST keep and
 * execute them separately to reproduce the explicit formula comparison.
 * ============================================================================
 * 
 * Compilation:
 *   g++ -O3 -fopenmp -std=c++17 -o mobius_tables mobius_tables.cpp
 * 
 * Execution:
 *   ./mobius_tables
 * 
 * Author: Vitor M. Pozza
 * Date: 2026
 * ============================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <cstdint>
#include <iomanip>
#include <algorithm>
#include <random>
#include <numeric>

using namespace std;

constexpr int8_t MOBIUS_ZERO = 0;
constexpr int8_t MOBIUS_POS = 1;
constexpr int8_t MOBIUS_NEG = -1;

const double PI = 3.14159265358979323846;
const double CATALAN = 0.915965594177219;

// Number of angular bins (360 = 1 degree per bin)
const int NUM_BINS = 360;

// Number of random pairs for multiplicativity test
const int NUM_MULTIPLICITY_TESTS = 10000;

// ============================================================================
// Gaussian integer representation
// ============================================================================
struct Gaussian {
    int a, b;
    Gaussian(int a_ = 0, int b_ = 0) : a(a_), b(b_) {}
    
    long long norm() const { return 1LL * a * a + 1LL * b * b; }
    
    bool operator<(const Gaussian& other) const {
        long long n1 = norm();
        long long n2 = other.norm();
        if (n1 != n2) return n1 < n2;
        if (a != other.a) return a < other.a;
        return b < other.b;
    }
    
    bool operator==(const Gaussian& other) const {
        return a == other.a && b == other.b;
    }
    
    Gaussian operator*(const Gaussian& other) const {
        return Gaussian(a * other.a - b * other.b, a * other.b + b * other.a);
    }
};

// ============================================================================
// Cornacchia's algorithm: solve x² + y² = p for p ≡ 1 mod 4
// ============================================================================
pair<int, int> cornacchia(int p) {
    // Find a square root of -1 modulo p
    int x0 = -1;
    for (long long a = 2; a < p; ++a) {
        if ((a * a) % p == p - 1) {
            x0 = (int)a;
            break;
        }
    }
    if (x0 == -1) return {0, 0};

    // Euclidean reduction
    int r0 = p, r1 = x0;
    while (r1 * r1 >= p) {
        int r2 = r0 % r1;
        r0 = r1;
        r1 = r2;
    }
    int x = r1;
    int y = (int)sqrt(p - x * x);
    if (x * x + y * y != p) return {0, 0};
    return {x, y};
}

// ============================================================================
// Generate all Gaussian primes with norm ≤ max_norm
// ============================================================================
vector<Gaussian> generate_gaussian_primes(int max_norm) {
    vector<Gaussian> primes;
    int R = max_norm;

    // Sieve for rational primes up to R
    vector<bool> is_prime(R + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= R; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= R; j += i)
                is_prime[j] = false;
        }
    }

    // Ramified prime (1+i) with norm 2
    if (2 <= max_norm) primes.push_back(Gaussian(1, 1));

    int split_count = 0, inert_count = 0;

    for (int p = 2; p <= R; ++p) {
        if (!is_prime[p]) continue;
        if (p == 2) continue;

        if (p % 4 == 3) {
            if (1LL * p * p <= max_norm) {
                primes.push_back(Gaussian(p, 0));
                inert_count++;
            }
        } else if (p % 4 == 1) {
            if (p <= max_norm) {
                auto [x, y] = cornacchia(p);
                if (x > 0 && y > 0) {
                    primes.push_back(Gaussian(x, y));
                    primes.push_back(Gaussian(x, -y));
                    split_count++;
                }
            }
        }
    }

    sort(primes.begin(), primes.end());

    cout << "Split primes (p ≡ 1 mod 4) added: " << split_count << " (each contributes 2 Gaussian primes)" << endl;
    cout << "Inert primes (p ≡ 3 mod 4) added: " << inert_count << endl;
    cout << "Ramified prime (1+i): 1" << endl;
    cout << "Total Gaussian primes: " << primes.size() << endl;

    return primes;
}

// ============================================================================
// Check divisibility: is (a+bi) divisible by (pa+pb i)?
// ============================================================================
bool is_divisible(int a, int b, int pa, int pb, long long p_norm) {
    long long num_real = 1LL * a * pa + 1LL * b * pb;
    long long num_imag = 1LL * b * pa - 1LL * a * pb;
    return (num_real % p_norm == 0) && (num_imag % p_norm == 0);
}

// ============================================================================
// MAIN
// ============================================================================
int main() {
    const long long MAX_NORM = 100000000;      // 10⁸
    const int MAX_PRIME_NORM = 10000;          // primes with norm ≤ 10000
    const int R = (int)ceil(sqrt(MAX_NORM));   // 10000

    cout << "========================================" << endl;
    cout << "Gaussian Integer Möbius Sieve (Complete Validation)" << endl;
    cout << "========================================" << endl;
    cout << "Maximum norm: " << MAX_NORM << endl;
    cout << "Prime norm limit: " << MAX_PRIME_NORM << endl;
    cout << "Radius: " << R << endl;
    cout << "Angular bins: " << NUM_BINS << endl;
    cout << "OpenMP threads: " << omp_get_max_threads() << endl;
    cout << "========================================\n" << endl;

    auto start_total = chrono::high_resolution_clock::now();

    // 1. Generate Gaussian primes
    cout << "Generating Gaussian primes..." << endl;
    vector<Gaussian> primes = generate_gaussian_primes(MAX_PRIME_NORM);

    // 2. Allocate arrays (OpenMP‑friendly types)
    int grid_size = 2 * R + 1;
    long long total_cells = 1LL * grid_size * grid_size;
    cout << "Grid size: " << grid_size << " x " << grid_size << endl;
    cout << "Estimated memory: " << (total_cells * (sizeof(short) + sizeof(char))) / (1024.0 * 1024.0) << " MB" << endl;

    vector<short> factor_count(total_cells, 0);
    vector<char> has_square(total_cells, 0);
    vector<int8_t> mu(total_cells, 0);

    // --- Angular accumulators (one per bin) ---
    vector<long long> angular_sum_mu(NUM_BINS, 0);
    vector<long long> angular_count(NUM_BINS, 0);

    auto index = [&](int a, int b) -> long long {
        return (a + R) * grid_size + (b + R);
    };

    // Helper to get bin index from angle (in radians)
    auto angle_to_bin = [&](double theta) -> int {
        if (theta < 0) theta += 2 * PI;
        int bin = (int)(theta / (2 * PI) * NUM_BINS);
        if (bin >= NUM_BINS) bin = NUM_BINS - 1;
        return bin;
    };

    // 3. Sieve: factor each point
    cout << "\nRunning sieve..." << endl;
    auto start_sieve = chrono::high_resolution_clock::now();

    const int progress_step = 1000;
    int next_progress_row = -R + progress_step;
    long long total_points = 0;

    #pragma omp parallel for collapse(2) reduction(+:total_points)
    for (int a = -R; a <= R; ++a) {
        for (int b = -R; b <= R; ++b) {
            if (a == 0 && b == 0) continue;
            long long norm = 1LL * a * a + 1LL * b * b;
            if (norm > MAX_NORM) continue;

            total_points++;

            int idx = index(a, b);
            int local_factor_count = 0;
            bool local_has_square = false;
            long long temp_a = a, temp_b = b;
            long long remaining_norm = norm;

            for (size_t pi = 0; pi < primes.size(); ++pi) {
                const Gaussian& p = primes[pi];
                long long p_norm = p.norm();

                if (p_norm > norm) break;
                if (p_norm > remaining_norm && remaining_norm > 1) continue;

                if (!is_divisible((int)temp_a, (int)temp_b, p.a, p.b, p_norm)) continue;

                int mult = 0;
                while (is_divisible((int)temp_a, (int)temp_b, p.a, p.b, p_norm)) {
                    mult++;
                    long long num_real = temp_a * p.a + temp_b * p.b;
                    long long num_imag = temp_b * p.a - temp_a * p.b;
                    temp_a = num_real / p_norm;
                    temp_b = num_imag / p_norm;
                    remaining_norm /= p_norm;
                }

                if (mult >= 2) {
                    local_has_square = true;
                    break;
                } else if (mult == 1) {
                    local_factor_count++;
                }

                if (remaining_norm == 1) break;
            }

            if (remaining_norm > 1) {
                local_factor_count++;
            }

            int8_t mu_val;
            if (local_has_square) {
                mu_val = MOBIUS_ZERO;
            } else {
                mu_val = (local_factor_count % 2 == 0) ? MOBIUS_POS : MOBIUS_NEG;
            }

            factor_count[idx] = local_factor_count;
            has_square[idx] = local_has_square ? 1 : 0;
            mu[idx] = mu_val;

            double theta = atan2((double)b, (double)a);
            int bin = angle_to_bin(theta);
            #pragma omp atomic
            angular_sum_mu[bin] += mu_val;
            #pragma omp atomic
            angular_count[bin]++;
        }

        if (a >= next_progress_row) {
            int percent = (a + R) * 100 / (2 * R);
            cout << "Progress: " << percent << "%" << endl;
            next_progress_row = a + progress_step;
        }
    }

    auto end_sieve = chrono::high_resolution_clock::now();
    auto sieve_time = chrono::duration_cast<chrono::seconds>(end_sieve - start_sieve).count();
    cout << "Sieve completed in " << sieve_time << " seconds." << endl;
    cout << "Total points: " << total_points << endl;

    // ================================================================
    // 4. VALIDATION: Conjugate Symmetry (CHECK ONLY, no fixing)
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "VALIDATION 1: Conjugate Symmetry" << endl;
    cout << "========================================" << endl;
    
    long long conj_violations = 0;
    long long conj_checked = 0;
    
    for (int a = -R; a <= R; ++a) {
        for (int b = 1; b <= R; ++b) {
            long long norm = 1LL * a * a + 1LL * b * b;
            if (norm > MAX_NORM) continue;
            conj_checked++;
            
            int idx1 = index(a, b);
            int idx2 = index(a, -b);
            
            if (mu[idx1] != mu[idx2]) {
                conj_violations++;
            }
        }
    }
    
    cout << "Points checked: " << conj_checked << endl;
    cout << "Violations found: " << conj_violations << endl;
    
    if (conj_violations == 0) {
        cout << "✓ PASS: Conjugate symmetry holds exactly for all checked points." << endl;
    } else {
        cout << "✗ FAIL: Conjugate symmetry violated for " << conj_violations << " points." << endl;
        cout << "  (This indicates a bug in the sieve implementation.)" << endl;
    }
    cout << "========================================\n" << endl;

    // ================================================================
    // 5. VALIDATION: Unit Symmetry (μ(α) == μ(iα) for all units)
    // ================================================================
    cout << "========================================" << endl;
    cout << "VALIDATION 2: Unit Symmetry" << endl;
    cout << "========================================" << endl;
    
    long long unit_violations = 0;
    long long unit_checked = 0;
    
    for (int a = 1; a <= R; ++a) {
        for (int b = 0; b <= R; ++b) {
            long long norm = 1LL * a * a + 1LL * b * b;
            if (norm > MAX_NORM) continue;
            if (a == 0 && b == 0) continue;
            
            int idx_alpha = index(a, b);
            int idx_ia = index(-b, a);   // i * (a + bi) = -b + ai
            int idx_ma = index(-a, -b);  // -1 * (a + bi) = -a - bi
            int idx_maia = index(b, -a); // -i * (a + bi) = b - ai
            
            int m1 = mu[idx_alpha];
            int m2 = mu[idx_ia];
            int m3 = mu[idx_ma];
            int m4 = mu[idx_maia];
            
            if (m1 != m2 || m1 != m3 || m1 != m4) {
                unit_violations++;
            }
            unit_checked++;
        }
    }
    
    cout << "Points checked: " << unit_checked << endl;
    cout << "Violations found: " << unit_violations << endl;
    
    if (unit_violations == 0) {
        cout << "✓ PASS: Unit symmetry holds exactly for all checked points." << endl;
    } else {
        cout << "✗ FAIL: Unit symmetry violated for " << unit_violations << " points." << endl;
        cout << "  (This indicates a bug in the sieve implementation.)" << endl;
    }
    cout << "========================================\n" << endl;

    // ================================================================
    // 6. VALIDATION: Multiplicativity (random coprime pairs)
    // ================================================================
    cout << "========================================" << endl;
    cout << "VALIDATION 3: Multiplicativity" << endl;
    cout << "========================================" << endl;
    
    // Collect only points with norm <= sqrt(MAX_NORM) = 10000
    // so that product of two norms fits into MAX_NORM
    const long long SMALL_LIMIT = (long long)sqrt(MAX_NORM); // 10000
    vector<Gaussian> small_nonzero_points;
    small_nonzero_points.reserve(total_points / 100); // rough estimate
    
    for (int a = -R; a <= R; ++a) {
        for (int b = -R; b <= R; ++b) {
            if (a == 0 && b == 0) continue;
            long long norm = 1LL * a * a + 1LL * b * b;
            if (norm > MAX_NORM) continue;
            if (norm > SMALL_LIMIT) continue; // only keep small points
            int idx = index(a, b);
            if (mu[idx] != 0) {
                small_nonzero_points.push_back(Gaussian(a, b));
            }
        }
    }
    
    cout << "Small non-zero points available (norm <= 10000): " << small_nonzero_points.size() << endl;
    
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<size_t> dist(0, small_nonzero_points.size() - 1);
    
    long long mult_violations = 0;
    long long mult_tested = 0;
    long long mult_skipped = 0;
    
    int max_tests = min(NUM_MULTIPLICITY_TESTS, (int)small_nonzero_points.size() / 2);
    
    cout << "Testing " << max_tests << " random coprime pairs (ensuring product norm <= MAX_NORM)..." << endl;
    
    for (int test = 0; test < max_tests; ++test) {
        size_t i1 = dist(gen);
        size_t i2 = dist(gen);
        
        Gaussian alpha = small_nonzero_points[i1];
        Gaussian beta = small_nonzero_points[i2];
        
        if (alpha.norm() <= 1 || beta.norm() <= 1) {
            mult_skipped++;
            continue;
        }
        
        long long n1 = alpha.norm();
        long long n2 = beta.norm();
        if (std::gcd(n1, n2) != 1) {
            mult_skipped++;
            continue;
        }
        
        // Product norm is guaranteed <= MAX_NORM because both norms <= 10000
        Gaussian product = alpha * beta;
        long long prod_norm = product.norm();
        // Safety check (should never fail)
        if (prod_norm > MAX_NORM) {
            mult_skipped++;
            continue;
        }
        
        int mu_alpha = mu[index(alpha.a, alpha.b)];
        int mu_beta = mu[index(beta.a, beta.b)];
        int mu_product = mu[index(product.a, product.b)];
        
        int expected = mu_alpha * mu_beta;
        if (mu_product != expected) {
            mult_violations++;
        }
        mult_tested++;
    }
    
    cout << "Tests performed: " << mult_tested << endl;
    cout << "Skipped (non-coprime): " << mult_skipped << endl;
    cout << "Violations found: " << mult_violations << endl;
    
    if (mult_violations == 0) {
        cout << "✓ PASS: Multiplicativity holds for all tested coprime pairs." << endl;
    } else {
        cout << "✗ FAIL: Multiplicativity violated for " << mult_violations << " pairs." << endl;
        cout << "  (This indicates a bug in the sieve implementation.)" << endl;
    }
    cout << "========================================\n" << endl;

    // 7. Verification for small norms
    cout << "\n=== Verifying Möbius values for small norms ===" << endl;

    auto check = [&](int a, int b, int expected, const string& name) {
        int val = (int)mu[index(a, b)];
        string status = (val == expected) ? "✓ PASS" : "✗ FAIL";
        cout << "μ(" << name << ") = " << val << " (expected " << expected << ") " << status << endl;
    };

    check(1, 0, 1, "1");
    check(1, 1, -1, "1+i");
    check(2, 0, 0, "2");
    check(3, 0, -1, "3");
    check(0, 3, -1, "3i");
    check(2, 1, -1, "2+i");
    check(1, 2, -1, "1+2i");
    check(2, -1, -1, "2-i");
    check(1, -2, -1, "1-2i");
    check(5, 0, 1, "5");
    check(3, 1, 1, "3+i");
    check(4, 0, 0, "4");
    check(6, 0, 0, "6");
    check(9, 0, 0, "9");
    check(10, 0, 0, "10");
    check(7, 0, -1, "7");
    check(1, 3, 1, "1+3i");

    cout << "=================================================\n" << endl;

    // 8. Statistics at checkpoints
    const long long checkpoints[4] = {100000, 1000000, 10000000, 100000000};
    long long S_vals[4] = {0, 0, 0, 0};
    long long total_counts[4] = {0, 0, 0, 0};
    long long count0[4] = {0, 0, 0, 0};
    long long count1[4] = {0, 0, 0, 0};
    long long countm1[4] = {0, 0, 0, 0};

    for (int a = -R; a <= R; ++a) {
        for (int b = -R; b <= R; ++b) {
            if (a == 0 && b == 0) continue;
            long long norm = 1LL * a * a + 1LL * b * b;
            if (norm > MAX_NORM) continue;
            int val = mu[index(a, b)];
            for (int k = 0; k < 4; ++k) {
                if (norm <= checkpoints[k]) {
                    S_vals[k] += val;
                    total_counts[k]++;
                    if (val == 0) count0[k]++;
                    else if (val == 1) count1[k]++;
                    else countm1[k]++;
                }
            }
        }
    }

    // 9. Autocorrelation C(X) for h = 1
    long long C_vals[4] = {0, 0, 0, 0};
    for (int a = -R; a <= R; ++a) {
        for (int b = -R; b <= R; ++b) {
            long long norm1 = 1LL * a * a + 1LL * b * b;
            if (norm1 == 0 || norm1 > MAX_NORM) continue;
            int mu1 = mu[index(a, b)];

            int a2 = a + 1;
            long long norm2 = 1LL * a2 * a2 + 1LL * b * b;
            if (norm2 == 0 || norm2 > MAX_NORM) continue;
            int mu2 = mu[index(a2, b)];

            long long prod = mu1 * mu2;
            long long max_norm_pair = (norm1 > norm2) ? norm1 : norm2;
            for (int k = 0; k < 4; ++k) {
                if (max_norm_pair <= checkpoints[k]) {
                    C_vals[k] += prod;
                }
            }
        }
    }

    // 10. Export tables
    cout << "\n=== Table I: Summatory function S(X) ===" << endl;
    cout << " X          S(X)          S(X)/√X          (log X)²" << endl;
    cout << "----------------------------------------------------" << endl;
    ofstream file1("table1_summatory.csv");
    file1 << "X,S(X),S(X)/sqrt(X),(log X)^2\n";
    for (int k = 0; k < 4; ++k) {
        double X = checkpoints[k];
        double S = S_vals[k];
        double ratio = S / sqrt(X);
        double log2 = pow(log(X), 2);
        cout << scientific << setprecision(4) << X << "  " << S << "  " << ratio << "  " << log2 << endl;
        file1 << X << "," << S << "," << ratio << "," << log2 << "\n";
    }
    cout << "----------------------------------------------------" << endl;
    file1.close();

    cout << "\n=== Table II: Autocorrelation Sums C(X) ===" << endl;
    cout << " X         C(X)          C(X)/X        C(X)/X^{0.85}" << endl;
    cout << "-------------------------------------------------" << endl;
    ofstream file2("table2_autocorr.csv");
    file2 << "X,C(X),C(X)/X,C(X)/X^0.85\n";
    for (int k = 0; k < 4; ++k) {
        double X = checkpoints[k];
        double C = C_vals[k];
        double ratio1 = C / X;
        double ratio2 = C / pow(X, 0.85);
        cout << scientific << setprecision(4) << X << "  " << C << "  " << ratio1 << "  " << ratio2 << endl;
        file2 << X << "," << C << "," << ratio1 << "," << ratio2 << "\n";
    }
    cout << "-------------------------------------------------" << endl;
    file2.close();

    // Table III: Distribution of μ
    double zeta_Q2 = (PI * PI) / 6 * CATALAN;
    double square_free_density = 1.0 / zeta_Q2;
    double theo0 = 1.0 - square_free_density;
    double theo1 = square_free_density / 2.0;

    cout << "\n=== Table III: Distribution of μ(α) ===" << endl;
    cout << " X          μ=0       μ=+1      μ=-1      Theoretical (0.336/0.332/0.332)" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    ofstream file3_full("table3_distribution_full.csv");
    file3_full << "X,mu0,mu1,mum1,theo0,theo1\n";
    for (int k = 0; k < 4; ++k) {
        double X = checkpoints[k];
        double p0 = (double)count0[k] / total_counts[k];
        double p1 = (double)count1[k] / total_counts[k];
        double pm1 = (double)countm1[k] / total_counts[k];
        cout << scientific << setprecision(4) << X << "  " << p0 << "  " << p1
             << "  " << pm1 << "  " << theo0 << "/" << theo1 << "/" << theo1 << endl;
        file3_full << X << "," << p0 << "," << p1 << "," << pm1 << ","
                   << theo0 << "," << theo1 << "\n";
    }
    cout << "----------------------------------------------------------------------" << endl;
    file3_full.close();

    // Final summary for graph
    ofstream file3_final("table3_distribution.csv");
    file3_final << "mu,Observed,Theoretical\n";
    int last = 3;
    double final_p0 = (double)count0[last] / total_counts[last];
    double final_p1 = (double)count1[last] / total_counts[last];
    double final_pm1 = (double)countm1[last] / total_counts[last];
    file3_final << "0," << final_p0 << "," << theo0 << "\n";
    file3_final << "+1," << final_p1 << "," << theo1 << "\n";
    file3_final << "-1," << final_pm1 << "," << theo1 << "\n";
    file3_final.close();

    // ================================================================
    // 11. VALIDATION: Density Convergence
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "VALIDATION 4: Density Convergence" << endl;
    cout << "========================================" << endl;
    
    cout << "Theoretical square-free density: 1/ζ_Q(i)(2) = " << square_free_density << endl;
    cout << "Empirical square-free density at X = " << checkpoints[last] << ": " << (1.0 - final_p0) << endl;
    
    double density_diff = fabs((1.0 - final_p0) - square_free_density);
    double expected_error = 1.0 / sqrt((double)checkpoints[last]) * log((double)checkpoints[last]);
    
    cout << "Difference: " << density_diff << endl;
    cout << "Theoretical error bound O(X^{-1/2} log X): " << expected_error << endl;
    
    if (density_diff < 10 * expected_error) {
        cout << "✓ PASS: Density convergence consistent with theoretical bound." << endl;
    } else {
        cout << "⚠ WARNING: Density deviation larger than expected." << endl;
        cout << "  (This may be due to finite-size effects or a bug.)" << endl;
    }
    cout << "========================================\n" << endl;

    // 12. Export angular data
    cout << "\n=== Angular data aggregation ===" << endl;
    ofstream ang_file("angular_data.csv");
    if (!ang_file.is_open()) {
        cerr << "Error: could not open angular_data.csv for writing." << endl;
    } else {
        ang_file << "angle_deg,sum_mu,count,avg_mu\n";
        for (int bin = 0; bin < NUM_BINS; ++bin) {
            double angle_deg = (double)bin / NUM_BINS * 360.0;
            long long sum = angular_sum_mu[bin];
            long long cnt = angular_count[bin];
            double avg = (cnt > 0) ? (double)sum / cnt : 0.0;
            ang_file << fixed << setprecision(4) << angle_deg << ","
                     << sum << ","
                     << cnt << ","
                     << avg << "\n";
        }
        ang_file.close();
        cout << "Angular data written to angular_data.csv (" << NUM_BINS << " bins)." << endl;
    }

    // 13. Export individual data for phase diagram
    cout << "\n=== Exporting individual data for phase diagram ===" << endl;
    ofstream data_file("mobius_data_1e8.csv");
    if (!data_file.is_open()) {
        cerr << "Error: could not open mobius_data_1e8.csv for writing." << endl;
    } else {
        data_file << "a,b,mu\n";
        long long exported = 0;
        for (int a = -R; a <= R; ++a) {
            for (int b = -R; b <= R; ++b) {
                if (a == 0 && b == 0) continue;
                long long norm = 1LL * a * a + 1LL * b * b;
                if (norm > MAX_NORM) continue;
                int idx = index(a, b);
                int val = (int)mu[idx];
                if (val != 0) {
                    data_file << a << "," << b << "," << val << "\n";
                    exported++;
                }
            }
        }
        data_file.close();
        cout << "Exported " << exported << " non-zero points to mobius_data_1e8.csv" << endl;
        cout << "(Zeros omitted; they do not contribute to the Möbius phase sum.)" << endl;
    }

    // ================================================================
    // 14. FINAL VALIDATION SUMMARY
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "VALIDATION SUMMARY" << endl;
    cout << "========================================" << endl;
    
    bool all_passed = true;
    
    cout << "1. Conjugate symmetry: " << (conj_violations == 0 ? "✓ PASS" : "✗ FAIL") << endl;
    if (conj_violations != 0) all_passed = false;
    
    cout << "2. Unit symmetry:      " << (unit_violations == 0 ? "✓ PASS" : "✗ FAIL") << endl;
    if (unit_violations != 0) all_passed = false;
    
    cout << "3. Multiplicativity:   " << (mult_violations == 0 ? "✓ PASS" : "✗ FAIL") << endl;
    if (mult_violations != 0) all_passed = false;
    
    cout << "4. Density convergence: " << (density_diff < 10 * expected_error ? "✓ PASS" : "⚠ WARNING") << endl;
    if (density_diff >= 10 * expected_error) all_passed = false;
    
    cout << "========================================" << endl;
    
    if (all_passed) {
        cout << "✓ ALL VALIDATIONS PASSED." << endl;
        cout << "  The computed Möbius values are numerically reliable." << endl;
    } else {
        cout << "✗ SOME VALIDATIONS FAILED." << endl;
        cout << "  Please review the sieve implementation." << endl;
    }
    
    // ================================================================
    // 15. NOTE ABOUT EXPLICIT FORMULA VALIDATION
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "EXPLICIT FORMULA VALIDATION (Table IV / Figure 4)" << endl;
    cout << "========================================" << endl;
    cout << "This validation is NOT performed in this C++ code." << endl;
    cout << "It is correctly handled by the external Python scripts:" << endl;
    cout << "  - generate_data_explicit.py" << endl;
    cout << "  - validate_explicit_formula.py" << endl;
    cout << "These scripts use the mpmath library for arbitrary precision" << endl;
    cout << "and are the appropriate tools for this analytic comparison." << endl;
    cout << "========================================\n" << endl;

    auto end_total = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::seconds>(end_total - start_total).count();

    cout << "\n========================================" << endl;
    cout << "EXECUTION SUMMARY" << endl;
    cout << "========================================" << endl;
    cout << "Total points processed: " << total_points << endl;
    cout << "Square-free density: " << 100.0 * (1.0 - final_p0) << "%" << endl;
    cout << "Total execution time: " << total_time << " seconds (" << total_time/60.0 << " minutes)" << endl;
    cout << "========================================" << endl;

    return 0;
}