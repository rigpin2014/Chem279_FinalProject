#include <iostream>             // For std::cout and std::endl
#include <fstream>              // For std::ifstream
#include <filesystem>           // For JSON config stuff
#include <string>               // For std::string
#include <iomanip>              // For std::fixed and std::setprecision
#include <cstdlib>              // For JSON config stuff
#include <stdexcept>            // For exceptions
#include <vector>               // For std::vector
#include <map>                  // For std::map
#include <utility>              // For std::pair
#include <cmath>                // For std::erf
#include <nlohmann/json.hpp>    // This is the JSON handling library
#include <armadillo>            // For eigenvalue stuff

// convenience definitions so the code is more readable
namespace fs = std::filesystem;
using json = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////////////////
//================================
// Define some custom data types
//================================
////////////////////////////////////////////////////////////////////////////////////////////

// Define the atom
struct Atom {
    int atomic_number;
    std::vector<double> coord;

    Atom(int atomic_number, double x, double y, double z)
        : atomic_number(atomic_number), coord{x, y, z} {}
};

// The basis functions are called contracted functions and are written as a linear
// combination of primitive gaussian functions
struct PG {             // Primitive Gaussian
    double alpha_k;     // alpha_k
    double d_k;         // d_uk
    double norm;        // N_k
    double atom_index;   // atom 0, 1, ...

    PG(double alpha, double d, double n = 0.0, int idx = -1)
    : alpha_k(alpha), d_k(d), norm(n), atom_index(idx) {}
};

struct BF {                         // Basis Function
    int atomic_number;
    std::vector<double> r;          // (x, y, z)
    double l, m, n;                 // angular momentum components
    std::vector<PG> pg;             // primitive gaussians


    BF(int atomic_number, std::vector<double> r, int l, int m, int n, std::vector<PG> pg)
        : atomic_number(atomic_number), r(r), l(l), m(m), n(n), pg(pg) {}
};

/////////////////////////////////////////////////////////////////////////////////////////
//================================
// Define the primitive gaussians
//================================
////////////////////////////////////////////////////////////////////////////////////////

// Basis input for H 1s AO
// | Exponents | 1s Contraction Ceoefficients | Normalization Constant |
const std::vector<PG> H_1s = {
    {3.42525091, 0.15432897, 0.0},
    {0.62391373, 0.53532814, 0.0},
    {0.16885540, 0.44463454, 0.0}
};

// Basis input for C 2s AO
// | Exponents | 2s Contraction Coefficients | Normalization Constant |
const std::vector<PG> C_2s = {
    {2.94124940, -0.09996723, 0.0},
    {0.68348310, 0.39951283, 0.0},
    {0.22228990, 0.70011547, 0.0}
};

// Basis input for C 2p AO
// | Exponents | 2p Contraction Coefficients | Normalization Constant |
const std::vector<PG> C_2p = {
    {2.94124940,  0.15591627, 0.0},
    {0.68348310,  0.60768372, 0.0},
    {0.22228990,  0.39195739, 0.0}
};

// Basis input for N 2s AO
// | Exponents | 2s Contraction Coefficients | Normalization Constant |
const std::vector<PG> N_2s = {
    {3.78045590, -0.09996723, 0.0},
    {0.87849660, 0.39951283, 0.0},
    {0.28571440, 0.70011547, 0.0}
};

// Basis input for N 2p AO
// | Exponents | 2p Contraction Coefficients | Normalization Constant |
const std::vector<PG> N_2p = {
    {3.78045590,  0.15591627, 0.0},
    {0.87849660,  0.60768372, 0.0},
    {0.28571440,  0.39195739, 0.0}
};

// Basis input for O 2s AO
// | Exponents | 2s Contraction Coefficients | Normalization Constant |
const std::vector<PG> O_2s = {
    {5.03315130, -0.09996723, 0.0},
    {1.16959610, 0.39951283, 0.0},
    {0.38038900, 0.70011547, 0.0}
};

// Basis input for N 2p AO
// | Exponents | 2p Contraction Coefficients | Normalization Constant |
const std::vector<PG> O_2p = {
    {5.03315130,  0.15591627, 0.0},
    {1.16959610,  0.60768372, 0.0},
    {0.38038900,  0.39195739, 0.0}
};

// Basis input for F 2s AO
// | Exponents | 2s Contraction Coefficients | Normalization Constant |
const std::vector<PG> F_2s = {
    {6.46480320, -0.09996723, 0.0},
    {1.50228120, 0.39951283, 0.0},
    {0.48858850, 0.70011547, 0.0}
};

// Basis input for N 2p AO
// | Exponents | 2p Contraction Coefficients | Normalization Constant |
const std::vector<PG> F_2p = {
    {6.46480320,  0.15591627, 0.0},
    {1.50228120,  0.60768372, 0.0},
    {0.48858850,  0.39195739, 0.0}
}; 

//////////////////////////////////////////////////////////////////////////////////
//=======================================
// Read in the coordinates of each atom
//=======================================
//////////////////////////////////////////////////////////////////////////////////

std::vector<Atom> read_atoms(const std::string& filename) {
    std::ifstream file (filename);

    int num_atoms;
    file >> num_atoms;

    std::string line;
    std::getline(file, line);
    std::getline(file, line);       // Skip the comment line

    int atomic_number;
    double x, y, z;
    std::vector<Atom> atoms;

    while (file >> atomic_number >> x >> y >> z) {
        atoms.push_back(Atom(atomic_number, x, y, z));
    }

    return atoms;
}

///////////////////////////////////////////////////////////////////////////////////
//===================================
// Build a list of the basis functions
//===================================
////////////////////////////////////////////////////////////////////////////////////

std::vector<BF> create_BFS_list (const std::vector<Atom>& atoms) {

    // Initialize the list of basis functions
    std::vector<BF> basis_functions;

    // Create a lambda function to tag the PG with the atom index
    auto tag_pg = [](const std::vector<PG>& shell, int atom_index) {
        std::vector<PG> tagged;

        // Loop through the shells and tag them
        for (const PG& pg: shell) {
            tagged.emplace_back(pg.alpha_k, pg.d_k, pg.norm, atom_index);
        }

        return tagged;
    };

    for (size_t atom_index = 0; atom_index < atoms.size(); atom_index++) {
        const Atom& atom = atoms[atom_index];
        // Define the center of the basis function
        std::vector<double> center = {atom.coord[0], atom.coord[1], atom.coord[2]};

        // Define the angular momentum components
        // Hydrogen
        if (atom.atomic_number == 1 ) {
            // 1s 
            basis_functions.push_back(BF(1, center, 0, 0, 0, tag_pg(H_1s, atom_index)));

        // Carbon
        } else if (atom.atomic_number == 6 ) {
            // 2s
            basis_functions.push_back(BF(6, center, 0, 0, 0, tag_pg(C_2s, atom_index)));
            // 2p
            basis_functions.push_back(BF(6, center, 1, 0, 0, tag_pg(C_2p, atom_index)));
            basis_functions.push_back(BF(6, center, 0, 1, 0, tag_pg(C_2p, atom_index)));
            basis_functions.push_back(BF(6, center, 0, 0, 1, tag_pg(C_2p, atom_index)));

        // Nitrogen
        } else if (atom.atomic_number == 7) {
            // 2s
            basis_functions.push_back(BF(7, center, 0, 0, 0, tag_pg(N_2s, atom_index)));
            // 2p
            basis_functions.push_back(BF(7, center, 1, 0, 0, tag_pg(N_2p, atom_index)));
            basis_functions.push_back(BF(7, center, 0, 1, 0, tag_pg(N_2p, atom_index)));
            basis_functions.push_back(BF(7, center, 0, 0, 1, tag_pg(N_2p, atom_index)));

        // Oxygen
        } else if (atom.atomic_number == 8) {
            // 2s
            basis_functions.push_back(BF(8, center, 0, 0, 0, tag_pg(O_2s, atom_index)));
            // 2p
            basis_functions.push_back(BF(8, center, 1, 0, 0, tag_pg(O_2p, atom_index)));
            basis_functions.push_back(BF(8, center, 0, 1, 0, tag_pg(O_2p, atom_index)));
            basis_functions.push_back(BF(8, center, 0, 0, 1, tag_pg(O_2p, atom_index)));

        // Flourine
        } else if (atom.atomic_number == 9) {
            // 2s
            basis_functions.push_back(BF(9, center, 0, 0, 0, tag_pg(F_2s, atom_index)));
            // 2p
            basis_functions.push_back(BF(9, center, 1, 0, 0, tag_pg(F_2p, atom_index)));
            basis_functions.push_back(BF(9, center, 0, 1, 0, tag_pg(F_2p, atom_index)));
            basis_functions.push_back(BF(9, center, 0, 0, 1, tag_pg(F_2p, atom_index)));
        }
    }

    return basis_functions;
}

////////////////////////////////////////////////////////////////////////////////////
//================================================
// Define functions to compute overlap integrals
//================================================
////////////////////////////////////////////////////////////////////////////////////

// Helper function to solve factorial (recursive)
long factorial(int f) {
    if (f < 0) {
        return 0;
    } else if (f == 0) {
        return 1;
    } else {
        return f * factorial(f - 1);
    }
}

// Helper function to solve binomials with coefficients (w, v)
int binomial(int w, int v){
    if (v > w) return 0;
    return factorial(w) / (factorial(v) * factorial(w - v));
}

// Helper function to solve double factorial
long double_factorial(int f) {
    if (f <= 1) {
        return 1;
    }
    long result = 1;
    for (int i = f; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

// Define a function to solve the 1-D overlap function analytically
double analytical_1D_overlap(const BF& bf1, const BF& bf2, const PG& pg1, const PG& pg2, int coord_dir) {

    // Extract the values from the basis functions and the primitive gaussians
    double Ri = bf1.r[coord_dir];
    double Rj = bf2.r[coord_dir];
    double alpha = pg1.alpha_k;
    double beta = pg2.alpha_k;
    int la = (coord_dir == 0) ? bf1.l : (coord_dir == 1) ? bf1.m : bf1.n;
    int lb = (coord_dir == 0) ? bf2.l : (coord_dir == 1) ? bf2.m : bf2.n;

    // Define the sub-parts of the function
    double exponential_prefactor = std::exp((-alpha * beta) / (alpha + beta) * std::pow(Ri - Rj, 2));
    double square_root = std::sqrt(M_PI / (alpha + beta));
    double Rp = ((alpha * Ri) + (beta * Rj)) / (alpha + beta);

    // Compute the summation terms
    double summation = 0.0;
    for (int i = 0; i <= la; i++) {
        for (int j = 0; j <= lb; j++) { 

            // Skip odd terms because they are equal to zero.
            if ((i + j) % 2 != 0) continue;
            
            // Obtain the binomial terms
            double binomial_terms = binomial(la, i) * binomial(lb, j);

            // Obtain the non-binomial terms
            long fa = double_factorial(i + j - 1);
            double fb = std::pow(Rp - Ri, la - i) * std::pow(Rp - Rj, lb - j);
            double fc = std::pow(2 * (alpha + beta), (i + j) / 2);

            // Add the combination of terms to the sum
            summation += binomial_terms * fa * fb / fc;
        }
    }
    /* This is print debugging stuff
    std::cout << "Ri: " << Ri << ", Rj: " << Rj << std::endl;
    std::cout << "Summation: " << summation << std::endl;
    std::cout << "Exponential prefactor: " << exponential_prefactor << ", Square root: " << square_root << std::endl;
    std::cout << "S_kl: " << exponential_prefactor * square_root * summation << ", coord dir: " << coord_dir << std::endl;
    */

    return exponential_prefactor * square_root * summation;
}

// Define a helper function to take the power of base 0
double safe_pow(double base, int exp) {
    if (base == 0.0 && exp < 0) return 0.0;
    return std::pow(base, exp);
}

// Define a function to evaluate the derivative of the 1D overlap
double derivative_1D_overlap(const BF& bf1, const BF& bf2, const PG& pg1, const PG& pg2, int coord_dir, int wrt) {

    // Extract the values from the basis functions and the primitive gaussians
    double Ri = bf1.r[coord_dir];
    double Rj = bf2.r[coord_dir];
    double alpha = pg1.alpha_k;
    double beta = pg2.alpha_k;
    double Rp = ((alpha * Ri) + (beta * Rj)) / (alpha + beta);
    int la = (coord_dir == 0) ? bf1.l : (coord_dir == 1) ? bf1.m : bf1.n;
    int lb = (coord_dir == 0) ? bf2.l : (coord_dir == 1) ? bf2.m : bf2.n;

    // Define FA1 and FA1p
    double c = (alpha * beta) / (alpha + beta);
    double FA1 = std::sqrt(M_PI / (alpha + beta)) * std::exp(-c * safe_pow(Ri - Rj, 2));
    double FA1p = -2 * c * (Ri - Rj) * FA1;
    FA1p = (wrt == 1) ? FA1p : -FA1p;           // wrt == 1 therefore dSuv/dRA, else dSuv/dRB

    // Define FA2 and FA2p
    double FA2 = 0.0;
    double FA2p = 0.0;

    // Compute the summation terms
    for (int i = 0; i <= la; i++) {
        for (int j = 0; j <= lb; j++) { 

            // Skip odd terms because they are equal to zero.
            if ((i + j) % 2 != 0) continue;
            
            // Obtain the binomial terms
            double binomial_terms = binomial(la, i) * binomial(lb, j);
            double iterative_terms = double_factorial(i + j - 1) / safe_pow(2 * (alpha + beta), (i + j) / 2);

            // Calculate FA2 this iteration
            FA2 += binomial_terms * iterative_terms * safe_pow(Rp - Ri, la - i) * safe_pow(Rp - Rj, lb - j);

            // Define derivatives of terms
            double cB = -beta / (alpha + beta);
            double cA = alpha / (alpha + beta);
            cB = (wrt == 1) ? cB : -cB;             // wrt == 1 therefore dSuv/dRA, else dSuv/dRB
            cA = (wrt == 1) ? cA : -cA;             // wrt == 1 therefore dSuv/dRA, else dSuv/dRB
            double FA2_3 = (la - i) * cB * safe_pow(Rp - Ri, la - i - 1);
            double FA2_4 = (lb - j) * cA * safe_pow(Rp - Rj, lb - j - 1);
            double FA2_2 = FA2_3 * safe_pow(Rp - Rj, lb - j) + FA2_4 * safe_pow(Rp - Ri, la - i);

            // Calculate FA2p this iteration
            FA2p += binomial_terms * iterative_terms * FA2_2;
        }
    }

    return FA1p * FA2 + FA2p * FA1;
    
}

// Define a function to evalute the 3D overlap function analytically
double analytical_3D_overlap(const BF& bf1, const BF& bf2, const PG& pg1, const PG& pg2) {

    // Compute the 3-D overlap integral by looping through the dimensions
    double integral = 1.0;
    for (int i = 0; i < 3; i++) {
    // Compute the integral in 1D
    double Si = analytical_1D_overlap(bf1, bf2, pg1, pg2, i);

    // Multiply the 1D components together
    integral *= Si;
    }

    return integral;
}

////////////////////////////////////////////////////////////////////////////////////////
//=============================================
// Function to compute nomalization constants
//=============================================
////////////////////////////////////////////////////////////////////////////////////////

// Define a function to compute the normalization constants for each primitive gaussian in each basis function
void compute_normalization_constants(std::vector<BF>& basis_functions){

    // Loop through all of the basis functions
    for (auto& bf : basis_functions) {
        // Loop through all of the primitive gaussians in the basis function
        for (auto& pg : bf.pg){
            // Compute the self overlap for the current basis function and primitive gaussian
            double S_AA = analytical_3D_overlap(bf, bf, pg, pg);

            // Compute and update the normalization constant for the current
            // basis function and primitive gaussian
            pg.norm = 1.0 / std::sqrt(S_AA); 
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////
//===============================================
// Functions to compute and print overlap matrix
//===============================================
////////////////////////////////////////////////////////////////////////////////////////

// Define a function to print armadillo matrices
void print_arma_matrix(const arma::mat& matrix) {

    // Loop through the rows
    for (size_t i = 0; i < matrix.n_rows; ++i) {
        // Loop through the columns
        for (size_t j = 0; j < matrix.n_cols; ++j) {
            std::cout << std::fixed << std::setprecision(6) << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Define a function to compute the 3D overlap matrix
std::vector<std::vector<double>> compute_contracted_overlap_matrix(const std::vector<BF>& basis_functions) {
    // Initialize an empty matrix
    int N = basis_functions.size();
    std::vector<std::vector<double>> overlap_matrix(N, std::vector<double>(N, 0.0));

    // Loop through all combinations of basis functions
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {

            // Define the value of the matrix at this iteration
            double S_wv = 0.0;

            // Loop over the primitive gaussians
            for (const auto& pg1 : basis_functions[i].pg) {
                for (const auto& pg2 : basis_functions[j].pg) {

                    // Compute the primitive unormalized overlap integral
                    double S_kl = analytical_3D_overlap(basis_functions[i], basis_functions[j], pg1, pg2);

                    /* This is print debugging stuff
                    std::cout << "Overlap matrix location (i, j): (" << i << ", " << j << ")" << std::endl;
                    std::cout << "BF1 center: (" << basis_functions[i].r[0] << ", " << basis_functions[i].r[1]
                              << ", " << basis_functions[i].r[2] << ")" << std::endl;
                    std::cout << "BF2 center: (" << basis_functions[j].r[0] << ", " << basis_functions[j].r[1]
                              << ", " << basis_functions[j].r[2] << ")" << std::endl;
                    std::cout << "BF1 alpha: " << pg1.alpha_k << std::endl;
                    std::cout << "BF2 alpha: " << pg2.alpha_k << std::endl;
                    std::cout << "pg1_d: " << pg1.d_k << "     pg2_d: " << pg2.d_k << std::endl;
                    std::cout << "pg1_norm: " << pg1.norm << "     pg2_norm: " << pg2.norm << std::endl;
                    std::cout << "S_kl: " << S_kl << std::endl;
                    

                    double yuio = pg1.d_k * pg2.d_k * pg1.norm * pg2.norm * S_kl;
                    std::cout << "Sub-total: " << yuio << std::endl;
                    std::cout << std::endl;
                    */

                    // Compute S_wv this iteration
                    S_wv += pg1.d_k * pg2.d_k * pg1.norm * pg2.norm * S_kl;
                }
            }

            // Update the contracted overlap matrix
            overlap_matrix[i][j] = S_wv;
            overlap_matrix[j][i] = S_wv;
        }
    }

    return overlap_matrix;
}

// Define a function to compute th partial derivative of the normalized contracted
// gaussian function (S) with respect to nulcear coordinates (Ra)
arma::cube compute_Suv_Ri(const std::vector<Atom>& atoms,
    const std::vector<BF>& basis_functions,
    int AorB) {

    // AorB : 1 == A, 2 == B   (The atom to take the derivative wrt)

    int N = basis_functions.size();

    // Allocate memory for the final matrix.
    arma::cube Suv_RA(3, N,  N, arma::fill::zeros);

    // Loop through all combinations of basis functions
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {

            // Create a vector that represents something
            arma::vec3 dSkl_tot = arma::zeros(3);

            // Loop over the primitive gaussians
            for (const auto& pg1 : basis_functions[i].pg) {
                for (const auto& pg2 : basis_functions[j].pg) {

                    // Compute the 1D overlaps
                    double S_klx = analytical_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 0);
                    double S_kly = analytical_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 1);
                    double S_klz = analytical_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 2);

                    double dS_klx_dXa = 0;
                    double dS_kly_dYa = 0;
                    double dS_klz_dZa = 0;

                    // Compute the derivative to the 1D overlap
                    if (basis_functions[i].r != basis_functions[j].r) {
                        dS_klx_dXa = derivative_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 0, AorB);
                        dS_kly_dYa = derivative_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 1, AorB);
                        dS_klz_dZa = derivative_1D_overlap(basis_functions[i], basis_functions[j], pg1, pg2, 2, AorB);
                    }
                    
                    // Define the normalized contraction coefficients
                    double c_norm = pg1.d_k * pg1.norm * pg2.d_k * pg2.norm;

                    // Compute the cartesian product
                    arma::vec3 dSkl = arma::zeros(3);
                    dSkl(0) = dS_klx_dXa * S_kly * S_klz;
                    dSkl(1) = S_klx * dS_kly_dYa * S_klz;
                    dSkl(2) = S_klx * S_kly * dS_klz_dZa;

                    dSkl_tot += c_norm * dSkl;
                }
            }

            // Store the results anti-symmetrically
            Suv_RA(0, i, j) = dSkl_tot(0);  // X component
            Suv_RA(1, i, j) = dSkl_tot(1);  // Y component
            Suv_RA(2, i, j) = dSkl_tot(2);  // Z component

            // By translational invariance or symmetry:
            Suv_RA(0, j, i) = -dSkl_tot(0);
            Suv_RA(1, j, i) = -dSkl_tot(1);
            Suv_RA(2, j, i) = -dSkl_tot(2);
        }
    }

    return Suv_RA;
}

// Create a function to convert nested vectors to armadillo matrices
arma::mat convert_to_arma(const std::vector<std::vector<double>>& matrix) {

    // Define an empty matrix
    arma::mat arma(matrix.size(), matrix[0].size());

    // Loop through the rows
    for (size_t i = 0; i < matrix.size(); ++i) {
        // Loop through the columns
        for (size_t j = 0; j < matrix[0].size();++j) {
            arma(i, j) = matrix[i][j];
        }
    }

    return arma;
}

////////////////////////////////////////////////////////////////////////////////////////
//=============================================
// Functions to compute fock matrix
//=============================================
/////////////////////////////////////////////////////////////////////////////////////////

double squared_distance(const std::vector<double> &r1, const std::vector<double> &r2) {
    double sum = 0.0;
    for (int i = 0; i < 3; i++) {
        sum += std::pow(r1[i] - r2[i], 2);
    }

    return sum;
}

double two_electron_integral(double sigma_A, double sigma_B, double R2_AB) {
    // Define constants
    double V2 = 1.0 / (sigma_A + sigma_B);
    double T = V2 * R2_AB;
    double U_A = std::pow(M_PI * sigma_A, 1.5);
    double U_B = std::pow(M_PI * sigma_B, 1.5);

    // Compute the integral
    if (R2_AB < 1e-10) {    // RA == RB therefore (RA - RB) = 0
        return U_A * U_B * std::sqrt(2 * V2) * std::sqrt(2 / M_PI);

    } else {                // RA != RB

        return (U_A * U_B / std::sqrt(R2_AB)) * std::erf(std::sqrt(T));
    }
}

double compute_gamma(const BF& bf1, const BF& bf2) {
    double gamma = 0.0;

    const std::vector<double>& RA = bf1.r;
    const std::vector<double>& RB = bf2.r;
    double R2_AB = squared_distance(RA, RB);

    //Compute the sum over all four primitive gaussian combinations
    for (const auto& pg_k : bf1.pg) {
        for (const auto& pg_kp : bf1.pg) {
            double sigma_A = 1.0 / (pg_k.alpha_k + pg_kp.alpha_k);
            double d_k_norm_A = pg_k.d_k * pg_k.norm;
            double d_kp_norm_A = pg_kp.d_k * pg_kp.norm;

            for (const auto& pg_l : bf2.pg) {
                for (const auto& pg_lp : bf2.pg) {
                    double sigma_B = 1.0 / (pg_l.alpha_k + pg_lp.alpha_k);
                    double d_l_norm_B = pg_l.d_k * pg_l.norm;
                    double d_lp_norm_B = pg_lp.d_k * pg_lp.norm; 

                    // Compute the two electron integral (in eV)
                    double integral = two_electron_integral(sigma_A, sigma_B, R2_AB) * 27.211;

                    // Apply the summation this iteration.
                    gamma += d_k_norm_A * d_kp_norm_A * d_l_norm_B * d_lp_norm_B * integral;
                }
            }
        }
    }

    return gamma;
}

arma::vec derivative_two_electron_integral(
    const double& sigma_A, 
    const double& sigma_B, 
    const std::vector<double>& RA,
    const std::vector<double>& RB) {
    
    // Define constants
    arma::vec RAB(3);
    for (int i = 0; i < 3; i++) {
        RAB[i] = RA[i] - RB[i];
    }
    double R2 = squared_distance(RA, RB);              
    double R = std::sqrt(R2);                       
    double V2 = 1.0 / (sigma_A + sigma_B);
    double T = V2 * R2;
    double U_A = std::pow(M_PI * sigma_A, 1.5);
    double U_B = std::pow(M_PI * sigma_B, 1.5);

    // Handle self interactions
    if (R2 < 1e-12) {
        return arma::zeros(3);
    }
    
    arma::vec FA = U_A * U_B * RAB / R2;
    double FB = (-std::erf(std::sqrt(T)) / R) + (2 * std::sqrt(V2) * std::exp(-T) / std::sqrt(M_PI));

    return FA * FB;
}

arma::vec compute_gamma_derivative(const BF& bf1, const BF& bf2) {

    // Define the derivative
    arma::vec gamma_grad = arma::zeros(3);

    // Extract centers
    std::vector<double> RA = bf1.r;       
    std::vector<double> RB = bf2.r;       

    //Compute the sum over all four primitive gaussian combinations
    for (const auto& pg_k : bf1.pg) {
        for (const auto& pg_kp : bf1.pg) {
            double sigma_A = 1.0 / (pg_k.alpha_k + pg_kp.alpha_k);
            double d_k_norm_A = pg_k.d_k * pg_k.norm;
            double d_kp_norm_A = pg_kp.d_k * pg_kp.norm;

            for (const auto& pg_l : bf2.pg) {
                for (const auto& pg_lp : bf2.pg) {
                    double sigma_B = 1.0 / (pg_l.alpha_k + pg_lp.alpha_k);
                    double d_l_norm_B = pg_l.d_k * pg_l.norm;
                    double d_lp_norm_B = pg_lp.d_k * pg_lp.norm;

                    arma::vec gamma_dRA = derivative_two_electron_integral(sigma_A, sigma_B, RA, RB);
                    double c_norm = d_k_norm_A * d_kp_norm_A * d_l_norm_B * d_lp_norm_B;

                    // Update the derivative
                    gamma_grad += c_norm * gamma_dRA;
                    
                }
            }
        }
    }

    return gamma_grad * 27.211; //eV
}

arma::mat build_gamma_matrix(const std::vector<Atom>& atoms,
    const std::vector<BF>& basis_functions) {
    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    int N_atoms = atoms.size();
    arma::mat gamma(N_atoms, N_atoms, arma::fill::zeros);

    // For each atom, pick the first AO centered on it (the s orbital)
    std::vector<int> atom_to_ao(N_atoms, -1);
    for (size_t i = 0; i < ao_to_atom.size(); ++i) {
        int atom_idx = ao_to_atom[i];
        if (atom_to_ao[atom_idx] == -1) {
        atom_to_ao[atom_idx] = i;  // Store the first AO for that atom
        }
    }

    // Loop over all pairs of atoms A and B
    for (int A = 0; A < N_atoms; ++A) {
        for (int B = 0; B <= A; ++B) {
            const BF& bf_A = basis_functions[atom_to_ao[A]];
            const BF& bf_B = basis_functions[atom_to_ao[B]];

            double gamma_AB = compute_gamma(bf_A, bf_B);
            gamma(A, B) = gamma_AB;
            gamma(B, A) = gamma_AB;  // Symmetric
        }
    }

    return gamma;
}

arma::cube build_derivative_gamma_matrix(
    const std::vector<Atom>& atoms,
    const std::vector<BF>& basis_functions) {

    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    int N_atoms = atoms.size();
    arma::cube gamma(3, N_atoms, N_atoms, arma::fill::zeros);

    // For each atom, pick the first AO centered on it (the s orbital)
    std::vector<int> atom_to_ao(N_atoms, -1);
    for (size_t i = 0; i < ao_to_atom.size(); ++i) {
        int atom_idx = ao_to_atom[i];
        if (atom_to_ao[atom_idx] == -1) {
        atom_to_ao[atom_idx] = i;  // Store the first AO for that atom
        }
    }

    // Loop over all pairs of atoms A and B
    for (int A = 0; A < N_atoms; ++A) {
        for (int B = 0; B <= A; ++B) {
            const BF& bf_A = basis_functions[atom_to_ao[A]];
            const BF& bf_B = basis_functions[atom_to_ao[B]];

            arma::vec3 d = arma::zeros(3);
            
            d = compute_gamma_derivative(bf_A, bf_B);

            // Store anti-symmetrically
            gamma.slice(A).col(B) = d;
            gamma.slice(B).col(A) = -d;
        }
    }

    return gamma;
}

std::vector<double> compute_total_density_on_atoms(
    const arma::mat& p_tot,
    const std::vector<int> ao_to_atom,
    int num_atoms) {

    // Obtain the number of basis functions
    int N = p_tot.n_rows;

    // Initialize the total density on the atoms as zero
    std::vector<double> p_tot_atoms(num_atoms, 0.0);

    // Loop through the diagonal values
    for (int mu = 0; mu < N; mu++) {
        // Determine which AO corresponds to this index
        int atom_idx = ao_to_atom[mu];
        p_tot_atoms[atom_idx] += p_tot(mu, mu);
    }

    return p_tot_atoms;
}

arma::mat compute_fock_matrix(
    const std::vector<BF>& basis_functions,
    const std::vector<Atom>& atoms,
    const arma::mat& p_tot,
    const arma::mat& p_exp) {

    // Create a map of ionization energies plus electron affinities to AO's
    std::map<std::pair<int, char>, double> I_plus_A_map = {
    {{1, 's'}, 7.176},                           // H
    {{6, 's'}, 14.051}, {{6, 'p'}, 5.572},       // C
    {{7, 's'}, 19.316}, {{7, 'p'}, 7.275},       // N
    {{8, 's'}, 25.390}, {{8, 'p'}, 9.111},       // O
    {{9, 's'}, 32.272}, {{9, 'p'}, 11.080}       // F
    };

    // Create a map of atomic bonding parameters to AO's
    std::map<int, double> beta_map = {
        {1, -9},       // H
        {6, -21},      // C
        {7, -25},      // N
        {8, -31},      // O
        {9, -39}       // F
    };

    // Map the atomic numbers to the number of valence electrons
    std::map<int, int> valence_electrons = {
        {1, 1},     // H
        {6, 4},     // C
        {7, 5},     // N
        {8, 6},     // O
        {9, 7}      // F
    };
    
    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
                ao_to_atom.push_back(j);
                
                // Stop comparing once a match has been found (Each AO belongs to one atom)
                break;
            }
        }
    }

    // Compute the total density on the atoms
    std::vector<double> p_tot_on_atoms = compute_total_density_on_atoms(p_tot, ao_to_atom, atoms.size());

    // Compute the overlap matrix
    arma::mat S = convert_to_arma(compute_contracted_overlap_matrix(basis_functions));

    // Compute the gamma matrix in this scope
    arma::mat gamma_matrix = build_gamma_matrix(atoms, basis_functions);

    // Initialize the fock matrix
    int N = basis_functions.size();
    arma::mat F(N, N, arma::fill::zeros);

    // Loop over matrix elements to fill them using CNDO/2 equations
    for (int mu = 0; mu < N; mu++) {
        const BF& bf_mu = basis_functions[mu];
        int A = ao_to_atom[mu];
        double p_AA = p_tot_on_atoms[A];
        bool is_s_orb  = (bf_mu.l == 0 && bf_mu.m == 0 && bf_mu.n == 0);
        char orb_type_mu = is_s_orb ? 's' : 'p';
        double I_plus_A = I_plus_A_map.at({bf_mu.atomic_number, orb_type_mu});
        double gamma_AA = gamma_matrix(A, A);

        for (int nu = 0; nu < N; nu++) {
            const BF& bf_nu = basis_functions[nu];
            int B = ao_to_atom[nu];
            double p_BB = p_tot_on_atoms[B];
            double gamma_AB = gamma_matrix(A, B);

            //std::cout << "Computing fock(" << mu << ", " << nu << ")" << std::endl;

            // Fill diagonal elements
            if (mu == nu) {
                double p_exp_mu_nu = p_exp(mu, nu);
                double fock = -I_plus_A + ((p_AA - valence_electrons.at(atoms[A].atomic_number)) - (p_exp_mu_nu - 0.5)) * gamma_AA;
                
                /* This is print debugging stuff
                std::cout << "Accessing I_plus_A_map at (" << bf_mu.atomic_number << ", " << orb_type_mu << ")" << std::endl;
                std::cout << "Diagonal Term: ";
                std::cout << -I_plus_A << " + ((" << p_AA << " - " << valence_electrons.at(atoms[A].atomic_number)
                    << ") - (" << p_exp_mu_nu << " - 0.5)) * " << gamma_AA << std::endl;
                */

                // Compute summation term
                for (int B_ = 0; B_ < atoms.size(); B_++) {
                    if (B_ == A) continue;
                    double p_BB_ = p_tot_on_atoms[B_];
                    fock += (p_tot_on_atoms[B_] - valence_electrons.at(atoms[B_].atomic_number)) * gamma_matrix(A, B_);
                    
                    /* This is print debugging stuff
                    std::cout << "Summation term: ";
                    std::cout << "(" << p_tot_on_atoms[B_] << " - " << valence_electrons.at(atoms[B_].atomic_number)
                        << ") * " << gamma_matrix(A, B_) << std::endl;
                    */
                }

                // Update the element value
                F(mu, nu) = fock;

            // Fill off diagonal elements    
            } else {
                double beta_A = beta_map.at(atoms[A].atomic_number);
                double beta_B = beta_map.at(atoms[B].atomic_number);
                double s_mu_nu = S(mu, nu);
                double fock = 0.5 * (beta_A + beta_B) * s_mu_nu - p_exp(mu, nu) * gamma_AB;

                /* This is print debugging stuff
                std::cout << "Off-diagonal term: ";
                std::cout << "0.5 * (" << beta_A << " + " << beta_B << ") * " << s_mu_nu
                    << " - " << p_exp(mu, nu) << " * " << gamma_AB << std::endl;
                */

                //Update the element value
                F(mu, nu) = fock;
            }
        }
    }

    return F;
}

//////////////////////////////////////////////////////////////////////////////////////
//==============================================
// Function to compute core hamiltonian matrix
//==============================================
//////////////////////////////////////////////////////////////////////////////////////

arma::mat compute_core_ham_matrix(const std::vector<BF>& basis_functions,
    const std::vector<Atom>& atoms) {

    // Create a map of ionization energies plus electron affinities to AO's
    std::map<std::pair<int, char>, double> I_plus_A_map = {
        {{1, 's'}, 7.176},                           // H
        {{6, 's'}, 14.051}, {{6, 'p'}, 5.572},       // C
        {{7, 's'}, 19.316}, {{7, 'p'}, 7.275},       // N
        {{8, 's'}, 25.390}, {{8, 'p'}, 9.111},       // O
        {{9, 's'}, 32.272}, {{9, 'p'}, 11.080}       // F
    };

    // Create a map of atomic bonding parameters to AO's
    std::map<int, double> beta_map = {
        {1, -9},       // H
        {6, -21},      // C
        {7, -25},      // N
        {8, -31},      // O
        {9, -39}       // F
    };

    // Map the atomic numbers to the number of valence electrons
    std::map<int, int> valence_electrons = {
        {1, 1},     // H
        {6, 4},     // C
        {7, 5},     // N
        {8, 6},     // O
        {9, 7}      // F
    };

    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    // Compute the overlap matrix
    arma::mat S = convert_to_arma(compute_contracted_overlap_matrix(basis_functions));

    // Compute the gamma matrix in this scope
    arma::mat gamma_matrix = build_gamma_matrix(atoms, basis_functions);

    // Initialize the core hamiltonian matrix
    int N = basis_functions.size();
    arma::mat H(N, N, arma::fill::zeros);

    // Loop over matrix elements to fill them using CNDO/2 equations
    for (int mu = 0; mu < N; mu++) {
        const BF& bf_mu = basis_functions[mu];
        int A = ao_to_atom[mu];
        bool is_s_orb  = (bf_mu.l == 0 && bf_mu.m == 0 && bf_mu.n == 0);
        char orb_type_mu = is_s_orb ? 's' : 'p';
        double I_plus_A = I_plus_A_map.at({bf_mu.atomic_number, orb_type_mu});
        double gamma_AA = gamma_matrix(A, A);

        for (int nu = 0; nu < N; nu++) {
            const BF& bf_nu = basis_functions[nu];
            int B = ao_to_atom[nu];
            double gamma_AB = gamma_matrix(A, B);

            //std::cout << "Computing fock(" << mu << ", " << nu << ")" << std::endl;

            // Fill diagonal elements
            if (mu == nu) {
                double ham = -I_plus_A  - (valence_electrons.at(atoms[A].atomic_number) - 0.5) * gamma_AA;

                /* This is print debugging stuff
                std::cout << "Accessing I_plus_A_map at (" << bf_mu.atomic_number << ", " << orb_type_mu << ")" << std::endl;
                std::cout << "Diagonal Term: ";
                std::cout << -I_plus_A << " + ((" << p_AA << " - " << valence_electrons.at(atoms[A].atomic_number)
                    << ") - (" << p_exp_mu_nu << " - 0.5)) * " << gamma_AA << std::endl;
                */

                // Compute summation term
                for (int B_ = 0; B_ < atoms.size(); B_++) {
                    if (B_ == A) continue;
                    ham -= valence_electrons.at(atoms[B_].atomic_number) * gamma_matrix(A, B_);

                    /* This is print debugging stuff
                    std::cout << "Summation term: ";
                    std::cout << "(" << p_tot_on_atoms[B_] << " - " << valence_electrons.at(atoms[B_].atomic_number)
                        << ") * " << gamma_matrix(A, B_) << std::endl;
                    */
                }

                // Update the element value
                H(mu, nu) = ham;

            // Fill off diagonal elements    
            } else {
                double beta_A = beta_map.at(atoms[A].atomic_number);
                double beta_B = beta_map.at(atoms[B].atomic_number);
                double s_mu_nu = S(mu, nu);
                double ham = 0.5 * (beta_A + beta_B) * S(mu, nu);

                /* This is print debugging stuff
                std::cout << "Off-diagonal term: ";
                std::cout << "0.5 * (" << beta_A << " + " << beta_B << ") * " << s_mu_nu
                << " - " << p_exp(mu, nu) << " * " << gamma_AB << std::endl;
                */

                //Update the element value
                H(mu, nu) = ham;
            }
        }
    }

    return H;
}

//////////////////////////////////////////////////////////////////////////////////////
//=============================================
// Function to compute CNDO/2 total energy
//=============================================
//////////////////////////////////////////////////////////////////////////////////////

double compute_CNDO_total_energy(const arma::mat& p_alpha,
    const arma::mat& p_beta,
    const arma::mat& H,
    const arma::mat& Fa,
    const arma::mat& Fb,
    const std::vector<Atom>& atoms) {

    // Compute the electron energy
    double alpha_term = 0.0;
    double beta_term = 0.0;
    int N = p_alpha.n_rows;

    for (int mu = 0; mu < N; mu++) {
        for(int nu = 0; nu < N; nu++) {
            alpha_term += p_alpha(mu, nu) * (H(mu, nu) + Fa(mu, nu));
            beta_term += p_beta(mu, nu) * (H(mu, nu) + Fb(mu, nu));
        }
    }

    double electron_energy = 0.5 * (alpha_term + beta_term);

    // Map the atomic numbers to the number of valence electrons
    std::map<int, int> valence_electrons = {
        {1, 1},     // H
        {6, 4},     // C
        {7, 5},     // N
        {8, 6},     // O
        {9, 7}      // F
    };

    // Compute the nuclear replusion energy
    double nre = 0.0;

    for (size_t A = 0; A < atoms.size() ; A++) {
        for (size_t B = 0; B < A; B++) {

            const std::vector<double>& RA = atoms[A].coord;
            const std::vector<double>& RB = atoms[B].coord;

            double R_AB = std::sqrt(squared_distance(RA, RB));

            nre += (valence_electrons.at(atoms[A].atomic_number)
             * valence_electrons.at(atoms[B].atomic_number)) / R_AB;
        }
    }

    double total_energy = electron_energy + nre * 27.211;

    std::cout << "Nuclear Repulsion Energy is " << nre * 27.211 << " eV" << std::endl;
    std::cout << "Electron Energy is " << electron_energy << " eV" << std::endl;

    return total_energy;
}

arma::mat compute_Xuv(
    const std::vector<Atom>& atoms,
    const std::vector<BF>& basis_functions,
    const arma::mat& p_tot) {

    int N = basis_functions.size();
    
    // Create a map of atomic bonding parameters to AO's
    std::map<int, double> beta_map = {
        {1, -9},       // H
        {6, -21},      // C
        {7, -25},      // N
        {8, -31},      // O
        {9, -39}       // F
    };

    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    // Define the matrix
    arma::mat Xuv(N, N, arma::fill::zeros);
    double p_tot_sum = arma::accu(p_tot);

    // Loop over matrix elements to fill them using CNDO/2 equations
    for (int mu = 0; mu < N; mu++) {
        for (int nu = 0; nu < N; nu++) {
            if (mu == nu) continue;
            int A = ao_to_atom[mu];
            int B = ao_to_atom[nu];
            double beta_A = beta_map.at(atoms[A].atomic_number);
            double beta_B = beta_map.at(atoms[B].atomic_number);

            Xuv(mu, nu) += (beta_A + beta_B) * p_tot(mu, nu);
        }
    }

    return Xuv;
}

arma::mat compute_yAB(
    const arma::mat& p_alpha,
    const arma::mat& p_beta,
    const std::vector<Atom>& atoms,
    const std::vector<BF>& basis_functions) {

    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    // Compute the total density on the atoms
    arma::mat p_tot = p_alpha + p_beta;
    std::vector<double> p_tot_on_atoms = compute_total_density_on_atoms(p_tot, ao_to_atom, atoms.size());
    p_alpha.print("p_alpha");
    p_beta.print("p_beta");
    p_tot.print("p_tot");

    int N = basis_functions.size();

    // Map the atomic numbers to the number of valence electrons
    std::map<int, int> valence_electrons = {
        {1, 1},     // H
        {6, 4},     // C
        {7, 5},     // N
        {8, 6},     // O
        {9, 7}      // F
    };

    // Compute yAB
    arma::mat yAB(atoms.size(), atoms.size(), arma::fill::zeros);
    for (int A = 0; A < atoms.size(); A++) {
        for (int B = 0; B < atoms.size(); B++) {
            double ZA = valence_electrons.at(atoms[A].atomic_number);
            double ZB = valence_electrons.at(atoms[B].atomic_number);

            double p_AA = p_tot_on_atoms[A];
            double p_BB = p_tot_on_atoms[B];

            std::cout << "P_AA(i, j):" << p_AA << "(" << A << ", " << B << ")" << std::endl;
            std::cout << "P_BB(i, j):" << p_BB << "(" << A << ", " << B << ")" << std::endl;
            if (A != B) {
                yAB(A, B) = p_AA * p_BB - ZB * p_AA - ZA * p_BB;
    
                for (int mu = 0; mu < N; mu++) {
                    if (ao_to_atom[mu] != A) continue;
                    for (int nu = 0; nu < N; nu++) {
                        if (ao_to_atom[nu] != B) continue;
                        yAB(A, B) -= std::pow(p_alpha(mu, nu), 2) + std::pow(p_beta(mu, nu), 2);
                    }
                }
            } else {
                yAB(A, A) = std::pow(p_AA, 2) - p_AA * (ZA + ZA);
    
                for (int mu = 0; mu < N; mu++) {
                    if (ao_to_atom[mu] != A) continue;
                    for (int nu = 0; nu < N; nu++) {
                        if (ao_to_atom[nu] != B) continue;
                        yAB(A, A) -= std::pow(p_alpha(mu, nu), 2) + std::pow(p_beta(mu, nu), 2);
                    }
                }
            }
        }
    }

    return yAB;
}

arma::mat compute_grad_elec(
    const arma::mat& Xuv,
    const arma::mat& yAB,
    const arma::cube& Suv_RA_cube,
    const arma::cube& gammaAB_RA_cube,
    const std::vector<BF>& basis_functions,
    const std::vector<Atom>& atoms){

    // Map the AO to the atom from the list of basis functions
    std::vector<int> ao_to_atom;
    for (size_t i = 0; i < basis_functions.size(); i++) {
        const BF& bf = basis_functions[i];

        // Loop through all of the atoms
        for (size_t j = 0; j < atoms.size(); j++) {

            // Compare the atom coordinates with the basis function (AO) centers
            if (atoms[j].coord == bf.r) {
            ao_to_atom.push_back(j);

            // Stop comparing once a match has been found (Each AO belongs to one atom)
            break;
            }
        }
    }

    arma::mat grad_elec(3, atoms.size());
    int N = basis_functions.size();

    // Loop over the atoms
    for (int A = 0; A < atoms.size(); A++) {

        // Compute the x-term
        for (int mu = 0; mu < N; mu++) {
            for (int nu = 0; nu < N; nu++) {
                int B = ao_to_atom[nu];
                if (mu == nu) continue;

                if (ao_to_atom[mu] != A && ao_to_atom[nu] == A) continue;
                grad_elec(0, A) += Xuv(mu, nu) * Suv_RA_cube(0, mu, nu);
                grad_elec(1, A) += Xuv(mu, nu) * Suv_RA_cube(1, mu, nu);
                grad_elec(2, A) += Xuv(mu, nu) * Suv_RA_cube(2, mu, nu);
        
            }
        }

        // Compute the y-term
        for (int a = 0; a < atoms.size(); a++) {
            for (int b = 0; b < atoms.size(); b++) {
                //if (b == atom) continue;
                if (a != A && b == A) {
                    grad_elec(0, A) += yAB(a, b) * gammaAB_RA_cube(0, a, b);
                    grad_elec(1, A) += yAB(a, b) * gammaAB_RA_cube(1, a, b);
                    grad_elec(2, A) += yAB(a, b) * gammaAB_RA_cube(2, a, b);
                }
            }
        }
    }

    return grad_elec;
}

arma::mat compute_grad_nuc(const std::vector<Atom>& atoms) {

    arma::mat grad_nuc(3, atoms.size(), arma::fill::zeros);

    // Map the atomic numbers to the number of valence electrons
    std::map<int, int> valence_electrons = {
        {1, 1},     // H
        {6, 4},     // C
        {7, 5},     // N
        {8, 6},     // O
        {9, 7}      // F
    };

    for (int A = 0; A < atoms.size(); A++) {
        const std::vector<double>& RA = atoms[A].coord;
        double ZA = valence_electrons.at(atoms[A].atomic_number);

        for (int B = 0; B < atoms.size(); B++) {
            if (A == B) continue;

            const std::vector<double>& RB = atoms[B].coord;
            double ZB = valence_electrons.at(atoms[B].atomic_number);

            arma::vec3 RA_minus_RB;
            for (int i = 0; i < 3; i++) {
                RA_minus_RB(i) = RA[i] - RB[i];
            }
            
            double RAB = std::sqrt(squared_distance(RA, RB));


            grad_nuc.col(A) -= ZA * ZB * RA_minus_RB / std::pow(RAB, 3);
        }
    }

    return grad_nuc * 27.211;
}

//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////MAIN/////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
    //==============================================================================
    // Instructor provided JSON config stuff
    //==============================================================================
    
    // check that a config file is supplied
    if (argc != 2){
        std::cerr << "Usage: " << argv[0] << " path/to/config.json" << std::endl; 
        return EXIT_FAILURE; 
    }
    
    // parse the config file 
    fs::path config_file_path(argv[1]);
    if (!fs::exists(config_file_path)){
        std::cerr << "Path: " << config_file_path << " does not exist" << std::endl; 
        return EXIT_FAILURE;
    }

    std::ifstream config_file(config_file_path); 
    json config = json::parse(config_file); 

    // extract the important info from the config file
    fs::path atoms_file_path = config["atoms_file_path"];
    fs::path output_file_path = config["output_file_path"];  
    int num_alpha_electrons = config["num_alpha_electrons"];
    int num_beta_electrons = config["num_beta_electrons"];

    //==============================================================================
    // Read the molecular information
    //==============================================================================
    std::vector<Atom> atoms = read_atoms(atoms_file_path);
    int num_atoms = atoms.size();
    
    std::cout << "Number of atoms: " << atoms.size() << std::endl;
    for (const auto& atom : atoms) {
        std::cout << "Atomic number: " << atom.atomic_number << ", Coordinates ("
                  << atom.coord[0] << ", " << atom.coord[1] << ", " << atom.coord[2]
                  << ")" << std::endl;
    }
    std::cout << std::endl;
    
    //==============================================================================
    // Create the basis functions and compute the normalization constants
    //==============================================================================

    // Create basis functions
    std::vector<BF> basis_functions = create_BFS_list(atoms);
    int num_basis_functions = basis_functions.size();

    // Compute normalization constants
    compute_normalization_constants(basis_functions);
     
    //===============================================================================
    // Compute updated density matrices using self-consistent field (SCF) method
    //===============================================================================

    // Define number of alpha and beta electrons
    int p = num_alpha_electrons;
    int q = num_beta_electrons;
    int N = basis_functions.size();

    // Define p_alpha and p_beta
    arma::mat p_alpha = arma::zeros(N, N);
    arma::mat p_beta = arma::zeros(N, N);

    // Compute the total density matrix
    arma::mat p_tot = p_alpha + p_beta;

    // Compute gamma matrix
    arma::mat gamma_matrix = build_gamma_matrix(atoms, basis_functions);
    gamma_matrix.print("gamma_matrix");

    // Compute the overlap matrix
    arma::mat S = convert_to_arma(compute_contracted_overlap_matrix(basis_functions));
    S.print("S");

    // Compute core hamiltonian matrix
    arma::mat H = compute_core_ham_matrix(basis_functions, atoms);

    // Initialize the fock matrices
    arma::mat f_alpha;
    arma::mat f_beta;

    // Define simulation variables
    const int max_n_steps = 150;
    const double convergence_tolerance = 1e-6;
    int step = 0;
    bool converged = false;

    // Perform the simulation loop until convergence (or until max iterations)
    while (step < max_n_steps && !converged) {

        // Build alpha and beta fock matrices
        f_alpha = compute_fock_matrix(basis_functions, atoms, p_tot, p_alpha);
        f_beta = compute_fock_matrix(basis_functions, atoms, p_tot, p_beta);      
        
        // Solve the eigenvalue problems to obtain new MO coefficients and eigenvalues
        arma::vec e_alpha;
        arma::mat c_alpha;
        arma::eig_sym(e_alpha, c_alpha, f_alpha);
        
        arma::vec e_beta;
        arma::mat c_beta;
        arma::eig_sym(e_beta, c_beta, f_beta);
        

        // Copy the old density matrices
        arma::mat p_alpha_old = p_alpha;
        arma::mat p_beta_old = p_beta;

        // Assemble new density matrices
        p_alpha.zeros();
        p_beta.zeros();

        for (int i = 0; i < p; i++) {                           // Loop through the number of alpha electrons
            p_alpha += c_alpha.col(i) * c_alpha.col(i).t();     // Multiply c_ui and c_vi (eq 1.1)
        }
        
        std::cout << "Pa, step " << step << std::endl;
        print_arma_matrix(p_alpha);

        for (int i = 0; i < q; i++) {                           // Loop through the number of beta electrons
            p_beta += c_beta.col(i) * c_beta.col(i).t();        // Multiply c_ui and c_vi (eq 1.2)
        }
        

        p_tot.zeros();
        p_tot = p_alpha + p_beta;                               // Recompute p_tot this iteration
        

        // Check for convergence
        double d_alpha = arma::abs(p_alpha - p_alpha_old).max();
        double d_beta = arma::abs(p_beta - p_beta_old).max();

        if (std::max(d_alpha, d_beta) < convergence_tolerance) {
            converged = true;
        }

        // Increment the counter
        step++;
    }
    
    // Compute total energy, print electron energy and nuclear repulsion energy
    double total_energy = compute_CNDO_total_energy(p_alpha, p_beta, H, f_alpha, f_beta, atoms);

    //==============================================================================
    // Compute Suv_RA
    //==============================================================================

    int num_3D_dims = 3; 
    arma::cube Suv_RA_cube(num_3D_dims, num_basis_functions,  num_basis_functions);

    // Compute the partial derivative of the normalized contracted guassian (S) with respect
    // to nuclear coordinates (RA)
    Suv_RA_cube = compute_Suv_Ri(atoms, basis_functions, 1);            // 1 == A, dSuv/dRA
    arma::mat Suv_RA(3, num_basis_functions * num_basis_functions);

    // Manual flattening
    for (int mu = 0; mu < num_basis_functions; mu++) {
        for (int nu = 0; nu < num_basis_functions; nu++) {
            int flat_index = mu * num_basis_functions + nu;
    
            // Grab the (x,y,z) derivative vector for this (mu, nu) pair
            for (int d = 0; d < 3; ++d) {
                Suv_RA(d, flat_index) = Suv_RA_cube(d, mu, nu);
            }
        }
    }

    //==============================================================================
    // Compute gammaAB_RA
    //==============================================================================

    arma::cube gammaAB_RA_cube(num_3D_dims, num_atoms,  num_atoms);

    gammaAB_RA_cube = build_derivative_gamma_matrix(atoms, basis_functions);
    arma::mat gammaAB_RA(num_3D_dims, num_atoms * num_atoms);

    // Manual flattening
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = 0; j < num_atoms; ++j) {
            int flat_index = i * num_atoms + j;
            gammaAB_RA.col(flat_index) = gammaAB_RA_cube.slice(i).col(j);
        }
    }

    //==============================================================================
    // Compute electronic gradient
    //============================================================================== 

    // Compute Xuv
    arma::mat Xuv = compute_Xuv(atoms, basis_functions, p_tot);
    Xuv.print("Xuv");

    // Compute yAB
    arma::mat yAB = compute_yAB(p_alpha, p_beta, atoms, basis_functions);
    yAB.print("yAB");

    // Compute Suv_RB
    //arma::cube Suv_RB_cube = compute_Suv_Ri(atoms, basis_functions, 2);            // 2 == AB, dSuv/dRB

    // Compute electronic gradient
    arma::mat gradient_electronic(num_3D_dims, num_atoms);
    gradient_electronic = compute_grad_elec(Xuv, yAB, Suv_RA_cube, gammaAB_RA_cube, basis_functions, atoms);

    //==============================================================================
    // Compute nuclear gradient
    //==============================================================================

    arma::mat gradient_nuclear(num_3D_dims, num_atoms);
    gradient_nuclear = compute_grad_nuc(atoms);

    //==============================================================================
    // Compute total energy gradient
    //==============================================================================
    arma::mat gradient(num_3D_dims, num_atoms);
    gradient = gradient_nuclear + gradient_electronic;


    // most of the code goes here 
   



    // You do not need to modify the code below this point 

    // Set print configs
    std::cout << std::fixed << std::setprecision(4) << std::setw(8) << std::right ; 

    // inspect your answer via printing
    Suv_RA.print("Suv_RA");
    gammaAB_RA.print("gammaAB_RA");
    gradient_nuclear.print("gradient_nuclear");
    gradient_electronic.print("gradient_electronic");
    gradient.print("gradient"); 

    // check that output dir exists
    if (!fs::exists(output_file_path.parent_path())){
        fs::create_directories(output_file_path.parent_path()); 
    }
    
    // delete the file if it does exist (so that no old answers stay there by accident)
    if (fs::exists(output_file_path)){
        fs::remove(output_file_path); 
    }

    // write results to file 
    Suv_RA.save(arma::hdf5_name(output_file_path, "Suv_RA", arma::hdf5_opts::append));
    gammaAB_RA.save(arma::hdf5_name(output_file_path, "gammaAB_RA", arma::hdf5_opts::append));
    gradient_nuclear.save(arma::hdf5_name(output_file_path, "gradient_nuclear", arma::hdf5_opts::append));
    gradient_electronic.save(arma::hdf5_name(output_file_path, "gradient_electronic", arma::hdf5_opts::append));
    gradient.save(arma::hdf5_name(output_file_path, "gradient", arma::hdf5_opts::append));
    
}  