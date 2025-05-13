#include "forceField.hpp"
#include "params.hpp"
#include <fstream>
#include <regex>
#include <cmath>
#include <iostream>
#include <sstream>
#include <armadillo>
#include <limits>

using namespace std;

// Conversion factor from mdynes/A to kcal/mol
double MDYNE_PER_A_TO_KCAL_PER_MOL = 143.9325;

// Map of element to atomic number
map<string, int> element_to_atomic_number = {
    {"H", 1},
    {"C", 6},
    {"N", 7},
    {"O", 8},
};

// ATOM CLASS FUNCTIONS
Atom::Atom(int atomic_number, double x, double y, double z, int atom_type)
    : atomic_number(atomic_number), x(x), y(y), z(z), atom_type(atom_type), bonded_atoms() {}

// Calculate bond length between two atoms
double atom_dist(const Atom& a1, const Atom& a2) {
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Calculate angle between three atoms (in degrees)
double angle_between(const Atom& a1, const Atom& a2, const Atom& a3) {
    arma::vec v1 = arma::vec({a1.x - a2.x, a1.y - a2.y, a1.z - a2.z});
    arma::vec v2 = arma::vec({a3.x - a2.x, a3.y - a2.y, a3.z - a2.z});
    double dot = arma::dot(v1, v2);
    double len1 = arma::norm(v1);
    double len2 = arma::norm(v2);
    double radians = acos(dot / (len1 * len2));
    double degrees = radians * 180.0 / M_PI;
    return degrees;
}

// BOND CLASS FUNCTIONS
BondParam::BondParam()
    : k(0.0), r0(0.0) {}

BondParam::BondParam(double k, double r0)
    : k(k), r0(r0) {}

Bond::Bond(int atom1, int atom2, BondParam param)
    : atom1(atom1), atom2(atom2), param(param) {}

// ANGLE CLASS FUNCTIONS
AngleParam::AngleParam()
    : k(0.0), theta0(0.0) {}

AngleParam::AngleParam(double k, double theta0)
    : k(k), theta0(theta0) {}

Angle::Angle(int atom1, int atom2, int atom3, AngleParam param)
    : atom1(atom1), atom2(atom2), atom3(atom3), param(param) {}

// TORSION CLASS FUNCTIONS
TorsionParam::TorsionParam()
    : v1(0.0), v2(0.0), v3(0.0) {}

TorsionParam::TorsionParam(double v1, double v2, double v3)
    : v1(v1), v2(v2), v3(v3) {}

Torsion::Torsion(int atom1, int atom2, int atom3, int atom4, TorsionParam param)
    : atom1(atom1), atom2(atom2), atom3(atom3), atom4(atom4), param(param) {}

// VDW CLASS FUNCTIONS
vdwParam::vdwParam()
    : alpha(0.0), N(0.0), A(0.0), G(0.0), type('X') {}

vdwParam::vdwParam(double alpha, double N, double A, double G, char type)
    : alpha(alpha), N(N), A(A), G(G), type(type) {}

// MOLECULE CLASS FUNCTIONS
Molecule::Molecule()
    : atoms(), bonds() {}

Molecule::Molecule(string name, vector<Atom> atoms, vector<Bond> bonds)
    : name(name), atoms(atoms), bonds(bonds) {}

void Molecule::add_atom(int atomic_number, double x, double y, double z, int atom_type) {
    atoms.push_back(Atom(atomic_number, x, y, z, atom_type));
}

void Molecule::add_bond(int atom1, int atom2, BondParam param) {
    bonds.push_back(Bond(atom1, atom2, param));
    atoms[atom1].bonded_atoms.push_back(atom2);
    atoms[atom2].bonded_atoms.push_back(atom1);
}

// Function to get the largest necessary dimension of a square/cube that can contain the molecule
double Molecule::get_volume_dimension() {
    double positive_infinity = std::numeric_limits<double>::infinity();
    double negative_infinity = -std::numeric_limits<double>::infinity();

    double min_x = positive_infinity;
    double min_y = positive_infinity;
    double min_z = positive_infinity;

    double max_x = negative_infinity;
    double max_y = negative_infinity;
    double max_z = negative_infinity;

    for (Atom atom : atoms) {
        if (atom.x < min_x) min_x = atom.x;
        if (atom.y < min_y) min_y = atom.y;
        if (atom.z < min_z) min_z = atom.z;
        if (atom.x > max_x) max_x = atom.x;
        if (atom.y > max_y) max_y = atom.y;
        if (atom.z > max_z) max_z = atom.z;
    }

    return max(max_x - min_x, max(max_y - min_y, max_z - min_z));
}

vector<Angle> Molecule::get_angles() {
    vector<Angle> angles;
    for (int i = 0; i < atoms.size(); i++) {
        Atom a1 = atoms[i];
        for (int a2 : a1.bonded_atoms) {
            for (int a3 : atoms[a2].bonded_atoms) {
                if (a3 == i) continue; // Not true angle between three atoms
                if (i < a3) { // Only add each angle once
                    // Type the angle
                    AngleParam angle_param = MMFF_Typer::get_angle_param(a1.atom_type, atoms[a2].atom_type, atoms[a3].atom_type);
                    angles.push_back(Angle(i, a2, a3, angle_param));
                }
            }
        }
    }
    return angles;
}

vector<Torsion> Molecule::get_torsions() {
    vector<Torsion> torsions;
    for (int i = 0; i < atoms.size(); i++) {
        Atom a1 = atoms[i];
        for (int a2 : a1.bonded_atoms) {
            for (int a3 : atoms[a2].bonded_atoms) {
                if (a3 == i) continue; // Not true torsion between four atoms
                for (int a4 : atoms[a3].bonded_atoms) {
                    if (a4 == i || a4 == a2) continue;
                    if (i < a4) { // Only add each torsion once
                        // Type the torsion
                        TorsionParam torsion_param = MMFF_Typer::get_torsion_param(a1.atom_type, atoms[a2].atom_type, atoms[a3].atom_type, atoms[a4].atom_type);
                        torsions.push_back(Torsion(i, a2, a3, a4, torsion_param));
                    }
                }
            }
        }
    }
    return torsions;
}

Atom& Molecule::get_atom(int index) {
    return atoms.at(index);
}

// FORCE FIELD CLASS FUNCTIONS
ForceField::ForceField()
    : molecules() {}

ForceField::ForceField(vector<Molecule> molecules)
    : molecules(molecules) {}

void ForceField::add_molecule(Molecule mol) {
    molecules.push_back(mol);
}

double ForceField::calculate_molecule_energy(Molecule mol) {
    double energy = 0.0;
    energy += calculate_bond_energy(mol);
    energy += calculate_angle_energy(mol);
    energy += calculate_torsion_energy(mol);
    energy += calculate_vdw_energy();
    return energy;
}

double ForceField::calculate_total_energy() {
    double energy = 0.0;
    for (Molecule mol : molecules) {
        energy += calculate_molecule_energy(mol);
    }
    return energy;
}

double ForceField::calculate_bond_energy(Molecule mol) {
    double energy = 0.0;
    for (Bond bond : mol.bonds) {
        double r = atom_dist(mol.atoms[bond.atom1], mol.atoms[bond.atom2]);
        double delta_r = r - bond.param.r0;
        energy += 0.5 * MDYNE_PER_A_TO_KCAL_PER_MOL * bond.param.k * pow(delta_r, 2) * (1 - (2.0 * delta_r) + (7.0/12.0 * (pow(-2, 2) * pow(delta_r, 2)))); // Quartic Morse potential
    }
    return energy;
}

double ForceField::calculate_angle_energy(Molecule mol) {
    double energy = 0.0;
    vector<Angle> angles = mol.get_angles();
    for (Angle angle : angles) {
        double theta = angle_between(mol.atoms[angle.atom1], mol.atoms[angle.atom2], mol.atoms[angle.atom3]);
        double delta_theta = theta - angle.param.theta0;
        double c1 = MDYNE_PER_A_TO_KCAL_PER_MOL * pow(M_PI / 180.0, 2);
        energy += 0.5 * c1 * angle.param.k * pow(delta_theta, 2) * (1.0 + (-0.007 * delta_theta)); // Harmonic potential with cubic bend expansion
    }
    return energy;
}

 // Compute torsion angle (angle between two normal vectors for the two planes defined by ijk and jkl atoms)
double compute_cos_torsion_angle(Torsion torsion, Molecule mol) {
    Atom i = mol.get_atom(torsion.atom1);
    Atom j = mol.get_atom(torsion.atom2);
    Atom k = mol.get_atom(torsion.atom3);
    Atom l = mol.get_atom(torsion.atom4);

    // Compute all bond vectors
    arma::vec r_ji = arma::vec({i.x - j.x, i.y - j.y, i.z - j.z});
    arma::vec r_jk = arma::vec({k.x - j.x, k.y - j.y, k.z - j.z});
    arma::vec r_kj = arma::vec({j.x - k.x, j.y - k.y, j.z - k.z});
    arma::vec r_kl = arma::vec({l.x - k.x, l.y - k.y, l.z - k.z});

    // Compute normal vectors for the two planes
    arma::vec n_ijk = arma::cross(r_ji, r_jk);
    arma::vec n_jkl = arma::cross(r_kj, r_kl);
    double len1 = arma::norm(n_ijk);
    double len2 = arma::norm(n_jkl);

    if (len1 == 0 || len2 == 0) return 0.0;

    // Compute torsion angle
    double dot = arma::dot(n_ijk, n_jkl);
    double cos_phi = dot / (len1 * len2);
    return cos_phi;
}

double ForceField::calculate_torsion_energy(Molecule mol) {
    double energy = 0.0;
    vector<Torsion> torsions = mol.get_torsions();
    for (Torsion torsion : torsions) {
        double cos_phi = compute_cos_torsion_angle(torsion, mol);
        double cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
        double cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);
        double v1 = torsion.param.v1;
        double v2 = torsion.param.v2;
        double v3 = torsion.param.v3;
        energy += 0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2_phi) + v3 * (1.0 - cos3_phi));
    }
    return energy;
}

double vdwParam::calculate_R_star_ii(int atom_type1) {
    vdwParam param = MMFF_Typer::get_vdw_param(atom_type1);
    return param.A * pow(param.alpha, 0.25);
}

double vdwParam::calculate_R_star_ij(int atom_type1, int atom_type2) {
    double R_star_ii = calculate_R_star_ii(atom_type1);
    double R_star_jj = calculate_R_star_ii(atom_type2);
    double gamma_ij = (R_star_ii - R_star_jj) / (R_star_ii * R_star_jj);
    return 0.5 * (R_star_ii + R_star_jj) * (1.0 + 0.2 * (1 - exp(-12 * pow(gamma_ij, 2))));
}

double vdwParam::calculate_epsilon_ij(int atom_type1, int atom_type2) {
    vdwParam param1 = MMFF_Typer::get_vdw_param(atom_type1);
    vdwParam param2 = MMFF_Typer::get_vdw_param(atom_type2);
    double R_star_ij = calculate_R_star_ij(atom_type1, atom_type2);
    return ((181.16 * param1.G * param2.G * param1.alpha * param2.alpha) / (sqrt(param1.alpha / param1.N) + sqrt(param2.alpha / param2.N))) * (1 / pow(R_star_ij, 6));
}

double ForceField::calculate_vdw_energy_pair_atoms(int atom_type1, int atom_type2, double r) {
    double energy = 0.0;
    double R_star_ij = vdwParam::calculate_R_star_ij(atom_type1, atom_type2);
    double epsilon_ij = vdwParam::calculate_epsilon_ij(atom_type1, atom_type2);
    return epsilon_ij * pow((1.07 * R_star_ij / (r + 0.07 * R_star_ij)), 7) * ((1.12 * pow(R_star_ij, 7)) / (pow(r, 7) + pow(R_star_ij, 7)) - 2);
}

bool Molecule::is_three_bond_away(int atom1, int atom2) {
    if (atom1 == atom2) return false;
    for (int a2 : atoms[atom1].bonded_atoms) {
        if (a2 == atom2) return false;
        for (int a3 : atoms[a2].bonded_atoms) {
            if (a3 == atom1) continue;
            if (a3 == atom2) return false;
            for (int a4 : atoms[a3].bonded_atoms) {
                if (a4 == atom2) return false;
            }
        }
    }
    return true;
}

double ForceField::calculate_vdw_energy() {
    double energy = 0.0;
    // Need to consider all pairs of atoms between molecules
    // and pairs of atoms within the same molecule that are separated by 3 or more bonds
    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules.size(); j++) {
            for (int a1 = 0; a1 < molecules[i].atoms.size(); a1++) {
                for (int a2 = 0; a2 < molecules[j].atoms.size(); a2++) {
                    if (i == j && !molecules[i].is_three_bond_away(a1, a2)) continue;
                    double r = atom_dist(molecules[i].atoms[a1], molecules[j].atoms[a2]);
                    energy += calculate_vdw_energy_pair_atoms(molecules[i].atoms[a1].atom_type, molecules[j].atoms[a2].atom_type, r);
                }
            }
        }
    }
    return energy;
}

Molecule read_sdf(const string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    string line;
    istringstream iss;
    Molecule mol = Molecule();
    string name;

    // Read molecule name
    getline(file, line);
    iss.str(line);
    iss >> name;
    mol.name = name;

    // Read header (next 2 lines)
    getline(file, line); // Info line
    getline(file, line); // Copyright line

    // Read number of atoms and bonds
    getline(file, line);
    iss.str(line);
    int num_atoms, num_bonds;
    iss >> num_atoms >> num_bonds;
    
    // Read atoms
    for (int i = 0; i < num_atoms; ++i) {
        getline(file, line);
        double x, y, z;
        string element;
        iss.str(line);
        iss >> x >> y >> z >> element;
        element = regex_replace(element, regex("^\\s+|\\s+$"), "");
        int atomic_number = element_to_atomic_number[element];
        int atom_type = MMFF_Typer::assign_atom_type(name, element);
        
        mol.add_atom(atomic_number, x, y, z, atom_type);
    }

    // Read bonds
    for (int i = 0; i < num_bonds; ++i) {
        getline(file, line);
        int atom1, atom2;
        iss.str(line);
        iss >> atom1 >> atom2;
        // Convert to 0-indexed
        atom1--;
        atom2--;
        BondParam param = MMFF_Typer::get_bond_param(mol.atoms[atom1].atom_type, mol.atoms[atom2].atom_type);

        mol.add_bond(atom1, atom2, param);
    }

    return mol;
}