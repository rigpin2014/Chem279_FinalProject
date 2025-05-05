#include "forceField.hpp"
#include "params.hpp"
#include <fstream>
#include <regex>
#include <cmath>
#include <iostream>
#include <sstream>
#include <armadillo>

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
// TorsionParam::TorsionParam(double k, double n, double delta)
//     : k(k), n(n), delta(delta) {}

// NONBONDED CLASS FUNCTIONS
// NonBondedParam::NonBondedParam(double epsilon, double sigma)
//     : epsilon(epsilon), sigma(sigma) {}

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

vector<Angle> Molecule::get_angles() {
    vector<Angle> angles;
    for (int i = 0; i < atoms.size(); i++) {
        Atom a1 = atoms[i];
        for (int a2 : a1.bonded_atoms) {
            for (int a3 : atoms[a2].bonded_atoms) {
                if (a3 == i) { // Not true angle between three atoms
                    continue;
                }
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

double ForceField::calculate_total_energy() {
    double energy = 0.0;
    for (Molecule mol : molecules) {
        energy += calculate_molecule_energy(mol);
    }
    return energy;
}

double ForceField::calculate_molecule_energy(Molecule mol) {
    double energy = 0.0;
    energy += calculate_bond_energy(mol);
    energy += calculate_angle_energy(mol);
    // energy += calculate_torsion_energy(mol);
    // energy += calculate_nonbonded_energy(mol);
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