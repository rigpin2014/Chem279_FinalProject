#include "forceField.hpp"
#include "params.hpp"
#include <fstream>
#include <regex>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

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

// Calculate angle between three atoms (in radians)
double angle_between(const Atom& a1, const Atom& a2, const Atom& a3) {
    double dx1 = a1.x - a2.x;
    double dy1 = a1.y - a2.y;
    double dz1 = a1.z - a2.z;

    double dx2 = a3.x - a2.x;
    double dy2 = a3.y - a2.y;
    double dz2 = a3.z - a2.z;

    // Normalize vectors
    double len1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
    double len2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
    dx1 /= len1;
    dy1 /= len1;
    dz1 /= len1;
    dx2 /= len2;
    dy2 /= len2;
    dz2 /= len2;

    double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
    return acos(dot);
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
TorsionParam::TorsionParam(double k, double n, double delta)
    : k(k), n(n), delta(delta) {}

// NONBONDED CLASS FUNCTIONS
NonBondedParam::NonBondedParam(double epsilon, double sigma)
    : epsilon(epsilon), sigma(sigma) {}

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
                if (a3 == i) {
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
        energy += 0.5 * bond.param.k * pow(r - bond.param.r0, 2);
    }
    return energy;
}

double ForceField::calculate_angle_energy(Molecule mol) {
    double energy = 0.0;
    vector<Angle> angles = mol.get_angles();
    for (Angle angle : angles) {
        double theta = angle_between(mol.atoms[angle.atom1], mol.atoms[angle.atom2], mol.atoms[angle.atom3]);
        energy += 0.5 * angle.param.k * pow(theta - angle.param.theta0, 2);
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

    // Read header (next 3 lines)
    getline(file, line); // Molecule name
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