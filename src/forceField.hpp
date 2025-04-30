#ifndef FORCEFIELD_HPP
#define FORCEFIELD_HPP

#include <vector>
#include <string>
#include <map>

using namespace std;

extern map<string, int> element_to_atomic_number;

class Atom {
    public:
        int atomic_number;
        double x;
        double y;
        double z;
        int atom_type;
        vector<int> bonded_atoms; // Stores indices of bonded atoms in the molecule atoms vector

        Atom();
        Atom(int atomic_number, double x, double y, double z, int atom_type);
};

class BondParam {
    public:
        double k; // Bond force constant [millidynes/angstrom]
        double r0; // Equilibrium bond length [Angstroms]

        BondParam();
        BondParam(double k, double r0);
};

class Bond {
    public:
        int atom1;
        int atom2;
        BondParam param;

        Bond();
        Bond(int atom1, int atom2, BondParam param);
};

class AngleParam {
    public:
        double k; // Angle force constant [millidynes * angstrom/radian]
        double theta0; // Equilibrium bond angle [radians]

        AngleParam();
        AngleParam(double k, double theta0);
};

class Angle {
    public:
        int atom1;
        int atom2;
        int atom3;
        AngleParam param;

        Angle();
        Angle(int atom1, int atom2, int atom3, AngleParam param);
};

class Molecule {
    public:
        string name;
        vector<Atom> atoms;
        vector<Bond> bonds;

        Molecule();
        Molecule(string name, vector<Atom> atoms, vector<Bond> bonds);

        void add_atom(int atomic_number, double x, double y, double z, int atom_type);
        void add_bond(int atom1, int atom2, BondParam param);

        vector<Angle> get_angles();
};

class ForceField {
    public:
        vector<Molecule> molecules;

        ForceField();
        ForceField(vector<Molecule> molecules);

        void add_molecule(Molecule mol);

        double calculate_total_energy();
        // double calculate_total_bond_energy();
        // double calculate_total_angle_energy();
        // double calculate_total_torsion_energy();
        // double calculate_total_nonbonded_energy();

        double calculate_molecule_energy(Molecule mol);
        double calculate_bond_energy(Molecule mol);
        double calculate_angle_energy(Molecule mol);
        // double calculate_torsion_energy(Molecule mol);
        // double calculate_nonbonded_energy(Molecule mol);
};

class TorsionParam {
    public:
        double k;
        double n;
        double delta;

        TorsionParam(double k, double n, double delta);
};

class NonBondedParam {
    public:
        double epsilon;
        double sigma;

        NonBondedParam(double epsilon, double sigma);
};

// Function to read a sdf file and return a Molecule object
Molecule read_sdf(string filename);
void assign_atom_types(Molecule& mol);

#endif // FORCEFIELD_HPP


