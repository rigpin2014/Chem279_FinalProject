#include "forceField.hpp"
#include <iostream>
#include "params.hpp"
#include <armadillo>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
    // Check that an sdf file was provided
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " path/to/sdf/file" << endl;
        return 1;
    }

    // Read input file
    Molecule mol = read_sdf(argv[1]);

    // Compute energy
    ForceField ff;
    ff.add_molecule(mol);
    double energy = ff.calculate_total_energy();
    cout << "Energy for " << mol.name << " is " << energy << " kcal/mol" << endl;
}