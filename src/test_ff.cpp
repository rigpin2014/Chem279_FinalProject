#include "forceField.hpp"
#include <iostream>

using namespace std;

int main() {
    Molecule mol = read_sdf("/work/molecules/water_3d.sdf");
    // Print atoms and coordinates of the molecule
    for (Atom atom : mol.atoms) {
        cout << atom.atomic_number << " " << atom.x << " " << atom.y << " " << atom.z << endl;
    }
    return 0;
}
