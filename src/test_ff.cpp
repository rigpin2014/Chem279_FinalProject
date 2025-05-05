#include "forceField.hpp"
#include <iostream>
#include "params.hpp"
#include <armadillo>
#include <fstream>

using namespace std;

int main() {
    Molecule water = read_sdf("/work/molecules/water_3d.sdf");

    ForceField ff;
    ff.add_molecule(water);

    // Get position vectors of atoms
    arma::vec o_pos = arma::vec({water.get_atom(1).x, water.get_atom(1).y, water.get_atom(1).z});
    arma::vec h1_pos = arma::vec({water.get_atom(0).x, water.get_atom(0).y, water.get_atom(0).z});
    arma::vec h2_pos = arma::vec({water.get_atom(2).x, water.get_atom(2).y, water.get_atom(2).z});

    // OH1 bond vector
    arma::vec oh1_vec = h1_pos - o_pos;
    double oh1_length = arma::norm(oh1_vec);
    arma::vec oh1_unit_vec = oh1_vec / oh1_length;

    // OH2 bond vector
    arma::vec oh2_vec = h2_pos - o_pos;
    double oh2_length = arma::norm(oh2_vec);

    // Original angle between OH1 and OH2
    double cos_orig_angle = arma::dot(oh1_unit_vec, oh2_vec) / (oh1_length * oh2_length);
    double orig_angle = acos(cos_orig_angle);

    // Create rotation matrix to rotate OH1 vector starting from OH2 vector
    arma::vec rot_axis = arma::cross(oh1_vec, oh2_vec);
    rot_axis = rot_axis / arma::norm(rot_axis);

    // Generate range of bond lengths (in Angstroms) for OH bond and calculate energy of water molecule at each bond length
    double start_bond_length = 0.7;
    double end_bond_length = 1.5;
    double step_bond_length = 0.01;
    
    vector<double> bond_energies;
    for (double bond_length = start_bond_length; bond_length < end_bond_length; bond_length += step_bond_length) {
        // Get new H1 position by scaling OH1 vector by bond_length
        arma::vec h1_pos_new = o_pos + oh1_unit_vec * bond_length;

        // Set new H1 position
        water.get_atom(0).x = h1_pos_new(0);
        water.get_atom(0).y = h1_pos_new(1);
        water.get_atom(0).z = h1_pos_new(2);

        // Calculate new energy
        double energy = ff.calculate_molecule_energy(water);
        bond_energies.push_back(energy);
    }
    // Reset H1 position
    water.get_atom(0).x = h1_pos(0);
    water.get_atom(0).y = h1_pos(1);
    water.get_atom(0).z = h1_pos(2);


    // Generate range of angles (in radians) for H-O-H angle and calculate energy of water molecule at each angle
    double start_angle = 1.4;
    double end_angle = 2.4 - 1e-10; // Avoid floating point precision issues
    double step_angle = 0.01;

    vector<double> angle_energies;
    for (double angle = start_angle; angle < end_angle; angle += step_angle) {
        // Get new H1 position by scaling OH1 vector by bond_length
        arma::vec h1_rot = oh2_vec * cos(-angle) + arma::cross(rot_axis, oh2_vec) * sin(-angle) + rot_axis * arma::dot(rot_axis, oh2_vec) * (1 - cos(-angle));
        h1_rot = ((h1_rot / arma::norm(h1_rot)) * oh1_length) + o_pos;

        // Set new H1 position
        water.get_atom(0).x = h1_rot(0);
        water.get_atom(0).y = h1_rot(1);
        water.get_atom(0).z = h1_rot(2);

        // Calculate new energy
        double energy = ff.calculate_molecule_energy(water);
        angle_energies.push_back(energy);
    }
    // Reset H1 position
    water.get_atom(0).x = h1_pos(0);
    water.get_atom(0).y = h1_pos(1);
    water.get_atom(0).z = h1_pos(2);

    // Output bond energies and angle energies to file
    ofstream outfile_bond("/work/data/bond_energies.txt");
    for (double energy : bond_energies) {
        outfile_bond << energy << endl;
    }
    outfile_bond.close();

    ofstream outfile_angle("/work/data/angle_energies.txt");
    for (double energy : angle_energies) {
        outfile_angle << energy << endl;
    }
    outfile_angle.close();
    
    return 0;
}
