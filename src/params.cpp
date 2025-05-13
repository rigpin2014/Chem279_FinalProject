#include "params.hpp"
#include <stdexcept>
#include <algorithm>
/*
Citations:
    Thomas A. Halgren, J. Comput. Chem., 17, 490-519 (1996).
    Thomas A. Halgren, J. Comput. Chem., 17, 520-552 (1996).
    Thomas A. Halgren, J. Comput. Chem., 17, 553-586 (1996).
    Thomas A. Halgren and Robert B. Nachbar, J. Comput. Chem., 17, 587-615 (1996).
    Thomas A. Halgren, J. Comput. Chem., 17, 616-641 (1996).
 */

map<string, int> MMFF_Typer::atomTypeMap = {
    {"OH2", 70}, // Oxygen in water
    {"HOH", 31}, // Hydrogen in water
    {"CR", 1}, // Allylic carbon
    {"C=C", 2}, // Vinyllic carbon
    {"HC", 5}, // Hydrogen attached to carbon
};

map<pair<int, int>, BondParam> MMFF_Typer::bondParams = {
    {make_pair(31, 70), BondParam(7.880, 0.969)}, // O-H bond in water
    {make_pair(1, 5), BondParam(4.776, 1.093)}, // C-H bond
    {make_pair(1, 1), BondParam(4.258, 1.508)}, // C-C single bond in alkane
    {make_pair(2, 2), BondParam(9.505, 1.333)}, // C=C double bond in alkene
    {make_pair(2, 5), BondParam(5.170, 1.083)}, // C-H bond with vinylic carbon
};

map<tuple<int, int, int>, AngleParam> MMFF_Typer::angleParams = {
    {make_tuple(31, 70, 31), AngleParam(0.658, 103.978)}, // O-H-O angle in water
    {make_tuple(5, 1, 5), AngleParam(0.516, 108.836)}, // H-C-H angle with allylic carbon
    {make_tuple(1, 1, 5), AngleParam(0.636, 110.549)}, // C-C-H angle with allylic carbons
    {make_tuple(5, 2, 5), AngleParam(0.365, 119.523)}, // H-C-H angle with vinylic carbon
    {make_tuple(2, 2, 5), AngleParam(0.535, 121.004)}, // C=C-H angle with vinylic carbon
};

map<tuple<int, int, int, int>, TorsionParam> MMFF_Typer::torsionParams = {
    {make_tuple(5, 1, 1, 5), TorsionParam(0.284, -1.386, 0.314)}, // H-C-C-H torsion with allylic carbons
    {make_tuple(5, 2, 2, 5), TorsionParam(-0.406, 1.767, 0.000)}, // H-C-C-H torsion with vinylic carbons
};

map<int, vdwParam> MMFF_Typer::vdwParams = {
    // Based on atom types.
    {1, vdwParam(1.050, 2.490, 3.890, 1.282, 'X')}, 
    {2, vdwParam(1.350, 2.490, 3.890, 1.282, 'X')}, 
    {5, vdwParam(0.250, 0.800, 4.200, 1.209, 'X')}, 
    {31, vdwParam(0.150, 0.800, 4.200, 1.209, 'D')},
    {70, vdwParam(0.870, 3.150, 3.890, 1.282, 'A')},
};

int MMFF_Typer::assign_atom_type(string molecule_name, string element) {
    transform(molecule_name.begin(), molecule_name.end(), molecule_name.begin(), ::tolower);

    if (molecule_name == "water") {
        if (element == "H") {
            return atomTypeMap["HOH"];
        } else if (element == "O") {
            return atomTypeMap["OH2"];
        }
    } else if (molecule_name == "ethylene") {
        if (element == "H") {
            return atomTypeMap["HC"];
        } else if (element == "C") {
            return atomTypeMap["C=C"];
        }
    } else if (molecule_name == "methane") {
        if (element == "H") {
            return atomTypeMap["HC"];
        } else if (element == "C") {
            return atomTypeMap["CR"];
        }
    } else if (molecule_name == "ethane") {
        if (element == "H") {
            return atomTypeMap["HC"];
        } else if (element == "C") {
            return atomTypeMap["CR"];
        }
    }
    return 0;
}

BondParam MMFF_Typer::get_bond_param(int atom_type1, int atom_type2) {
    if (bondParams.find(make_pair(atom_type1, atom_type2)) != bondParams.end()) {
        return bondParams[make_pair(atom_type1, atom_type2)];
    } else if (bondParams.find(make_pair(atom_type2, atom_type1)) != bondParams.end()) {
        return bondParams[make_pair(atom_type2, atom_type1)];
    } else {
        throw invalid_argument("Bond parameter not found for atom types " + to_string(atom_type1) + " and " + to_string(atom_type2));
    }
}

AngleParam MMFF_Typer::get_angle_param(int atom_type1, int atom_type2, int atom_type3) {
    if (angleParams.find(make_tuple(atom_type1, atom_type2, atom_type3)) != angleParams.end()) {
        return angleParams[make_tuple(atom_type1, atom_type2, atom_type3)];
    } else if (angleParams.find(make_tuple(atom_type3, atom_type2, atom_type1)) != angleParams.end()) {
        return angleParams[make_tuple(atom_type3, atom_type2, atom_type1)];
    } else {
        throw invalid_argument("Angle parameter not found for atom types " + to_string(atom_type1) + ", " + to_string(atom_type2) + ", and " + to_string(atom_type3));
    }
}

TorsionParam MMFF_Typer::get_torsion_param(int atom_type1, int atom_type2, int atom_type3, int atom_type4) {
    if (torsionParams.contains(make_tuple(atom_type1, atom_type2, atom_type3, atom_type4))) {
        return torsionParams[make_tuple(atom_type1, atom_type2, atom_type3, atom_type4)];
    } else if (torsionParams.contains(make_tuple(atom_type4, atom_type3, atom_type2, atom_type1))) {
        return torsionParams[make_tuple(atom_type4, atom_type3, atom_type2, atom_type1)];
    } else {
        throw invalid_argument("Torsion parameter not found for atom types " + to_string(atom_type1) + ", " + to_string(atom_type2) + ", " + to_string(atom_type3) + ", and " + to_string(atom_type4));
    }
}

vdwParam MMFF_Typer::get_vdw_param(int atom_type) {
    if (vdwParams.contains(atom_type)) {
        return vdwParams[atom_type];
    } else {
        throw invalid_argument("VDW parameter not found for atom type " + to_string(atom_type));
    }
}
