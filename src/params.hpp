#pragma once

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include "forceField.hpp"

using namespace std;

class MMFF_Typer {
    private:
        static map<string, int> atomTypeMap;
        static map<pair<int, int>, BondParam> bondParams;
        static map<tuple<int, int, int>, AngleParam> angleParams;
        static map<tuple<int, int, int, int>, TorsionParam> torsionParams;
        static map<int, vdwParam> vdwParams;

    public:
        static int assign_atom_type(string molecule_name, string element);
        static BondParam get_bond_param(int atom_type1, int atom_type2);
        static AngleParam get_angle_param(int atom_type1, int atom_type2, int atom_type3);
        static TorsionParam get_torsion_param(int atom_type1, int atom_type2, int atom_type3, int atom_type4);
        static vdwParam get_vdw_param(int atom_type);
};

