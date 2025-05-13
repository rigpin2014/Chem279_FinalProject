# Chem279_FinalProject

To run our project first clone the project. Then from the root directory, launch the docker container through:
`./interactive.sh`
Then to build our C++ code run:
`./clean.sh`
`./build.sh`
To run the simulation code run:
`./run.sh`

The MMFF94 implementation code can also be experimented with. To output the data that I collected from varying the bond lengths and angles of a water molecule, you can run:
`./run_water_plot_data.sh`
This will write the energy data to the files `data/angle_energies.txt` and `data/bond_energies.txt`. Once those are populated, you can run `./install_python.sh` to install the necessary Python packages and then run `python rdkit_comp.py` to see how I compared the RDKit implementation of MMFF94 to my own. This is how I generated the comparison plots also located in the `data` directory.

Additionally, you can simply calculate the energy of any of the molecules that are in the `molecules` directory, by running `./run_mol_energy.sh`. Based on which molecule you want to get the energy for, you can change the file path in the `run_mol_energy.sh` file.

Note: If any scripts don't work, try `chmod +x /path_to_script.sh` to make sure they are executable first.
