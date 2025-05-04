#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>      // std::isstringstream
#include <vector>       // std::vector
#include <string>       // std::string
#include <regex>        // std::regex_replace
#include <algorithm>    // std::sort
#include <random>       // For generating random numbers
#include <cmath>        // std::sqrt, std::legendre
#include <numeric>      // std::inner_product
#include <nlohmann/json.hpp>    // This is the JSON handling library

// convenience definitions so the code is more readable
namespace fs = std::filesystem;
using json = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////////////////
//================================
// Define some custom data types
//================================
////////////////////////////////////////////////////////////////////////////////////////////

// Define the atom
struct Atom {
    int atomic_number;
    std::vector<double> coord;

    Atom(int atomic_number, double x, double y, double z)
        : atomic_number(atomic_number), coord{x, y, z} {}
};

class Molecule {
    public:

        // Constructor
        Molecule(const std::vector<Atom>& input_atoms) : atoms(input_atoms) {
            assert(atoms.size() == 3);              // Make sure this only works for three atoms systems
            assert(atoms[0].atomic_number == 8);    // Make sure the first atom is Oxygen
            update_geometry();
            update_decomposition_status();
        }

        std::vector<double> compute_bond_energies() const {
            std::vector<double> bond_energies(bond_lengths.size());
            for (int i = 0; i < bond_lengths.size(); i++) {
                bond_energies[i] = 0.5 * k_bond * std::pow(bond_lengths[i] - r_eq, 2);
            }
            return bond_energies;
        }

        double compute_bond_angle_energy() const {
            return 0.5 * k_angle * std::pow(bond_angle - theta_eq, 2);
        }

        void displace_atom(int atom_index, double recoil_energy, double step = 0.0001) {
            
            while (!decomposed) {
                atoms[atom_index].coord[0] += step;
                update_geometry();
                update_decomposition_status();

                // Get the energies
                std::vector<double> bond_energies = compute_bond_energies();
                double bond_angle_energy = compute_bond_angle_energy();

                // Define the evaluated energy
                double evaluated_energy;
                if (atoms[atom_index].atomic_number == 8) {     // If Oxygen is the incident atom
                    evaluated_energy = bond_energies[0] + bond_energies[1] + bond_angle_energy;
                } else {                                        // If Hydrogen is the incident atom
                    evaluated_energy = bond_energies[0] + bond_angle_energy;
                }

                // Stop iterating once all of the imparted energy has been accounted for
                if (evaluated_energy >= recoil_energy) break;
            }
        }

        bool get_decomposed() const {return decomposed;}

    private:
        // Define force field variables
        std::vector<Atom> atoms;
        std::vector<double> bond_lengths;
        double bond_angle;
        bool decomposed = false;

        // Define force field constants
        const double k_bond = 450.0;        // Bond stiffness (kcal/mol/A^2)
        const double k_angle = 55.0;        // Angle stiffness (kcal/mol/rad^2)
        const double r_eq = 0.9572;         // Equilibrium OH bond length (A)
        const double theta_eq = 1.824;      // Equilibrium HOH bond angle (rad)

        // Decomposition thresholds
        const double OH_bond_critical_energy = 117.5;    // kcal/mol
        const double max_bond_angle = 2.27;              // rad

        // Private member functions
        void update_geometry() {
            bond_lengths.clear();
            for (int i = 1; i < atoms.size(); i++) {
                bond_lengths.push_back(compute_distance(atoms[0], atoms[i]));
            }
            bond_angle = compute_bond_angle();
        }

        double compute_distance(const Atom& a1, const Atom& a2) {
            double sum = 0.0;
            for (int i = 0; i < 3; i++){
                sum += std::pow(a1.coord[i] - a2.coord[i], 2);
            }
            return std::sqrt(sum);
        }

        double compute_bond_angle() const {
            // Obtain vectors in the direction of the bonds
            std::vector<double> vec1(3), vec2(3);
            for (int i = 0; i < 3; i++) {
                vec1[i] = atoms[1].coord[i] - atoms[0].coord[i];
                vec2[i] = atoms[2].coord[i] - atoms[0].coord[i];
            }

            // Compute the dot product and the norms
            double dot = std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0);
            double norm1 = std::sqrt(std::inner_product(vec1.begin(), vec1.end(), vec1.begin(), 0.0));
            double norm2 = std::sqrt(std::inner_product(vec2.begin(), vec2.end(), vec2.begin(), 0.0));

            // Return the bond angle
            return std::acos(dot / (norm1 * norm2));
        }

        void update_decomposition_status() {
            std::vector<double> bond_energies = compute_bond_energies();

            // Update the flag if the decomposition criteria are met
            if (bond_energies[0] >= OH_bond_critical_energy ||
                bond_energies[1] >= OH_bond_critical_energy ||
                bond_angle >= max_bond_angle) {
                    decomposed = true;
            } 
        }
};

//////////////////////////////////////////////////////////////////////////////////
//=======================================
// Read in the coordinates of each atom
//=======================================
//////////////////////////////////////////////////////////////////////////////////

std::vector<Atom> read_atoms(const std::string& filename) {
    std::ifstream file (filename);

    int num_atoms;
    file >> num_atoms;

    std::string line;
    std::getline(file, line);
    std::getline(file, line);       // Skip the comment line

    int atomic_number;
    double x, y, z;
    std::vector<Atom> atoms;

    while (file >> atomic_number >> x >> y >> z) {
        atoms.push_back(Atom(atomic_number, x, y, z));
    }

    return atoms;
}


// Create a function to convert the numbers to scientific notation
std::string convert_to_sci(
    const std::string& s
){
    // Insert a literal E between the mantissa and the exponent.
    // Insert a literal space after the exponent to fix formatting issue in raw data.
    std::string E = std::regex_replace(s, std::regex("([0-9])([-+][0-9]{1,2})"), "$1E$2 ");
    return E;
}

// Create a function to read angular distributions based on the number of Legendre Coefficients
void read_angular_distributions(
    const std::string& file,
    std::vector<std::pair<double, std::vector<double>>>& data
){
    // Open the file
    std::ifstream infile(file);

    std::string energy_num_coeffs_line;
    while (std::getline(infile, energy_num_coeffs_line)) {          // Move the pointer to the next line
        std::istringstream iss1(energy_num_coeffs_line);
        std::string _, NeutronEnergyRaw;
        int __, ___, num_coeffs;

        // Read the line
        iss1 >> _ >> NeutronEnergyRaw;
        NeutronEnergyRaw = convert_to_sci(NeutronEnergyRaw);
        double NeutronEnergy = std::stod(NeutronEnergyRaw);

        iss1 >> __ >> ___ >> num_coeffs;

        // Define the number of coefficients on line one and line two
        int num_coeffs_line1, num_coeffs_line2;
        if (num_coeffs < 6){
            num_coeffs_line1 = num_coeffs;
            num_coeffs_line2 = 0;
        } else {
            num_coeffs_line1 = 6;
            num_coeffs_line2 = num_coeffs - 6;
        }

        // Define the coefficient vector
        std::vector<double> Legendre_Coefficients;

        // Read the first coefficient line
        std::string coeff_line1;
        if (!std::getline(infile, coeff_line1)){                    // Move the pointer to the next line
            std::cerr << "End of file" << std::endl;
            return;
        }
        coeff_line1 = convert_to_sci(coeff_line1);
        std::istringstream iss2(coeff_line1);

        // Append the first line of coefficients
        for (int i = 0; i < num_coeffs_line1; i++) {
            std::string token;
            iss2 >> token;
            Legendre_Coefficients.push_back(std::stod(token));
        }

        // Read the second coefficient line if there is one
        if (num_coeffs > 6) {
            std::string coeff_line2;
            if (!std::getline(infile, coeff_line2)){                // Move the pointer to the next line
                std::cerr << "End of file" << std::endl;
                return;
            }
            coeff_line2 = convert_to_sci(coeff_line2);
            std::istringstream iss3(coeff_line2);

            // Append the second line of coefficients
            for (int i = 0; i < num_coeffs_line2; i++) {
                std::string token;
                iss3 >> token;
                Legendre_Coefficients.push_back(std::stod(token));
            }

        }
        
        // Append the pair to the data vector
        data.push_back(std::make_pair(NeutronEnergy, Legendre_Coefficients));
        
    }
}

// Define a helper function to determine if a value being read needs formatting
bool needs_formatting(
    const std::string& s
){
    // Check if the string contains '+' or '-' indicating it has an exponent
    for (size_t i = 1; i < s.size(); i++) {
        if (s[i] == '+' || s[i] == '-'){
            return true;
        }
    }

    return false;
}

// Define a function to read the microscopic cross section data
void read_cross_sections(
    const std::string& file,
    std::vector<std::pair<double, double>>& data
){
    // Open the file
    std::ifstream infile(file);

    // Read the file one line at a time
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss1(line);

        // Loop through the three pairs of data on this line
        for (int i = 0; i < 3; i++){
            std::string NeutronEnergyRaw, CrossSectionRaw;
            iss1 >> NeutronEnergyRaw >> CrossSectionRaw;

            // Check if the number needs to be formatted
            if (needs_formatting(NeutronEnergyRaw)) {
                NeutronEnergyRaw = convert_to_sci(NeutronEnergyRaw);
            }
            
            if (needs_formatting(CrossSectionRaw)) {
                CrossSectionRaw = convert_to_sci(CrossSectionRaw);
            }

            // Convert strings to doubles
            double NeutronEnergy = std::stod(NeutronEnergyRaw);
            double CrossSection = std::stod(CrossSectionRaw);

            // Append the data
            data.push_back(std::make_pair(NeutronEnergy, CrossSection));

            // Sort the data by NeutronEnergy in ascending order
            std::sort(
                data.begin(),
                data.end(),
                [](const auto &a, const auto &b){
                    return a.first < b.first;
                }
            );
        }
    }
}

// Define a function for the prompt fission spectrum
double chi(double E) {
    return 0.453 * std::exp(-1.036 * E) * std::sinh(std::sqrt(2.29 * E));
}

// Define a function to compute the cumulative distribution function (CDF)
void compute_chi_cdf(
    const std::vector<std::pair<double, double>>& energy_map,
    std::vector<double>& energies,
    std::vector<double>& cdf
){
    // Size the vectors
    size_t N = energy_map.size();
    energies.resize(N);
    cdf.resize(N);

    // Set the initial values of the vectors
    cdf[0] = 0.0;
    energies[0] = energy_map[0].first;

    // Compute the CDF
    for (int i = 1; i < N; i++) {
        double E_prev = energy_map[i-1].first;
        double E_curr = energy_map[i].first;
        double chi_prev = chi(E_prev);
        double chi_curr = chi(E_curr);
        double dE = E_curr - E_prev;

        cdf[i] = cdf[i-1] + 0.5 * (chi_prev + chi_curr) * dE;
        energies[i] = E_curr;
    }

    // Normalize the CDF between [0, 1]
    double max = cdf.back();
    for (double& val : cdf) {
        val /= max;
    }
}

// Define a function to sample the CDF to obtain a neutron energy randomly.
double sample_from_cdf(
    const std::vector<double>& data,
    const std::vector<double>& cdf
){
    // Generate a random number uniformly between zero and one.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double r = dist(gen);

    // Define an iterator where cdf[i] >= r
    auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
    if (it == cdf.begin()) return data.first();
    if (it == cdf.end()) return data.back();

    // Find the index of the iterator
    size_t i = std::distance(cdf.begin(), it);

    // Define interpolation variables
    double cdf_low = cdf[i-1], cdf_high = cdf[i];
    double data_low = data[i-1], data_high = data[i];
    double f = (r - cdf_low) / (cdf_high = cdf_low);

    // Interpolate
    return data_low + f * (data_high - data_low);
}

// Define a function to extract the Legendre polynomials given the neutron energy
std::vector<double> extract_polynomials(
    const double& sampled_energy,
    const std::vector<std::pair<double, std::vector<double>>>& angular_distributions
){
    // Compute the minimum difference between the sampled energy and the table of energies
    double  min = std::abs(angular_distributions[0].first - sampled_energy);
    size_t min_idx = 0;

    // Loop through the angular distributions and find the index with the minimum difference
    for (int i = 1; i < angular_distributions.size(); i++) {
        double diff = std::abs(angular_distributions[i].first - sampled_energy);
        if (diff < min) {
            min = diff;
            min_idx = i;
        }
    }

    return angular_distributions[min_idx].second;
}

// Define a function to compute the probability density function over the angular distributions
double compute_angular_pdf(
    double mu,
    const std::vector<double>& legendre_coeffs
){
    double result = 0.5 * std::legendre(0, mu);

    for (int l = 1; l < legendre_coeffs.size(); l++){
        result += (2 * l + 1) * 0.5 * legendre_coeffs[l - 1] * std::legendre(l, mu);
    }

    return result;
}

// Define a function to compute the most probable scattering angle
double compute_most_probable_scattering_angle(
    const std::vector<double>& legendre_coeffs
){
    // Create a differential domain for mu on [-1, 1]
    int N = 1000;
    std::vector<double> mu_domain(N);
    double dmu = 2 / (N - 1);

    // Define the maximum values
    double pdf_max = -1.0;
    double mu_max = -1.0;

    // Determine the most probable scattering angle by iterating through the
    // probability distribution function
    for (int i = 0; i < N; i++) {
        double mu = -1.0 + i * dmu;
        double pdf_curr = compute_angular_pdf(mu, legendre_coeffs);

        // Update the maximum value 
        if (pdf_curr > pdf_max) {
            pdf_max = pdf_curr;
            mu_max = mu;
        }
    }

    return std::acos(mu_max);
}


int main(int argc, char** argv){

    //==============================================================================
    // Instructor provided JSON config stuff
    //==============================================================================
    
    // check that a config file is supplied
    if (argc != 2){
        std::cerr << "Usage: " << argv[0] << " path/to/config.json" << std::endl; 
        return EXIT_FAILURE; 
    }
    
    // parse the config file 
    fs::path config_file_path(argv[1]);
    if (!fs::exists(config_file_path)){
        std::cerr << "Path: " << config_file_path << " does not exist" << std::endl; 
        return EXIT_FAILURE;
    }

    std::ifstream config_file(config_file_path); 
    json config = json::parse(config_file); 

    // extract the important info from the config file
    fs::path atoms_file_path = config["atoms_file_path"];
    fs::path output_file_path = config["output_file_path"];  
    int num_alpha_electrons = config["num_alpha_electrons"];
    int num_beta_electrons = config["num_beta_electrons"];

    //===============================================================================
    // Read the tables
    //===============================================================================

    // Read in the H-1 angular distributions
    std::vector<std::pair<double, std::vector<double>>> H1_angular_distributions;
    read_angular_distributions("../../tables/H1_angular_distributions.txt", H1_angular_distributions);

    // Verify an entry
    std::cout << "H-1 Angular distribution check: " << std::endl;
    std::cout << "Neutron Energy: " << H1_angular_distributions[0].first << std::endl;
    std::cout << "Legendre Coefficients: " << std::endl;
    for (const auto& coeff: H1_angular_distributions[0].second) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl << std::endl;

    // Read the O-16 angular distributions
    std::vector<std::pair<double, std::vector<double>>> O16_angular_distributions;
    read_angular_distributions("../../tables/O16_angular_distributions.txt", O16_angular_distributions);

    // Verify an entry
    std::cout << "O-16 Angular distribution check: " << std::endl;
    std::cout << "Neutron Energy: " << O16_angular_distributions[7].first << std::endl;
    std::cout << "Legendre Coefficients: " << std::endl;
    for (const auto& coeff: O16_angular_distributions[7].second) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl << std::endl;

    // Read the H-1 cross sections
    std::vector<std::pair<double, double>> H1_cross_sections;
    read_cross_sections("../../tables/H1_cross_sections.txt", H1_cross_sections);

    // Verify an entry
    std::cout << "H-1 Cross section check: " << std::endl;
    std::cout << "Neutron Energy: " << H1_cross_sections[0].first << std::endl;
    std::cout << "Cross Section: " << H1_cross_sections[0].second << std::endl << std::endl;

    //Read the O-16 cross sections
    std::vector<std::pair<double, double>> O16_cross_sections;
    read_cross_sections("../../tables/O16_cross_sections.txt", O16_cross_sections);

    // Verify an entry
    std::cout << "O-16 Cross section check: " << std::endl;
    std::cout << "Neutron Energy: " << O16_cross_sections[0].first << std::endl;
    std::cout << "Cross Section: " << O16_cross_sections[0].second << std::endl << std::endl;

    //===============================================================================
    // MC Simulation
    //===============================================================================

    // Define the neutron energies and the CDF
    std::vector<double> energies;
    std::vector<double> cdf;
    compute_chi_cdf(O16_cross_sections, energies, cdf);         // Modify in-place

    // Obtain a neutron energy from the prompt fission spectra
    double neutron_energy = sample_from_cdf(energies, cdf);

    // Extract the Legendre coefficients for the sampled neutron energy.
    std::vector<double> legendre_coeffs = extract_polynomials(neutron_energy, H1_angular_distributions);

    // Compute the most probable scattering angle in radians
    double most_probable_theta = compute_most_probable_scattering_angle(legendre_coeffs);
    return 0;
}