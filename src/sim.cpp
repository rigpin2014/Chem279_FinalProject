#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>      // std::isstringstream
#include <vector>       // std::vector
#include <string>       // std::string
#include <regex>        // std::regex_replace
#include <algorithm>    // std::sort

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

int main(int argc, char** argv){

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
}