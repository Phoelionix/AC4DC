/**
 * @file GridSpacing.hpp
 * 
 */
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#include <vector>
#include <iostream>
#ifndef GRIDSPACING_CXX_H
#define GRIDSPACING_CXX_H

struct GridSpacing {
    // GridSpacing(){
    //     mode= linear;
    // }
    // GridSpacing(char c){
    //     mode= c;
    // }
    const static char linear = 0;
    const static char quadratic = 1;
    const static char exponential = 2;
    const static char hybrid = 3;
    const static char powerlaw = 4;
    const static char test = 5;
    const static char hybrid_power = 6;
    const static char unknown = 101;
    char mode = unknown;
    size_t num_low = 0; // Used for hybrid spec only
    double transition_e = 0; // Also used to estimate Coulomb logarithm. THis should be split into a separate class TODO
    unsigned zero_degree_0 = 0; // Number of derivatives to set to zero at the origin
    unsigned zero_degree_inf = 0; // Number of derivatives to set to zero at infinity
    double min_coulomb_density=0; // Density below which to ignotre Coulomb interactions 
    // (NOTE: This should have its own class, along with coulomb estimation cutoff)
};

// currently for hybrid_power only.
struct GridBoundaries {
    // Quick and dirty attachment parameters TODO this should replace num_low, and should fill the role of 
    // transition_e so it can be split off (as the coulomb log cutoff) along with the coulomb density.
    // Should also replace num_elec_points, min_elec_e, max_elec_e
    // Note this current implementation includes the energy/index of last grid point.
    std::vector<int> start ={-1};
    bool start_parsed = false;
    std::vector<double> E_min ={-1.};  // In eV
};

namespace {
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, GridSpacing gs) {
        switch (gs.mode)
        {
        case GridSpacing::linear:
            os << "linear";
            break;
        case GridSpacing::quadratic:
            os << "quadratic";
            break;
        case GridSpacing::exponential:
            os << "exponential";
            break;
        case GridSpacing::hybrid:
            os << "hybrid linear-linear grid, transition at M="<<gs.num_low<<", e="<<gs.transition_e;
            break;
        case GridSpacing::powerlaw:
            os << "power-law grid going through M="<<gs.num_low<<", e="<<gs.transition_e;
            break;
        case GridSpacing::test:
            os << "experimental grid going through M="<<gs.num_low<<", e="<<gs.transition_e;
            break;
        default:
            os << "Unknown grid type";
            break;
        }
        return os;
    }

    /// Sets mode of the grid style aka grid type
    [[maybe_unused]] std::istream& operator>>(std::istream& is, GridSpacing& gs) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No grid type provided, defaulting to linear..."<<std::endl;
            gs.mode = GridSpacing::linear;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'l':
            gs.mode = GridSpacing::linear;
            break;
        case 'q':
            gs.mode = GridSpacing::quadratic;
            break;
        case 'e':
            gs.mode = GridSpacing::exponential;
            break;
        case 'h':
            gs.mode = GridSpacing::hybrid;
            break;
        case 'p':
            gs.mode = GridSpacing::powerlaw;
            break;
        case 't':
            gs.mode = GridSpacing::test;
            break;
        case 'w':
            gs.mode = GridSpacing::hybrid_power;
            break;
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            gs.mode = GridSpacing::linear;
            break;
        }
        return is;
    }

    /// Sets region boundaries
    [[maybe_unused]] std::istream& operator>>(std::istream& is, GridBoundaries& gb) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No boundaries provided, default case not implemented yet...!"<<std::endl;
            return is;
        }

        // Get vector from string. I'm not proud of this.
        std::vector<double> target_vector;
        std::string delimiter_char = "|";
        std::string skip_char = " ";
        size_t pos = 0;
        //size_t space_pos = 0;
        double token;
        std::string token_str;  
        while ((pos = tmp.find(delimiter_char)) != std::string::npos) {
            token_str = tmp.substr(0, pos);
            // // get rid of spaces 
            // std::string::iterator end_pos = std::remove(token_str.begin(), token_str.end(), ' ');
            // token_str.erase(end_pos, token_str.end());     
            // space_pos = token_str.find(' ');
            // token_str = tmp.substr(0, space_pos);

            // Add num to vector
            token =  stod(token_str);
            target_vector.push_back(token);
            // Erase input string up to next delimiter
            tmp.erase(0, pos + delimiter_char.length());
        }        
        if(tmp.size() > 0){
            std::string::iterator end_pos = std::remove(tmp.begin(), tmp.end(), ' ');
            tmp.erase(end_pos, tmp.end());
            target_vector.push_back(stod(tmp));
        }

        if(!gb.start_parsed){
            gb.start.resize(target_vector.size());
            for(int i = 0; i < target_vector.size(); i++){
                target_vector[i] += 0.5;
                gb.start[i] = (int)target_vector[i];
                gb.start_parsed = true;
            }
        }
        else{
            gb.E_min = target_vector;
        }
        std::cout << "elem of grid region inputs:" << std::endl;
        std::cout << gb.E_min[0] << std::endl; 
        std::cout << gb.start[0] << std::endl;
        return is;
    }    
}

#endif