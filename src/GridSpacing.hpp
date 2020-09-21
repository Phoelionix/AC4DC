#include <vector>
#include <iostream>
#ifndef GRIDSPACING_CXX_H
#define GRIDSPACING_CXX_H

struct GridSpacing {
    GridSpacing(){
        mode= linear;
    }
    GridSpacing(char c){
        mode= c;
    }
    const static char linear = 0;
    const static char quadratic = 1;
    const static char exponential = 2;
    const static char hybrid = 3;
    const static char unknown = 101;
    char mode = 5;
    size_t num_exp = 0; // Used for hybrid spec only
};

namespace{
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
            os << "hybrid exponential-linear grid";
            break;
        default:
            os << "Unknown grid type";
            break;
        }
        return os;
    }

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
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            gs.mode = GridSpacing::linear;
            break;
        }
        return is;
    }
}

#endif