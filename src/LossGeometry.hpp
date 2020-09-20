#include <vector>
#include <iostream>
#ifndef LOSS_GEOMETRY_CXX_H
#define LOSS_GEOMETRY_CXX_H

// Overkill, to be sure

struct LossGeometry {
    const static char sphere = 0;
    const static char cylinder = 1;
    const static char plane = 2;
    int mode;
    double L0;
    double factor(){
        return _C[mode]/L0;
    }
    private:
    const double _C[3] = {2, 3, 4};
};

namespace{
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, LossGeometry g) {
        switch (g.mode)
        {
        case LossGeometry::sphere:
            os << "Sphere";
            break;
        case LossGeometry::cylinder:
            os << "Thin Cylinder";
            break;
        case LossGeometry::plane:
            os << "Thin Plane";
            break;
        default:
            os << "Unknown shape";
            break;
        }
        return os;
    }

    [[maybe_unused]] std::istream& operator>>(std::istream& is, LossGeometry& lg) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No geometry specifier provided, defaulting to spherical..."<<std::endl;
            lg.mode = LossGeometry::sphere;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 's':
            lg.mode = LossGeometry::sphere;
            break;
        case 'c':
            lg.mode = LossGeometry::cylinder;
            break;
        case 'p':
            lg.mode = LossGeometry::plane;
            break;
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            lg.mode = LossGeometry::sphere;
            break;
        }
        return is;
    }
}

#endif