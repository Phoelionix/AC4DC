#include <vector>
#include <iostream>
#include <cassert>
#ifndef LOSS_GEOMETRY_CXX_H
#define LOSS_GEOMETRY_CXX_H

// Overkill, to be sure

struct LossGeometry {
    const static int sphere = 0;
    const static int cylinder = 1;
    const static int plane = 2;
    int mode = 0;
    double L0 = 1;
    double factor() const {
        assert(mode >=0 && mode <=2);
        return 1.*constants[mode]/L0;
    }
    private:
    const double constants[3] = {2, 3, 4};
};

namespace{
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, const LossGeometry& g) {
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