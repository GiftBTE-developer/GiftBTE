//
// Created by yuehu on 2021/9/10.
//

#ifndef STATICBTESOLVER_BTEANGLE_H
#define STATICBTESOLVER_BTEANGLE_H

#include "utility/utility.h"

class BTEAngle {   
public:
    double SolidAngle;
    std::vector<vec> direction;
    std::vector<double> weight;
    std::vector<vec> control_direction;

    BTEAngle() = default;
    BTEAngle(int ntheta, int nphi, int Dimension_material, int Angle_method);
    BTEAngle(const BTEAngle &angle) = delete;
    // void copyVector();
    ~BTEAngle();
};

#endif // STATICBTESOLVER_BTEANGLE_H
