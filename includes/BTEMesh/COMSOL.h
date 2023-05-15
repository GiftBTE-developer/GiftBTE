//
// Created by yuehu on 2023/3/6.
//

#ifndef STATICBTESOLVER_COMSOL_H
#define STATICBTESOLVER_COMSOL_H
#include "utility/utility.h"
#include <fstream>

class COMSOL
        {
public:
    std::vector<double> nodeX;
    std::vector<double> nodeY;
    std::vector<double> nodeZ;
    std::vector<std::vector<int>> volumeElements;
    std::vector<int> volumeElementIndex;
    std::vector<std::vector<int>> boundaryElements;
    std::vector<int> boundaryElementIndex;

    COMSOL(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z);

    ~COMSOL();
        };

#endif //STATICBTESOLVER_COMSOL_H
