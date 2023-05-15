//
// Created by yuehu on 2023/3/27.
//

#ifndef STATICBTESOLVER_ONED_H
#define STATICBTESOLVER_ONED_H
#include "utility/utility.h"
#include <fstream>

class OneD
        {
        public:
            std::vector<double> nodeX;
            std::vector<double> nodeY;
            std::vector<double> nodeZ;
            std::vector<std::vector<int>> volumeElements;
            std::vector<int> volumeElementIndex;
            std::vector<std::vector<int>> boundaryElements;
            std::vector<int> boundaryElementIndex;

            OneD(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z);
            ~OneD();
            //~COMSOL();
        };

#endif //STATICBTESOLVER_ONED_H
