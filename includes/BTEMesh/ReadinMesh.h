//
// Created by yuehu on 2023/4/1.
//

#ifndef STATICBTESOLVER_READINMESH_H
#define STATICBTESOLVER_READINMESH_H
#include "utility/utility.h"
#include <fstream>
class ReadIn
        {
public:
            std::vector<double> nodeX;
            std::vector<double> nodeY;
            std::vector<double> nodeZ;
            std::vector<std::vector<int>> volumeElements;
            std::vector<int> volumeElementIndex;
            std::vector<std::vector<int>> boundaryElements;
            std::vector<int> boundaryElementIndex;
            ReadIn(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z, std::string mesh_type);
            void OneD(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z);
            void COMSOL(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z);
            void MSH(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z);
            ~ReadIn();
        };
#endif //STATICBTESOLVER_READINMESH_H
