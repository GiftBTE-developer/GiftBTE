//
// Created by yuehu on 2023/3/6.
//

#ifndef STATICBTESOLVER_DISTRIBUTE_H
#define STATICBTESOLVER_DISTRIBUTE_H
#include "utility/utility.h"
#include "BTEMesh/BTEMesh.h"
#include "ReadinMesh.h"
using namespace std;
class DistributeMesh{
public:
    std::vector<double> nodeXFourier;
    std::vector<double> nodeYFourier;
    std::vector<double> nodeZFourier;
    std::vector<std::vector<int>> volumeElementsFourier;
    std::vector<int> volumeElementIndexFourier;
    std::vector<std::vector<int>> boundaryElementsFourier;
    std::vector<int> boundaryElementIndexFourier;

    std::vector<double> nodeXBTE;
    std::vector<double> nodeYBTE;
    std::vector<double> nodeZBTE;
    std::vector<std::vector<int>> volumeElementsBTE;
    std::vector<int> volumeElementIndexBTE;
    std::vector<std::vector<int>> boundaryElementsBTE;
    std::vector<int> boundaryElementIndexBTE;
    std::vector<int> nodeProjection;
    std::vector<int> elementProjection;
    std::vector<int> boundaryProjection;




    BTEMesh FourierMeshes;
    BTEMesh BTEMeshes;
    std::vector<BTEMesh> mesh;

    DistributeMesh(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z,string mesh_type,
                   BTEBand * bands, BTEBoundaryCondition *bcs, std::ifstream &inHeat,double Uniform_Heat, string heat_type); //yufei
    void _build_BTEMesh(int Dimension_Geometry,double L_x, double L_y, double L_z,BTEBand * bands, BTEBoundaryCondition *bcs, std::ifstream &inHeat,double Uniform_Heat,std::string mesh_solver_file);
    ~DistributeMesh();
};
#endif //STATICBTESOLVER_DISTRIBUTE_H
