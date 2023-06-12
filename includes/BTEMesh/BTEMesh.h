//
// Created by yuehu on 2021/9/8.
//

#ifndef STATICBTESOLVER_BTEMESH_H
#define STATICBTESOLVER_BTEMESH_H
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <cassert>

#include "utility/utility.h"
#include "BTEBand/BTEBand.h"
#include "BTEBoundaryCondition/BTEBoundaryCondition.h"


class BTEMesh {
public:
    int nedge=2;
    int Dimension=1;
    double L_x;
    double L_y;
    double L_z;
    std::vector<Cell> Elements;
    std::vector<Point> Nodes;
    std::vector<Point> Heatnodes; //yufei:add Heatnodes here
    std::vector<Point> Tempnodes; //jiaxuan: for initial temperature
    std::vector<Boundary> Boundaries;
    std::vector<int> Boundnodes;
    std::vector<std::vector<int>> Boundnodes_cell;
    std::vector<int> Boundnodesconnect;
    std::vector<int> Boundnodesindex;
    std::vector<faces> Element_Faces;
    std::vector<std::vector<int>> Boundaryindex;

    BTEMesh() = default;
    BTEMesh(int Dimension_Geometry,std::ifstream& inFile, double L_x, double L_y , double L_z,BTEBand * bands, BTEBoundaryCondition *bcs,std::string mesh_type); 
    //BTEMesh(const BTEMesh& mesh) = delete;
    BTEMesh(int Dimension_Geometry,double L_x,double L_y,double L_z, std::vector<double> &nodeX,std::vector<double> &nodeY, std::vector<double> &nodeZ, std::vector<std::vector<int>> &volumeElements
            , std::vector<int> &volumeElementIndex, std::vector<std::vector<int>> &boundaryElements, std::vector<int> &boundaryElementIndex) ;


    void setMeshParams(BTEBand * bands);
    void setMeshParams1(BTEBoundaryCondition *bcs);

    void BTEMesh_heatin( std::ifstream &inHeat,double Uniform_Heat, std::string heat_type ); //yufei
    void BTEMesh_initialTemp( std::ifstream &initialTemp,double Tref); //jiaxuan
    ~BTEMesh();
};
#endif //STATICBTESOLVER_BTEMESH_H
