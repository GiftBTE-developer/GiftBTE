//
// Created by huyue on 2023/1/17.
//

#ifndef STATICBTESOLVER_STATICFOURIER_H
#define STATICBTESOLVER_STATICFOURIER_H
//
// Created by yuehu on 2022/12/12.
//

#include "BTEMesh/BTEMesh.h"
#include "BTEBand/BTEBand.h"
#include "BTEBoundaryCondition/BTEBoundaryCondition.h"
#include "BTEAngle/BTEAngle.h"

#include "Eigen/SparseLU"
#include "Eigen/Sparse"

typedef Eigen::Triplet<double> Tri;
class StaticFourier
{
public:
    int numProc;
    int worldRank;
    double L_x;
    double L_y;
    double L_z;

    double *temperature;
    double *temperatureOld;
    double *temperatureVertex;
    int numCell;
    int numBound;

    double *kappa;
    int dimension;
    int *elementFaceSize;
    double *elementCenterX;
    double *elementCenterY;
    double  *elementCenterZ;
    double *elementFaceCenterX;
    double *elementFaceCenterY;
    double *elementFaceCenterZ;
    int *elementFaceNeighobr;
    int *elementFaceBound;
    double *elementHeatSource;
    double *elementVolume;
    double *elementFaceArea;
    int *elementFaceVertexesSize;
    int *elementFaceVertexes;
    double *elementFaceNormX;
    double *elementFaceNormY;
    double *elementFaceNormZ;

    int *elementVertexesSize;
    int *elementVertexes;



    double *nodeX;
    double *nodeY;
    double *nodeZ;
    int numnode;

    double *heatFluxXGlobal;
    double *heatFluxYGlobal;
    double *heatFluxZGlobal;

    int *boundaryType;
    double *boundaryTemp;
    double *boundaryThermal;
    double *boundaryResis;
    double *boundaryFlux;
    int *boundaryCell;
    int *boundaryFace;
    int *boundaryConnect;
    std::vector<std::vector<int>> elementNeighborList;

    Eigen::SparseMatrix<double> stiffMatrix;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> cg_macro;
    Eigen::VectorXd Re1;
    Eigen::VectorXd sol1;

    unsigned int *csrRowPtr;
    unsigned int * csrColInd;
    double *csrVal;
    double *Re;

    double ***CellMatrix;

    StaticFourier(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc, int world_rank);

    void _set_initial();
    void _set_vertex_temperature();
    void _get_coefficient();
    void _get_Re();
    void _get_bound_temp() const;
    void _get_CellMatrix() const;
    void _get_heat_flux() const;
    void _get_bound_flux();
    void _print_out();
    void _solve();

    ~StaticFourier();

    void _get_coefficient_macro();
};


#endif //STATICBTESOLVER_STATICFOURIER_H
