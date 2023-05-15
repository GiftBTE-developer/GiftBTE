//
// Created by huyue on 2023/1/17.
//

#ifndef STATICBTESOLVER_TRANSIENTFOURIER_H
#define STATICBTESOLVER_TRANSIENTFOURIER_H
#include "BTEMesh/BTEMesh.h"
#include "BTEBand/BTEBand.h"
#include "BTEBoundaryCondition/BTEBoundaryCondition.h"
#include "BTEAngle/BTEAngle.h"
class TransientFourier
{
public:
    double deltaT;
    double totaltime;
    double *capacity;

    int Dimension;
    double *temperaturen;
    double *temperaturen_1;
    double *temperatureTemp;
    double *temperatureVertex;
    int numcell;
    int numBound;
    double *kappa;

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

    int *boundaryType;
    double *boundaryTemp;

    unsigned int *csrRowPtr;
    unsigned int * csrColInd;
    double *csrVal;
    double *Re;

    TransientFourier(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc, int world_rank);
    void _set_initial();
    void _get_coefficient();
    void _set_vertex_temperature();
    void _get_Re();
    void _print_out();
    void _solve();
};

#endif //STATICBTESOLVER_TRANSIENTFOURIER_H
