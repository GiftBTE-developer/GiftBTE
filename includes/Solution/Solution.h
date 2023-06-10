//
// Created by yuehu on 2023/3/6.
//

#ifndef STATICBTESOLVER_SOLUTION_H
#define STATICBTESOLVER_SOLUTION_H
#include "utility/utility.h"
#include "StaticBTESynthetic/StaticBTESynthetic.h"
#include "StaticFourier/StaticFourier.h"
#include "TransientBTE/transient.h"
#include "BTEMesh/Distribute.h"

class Solution{
public:
    std::string solverType;
    int numCell;
    int numBandLocal;
    int numDirectionLocal;
    int numBound;
    int numDirection;
    int numBand;
    double ***energyDensity;
    double *ebound;
    double *ebound_old;
    double *Temperature;
    double *Temperature_old;
    double *heatFluxX;
    double *heatFluxY;
    double *heatFluxZ;
    double *heatFluxX_old;
    double *heatFluxY_old;
    double *heatFluxZ_old;
    double *boundaryTemp;
    double *boundaryFlux;
    StaticBTESynthetic *BTEsolver;
    StaticFourier *Fouriersolver;
};

class SolutionAll{
public:
    Solution solution;
    std::vector<Solution> solutions;
    SolutionAll(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles,int num_proc, int world_rank);

    void _Fourier_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc, int world_rank) const;

    void _BTE_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank
                     ,int Use_Backup, int Num_Max_Iter, string Order, string Method, string Matrix_solver, double error_temp_limit, double error_flux_limit) const;
    void _Transient_BTE_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank,
                               int Use_Backup, string Order, double error_temp_limit, double error_flux_limit, double deltaT, double totalT);

    void _set_initial(DistributeMesh *mesh,BTEBand *bands, BTEAngle *angles);
    void _update_bound(int is, DistributeMesh *mesh,BTEBand *bands, BTEAngle *angles);
    void _print_out(DistributeMesh *mesh);
    void _get_solver(int is, DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank,int Use_Backup, int Num_Max_Iter, int Use_Sythetic, int Use_Limiter, double error_temp_limit, double error_flux_limit);
    void _solve(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank,int Use_Backup, int Num_Max_Iter, int Use_Sythetic, int Use_Limiter, double error_temp_limit, double error_flux_limit);
    void _transfer_solution(DistributeMesh *mesh);

    //void _get_convergent();
    ~SolutionAll();
};

#endif //STATICBTESOLVER_SOLUTION_H
