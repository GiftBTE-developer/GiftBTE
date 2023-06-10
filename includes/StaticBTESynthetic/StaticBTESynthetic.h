//
// Created by yuehu on 2021/9/8.
//

#ifndef STATICBTESOLVER_STATICBTESYNTHETIC_H
#define STATICBTESOLVER_STATICBTESYNTHETIC_H
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>

#include "BTEMesh/BTEMesh.h"
#include "BTEBand/BTEBand.h"
#include "BTEBoundaryCondition/BTEBoundaryCondition.h"
#include "BTEAngle/BTEAngle.h"
#include "utility/utility.h"

#include "Eigen/SparseLU"
#include "Eigen/Sparse"

typedef Eigen::Triplet<double> Tri;


class StaticBTESynthetic
{
public:

    double L_x;
    double L_y;
    double L_z;
    //BTEMesh *mesh;
    //BTEBoundaryCondition *bcs;
    //BTEBand *bands;
    //BTEAngle *angles;

    int dimension;
    int numProc;
    int worldRank;
    int numDirection;
    int numCell;
    int numBand;
    int numEdge;
    int numDirectionLocal;
    int numBandLocal;
    int numCellLocal;
    int *numCellLocalList;
    int *startCellList;
    int startCell;
    int numNode;
    int numBound;
    int numFace;

    int * matter;

    double ***energyDensity;
    double *energyDensityVertex;
    double *eboundLocal;
    double *ebound;
    double *temperature;
    double *temperatureOld;
    double *temperatureLocal;
    double *temperatureVertex;
    double *temperature1;

    double *totalEnergy;
    double *totalEnergyLocal;
    double *heatFluxXLocal;
    double *heatFluxYLocal;
    double *heatFluxZLocal;
    double *heatFluxXGlobal;
    double *heatFluxYGlobal;
    double *heatFluxZGlobal;
    double *heatFluxXOld;
    double *heatFluxYOld;
    double *heatFluxZOld;
    double *gradientX;
    double *gradientY;
    double *gradientZ;
    double *gradientFaceX;
    double *gradientFaceY;
    double *gradientFaceZ;
    double *limit;

    double **boundaryTempMacro;



    double ***FaceMatrix;
    double ***CellMatrix;

    double error;
    double error_heat;
    double errorOld;
    double errorOld_heat;
    double error_old;
    int errorIncreaseTime;
    int intererrorIncreaseTime;

    int *elementIndex;
    int *elementFaceSize;
    double *elementCenterX;
    double *elementCenterY;
    double  *elementCenterZ;
    double *elementFaceCenterX;
    double *elementFaceCenterY;
    double *elementFaceCenterZ;
    int *elementFaceNeighobr;
    std::vector<std::vector<int>> elementNeighborList;
    int *elementFaceBound;
    double *elementHeatSource;
    double *elementVolume;
    double *elementFaceArea;
    int *elementFaceVertexesSize;
    int *elementFaceVertexes;
    double *elementFaceNormX;
    double *elementFaceNormY;
    double *elementFaceNormZ;
    int *elementFaceIndex;
    int *elementVertexesSize;
    int *elementVertexes;
    double *nodeX;
    double *nodeY;
    double *nodeZ;
    int *boundaryType;
    double *boundaryTemp;
    double *boundaryThermal;
    double *boundaryFlux;
    int **boundaryCell;
    int **boundaryFace;

    std::vector<std::vector<int>> boundaryNeighbors;

    //int *boundaryIndex;
    int *boundaryConnect;
    double ***boundaryTrans;
    int *faceBound;
    int **faceCell;
    int **faceFace;
    int numofMatter;
    double *kappaBulk;
    double *totalWeight;

    double ***groupVelocityX;
    double ***groupVelocityY;
    double ***groupVelocityZ;
    double ***relaxationTime;
    double ***heatCapacity;
    double ***latticeRatio;
    double ***heatRatio;
    double ***modeWeight;
    double *capacityBulk;

    double *directionX;
    double *directionY;
    double *directionZ;

    Eigen::SparseMatrix<double> stiffMatrix;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> **cg;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> cg_macro;
    Eigen::VectorXd Re1;
    Eigen::VectorXd sol1;

    double *Re;
    double *ReMacro;
    double *ReMacroLocal;

    unsigned int ***csrRowPtr;
    unsigned int *** csrColInd;
    double ***csrVal;

    unsigned int *csrRowPtr_macro;
    unsigned int * csrColInd_macro;
    double *csrVal_macro;

    StaticBTESynthetic(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank);

    void _set_matrix(const std::string& Matrix_type);

    void _delete_matrix(const std::string& Matrix_type) const;

    void _set_initial(int Use_Backup) const;

    void _get_coefficient(int iband_local, int inf_local);
    void _get_coefficient_macro();

    void _get_coefficient_Iterative(int iband_local, int inf_local);
    void _get_coefficient_macro_Iterative() const;

    void _get_Re(int iband_local, int inf_local);
    void _get_Re_firstorder(int iband_local, int inf_local);

    void _set_vertex_energydensity(int iband_local, int inf_local) const;

    void _get_face_matrix() const;
    void _get_face_gradient(int iband_local, int inf_local) const;
    void _get_non_fourier(int iband_local, int inf_local) const;
    void  _q_non_fourier_2D(int iband_local,int inf_local) const;

    void _set_vertex_temperature() const;
    void _get_Re_macro();

    void _get_CellMatrix() const;
    void _get_gradient(int Use_limiter, int iband_local, int inf_local) const;
    void _get_bound_ee(int iband_local, int inf_local) const;
    void _set_bound_ee_1() const;
    void _set_macro_boundary_condition() const;

    void _recover_temperature(int iband, int inf_local) const;
    void _get_total_energy(int iband, int inf_local) const;
    void _get_heat_flux(int iband, int inf_local) const;

    void _get_bound_temp() const;
    void _get_bound_flux() const;

    bool _get_magin_check_error(int nt, double error_temp_limit, double error_flux_limit);
    void _print_out() const;


    void solve(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_Iterative(int Use_Backup, int max_iter,  int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_firstorder(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_firstorder_Iterative(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_DOM(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_DOM_Iterative(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_DOM_firstorder(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void solve_DOM_firstorder_Iterative(int Use_Backup, int max_iter, int USE_limiter, double error_temp_limit, double error_flux_limit);

    void copy() const;


    StaticBTESynthetic(const StaticBTESynthetic &solver) = delete;
 ~StaticBTESynthetic();

    void _get_bound_ee_firstorder(int iband_local, int inf_local) const;


    void _marco_solution();

    void _set_face_matrix();

    void _set_cell_matrix();

    void _delete_cell_matrix() const;

    void _delete_face_matrix() const;

    void _set_cell_matrix_larger();

    void _get_CellMatrix_larger() const;

    void _get_gradient_larger(int Use_limiter, int iband_local, int inf_local) const;
};

#endif // STATICBTESOLVER_STATICBTESYNTHETIC_H
