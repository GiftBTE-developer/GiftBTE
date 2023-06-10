//
// Created by yuehu on 2023/5/21.
//

#ifndef STATICBTESOLVER_TRANSIENT_H
#define STATICBTESOLVER_TRANSIENT_H
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

using namespace std;

class Transient
        {
        public:
            double deltaT;
            double totalT;
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


            double *Re;


            Transient(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank,double deltaT,double totalT);

            void _set_initial(int Use_Backup) const;

             void _get_explicit_Re(int itime, int spatial_order, int Use_limiter,int iband_local, int inf_local,double deltaTime);

            void _get_bound_ee(int iband_local, int inf_local) const;
            void _set_bound_ee_1() const;

            void _recover_temperature(int iband, int inf_local) const;
            void _get_total_energy(int iband, int inf_local) const;
            void _get_heat_flux(int iband, int inf_local) const;



            bool _get_magin_check_error(int nt, double error_temp_limit, double error_flux_limit);
            void _print_out() const;


            void copy() const;

            void _delete_cell_matrix() const;

            void _set_cell_matrix_larger();

            void _get_CellMatrix_larger() const;

            void _get_gradient_larger(int Use_limiter, int iband_local, int inf_local) const;


            void solve(int Use_Backup, double error_temp_limit, double error_flux_limit,double deltaT,double totalT);

            void solve_first_order(int Use_Backup, int Use_Limiter, double error_temp_limit, double error_flux_limit,double deltaT,double totalT);



            ~Transient();
        };
#endif //STATICBTESOLVER_TRANSIENT_H
