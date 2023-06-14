//
// Created by yuehu on 2023/2/16.
//

//
// Created by yuehu on 2021/9/8.
//

#include "StaticBTESynthetic/StaticBTESynthetic.h"
#include <algorithm>
#include <iomanip>
#include <chrono>


#ifdef USE_GPU
#include "StaticBTESolver/GPU_kernel.h"
#endif

#ifdef USE_TIME
#include <chrono>
#endif

#ifdef USE_GPU
#define VIENNACL_WITH_CUDA 1
#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#endif
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/amg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/tools/timer.hpp"
#include <fstream>
viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
#endif

using namespace std;

StaticBTESynthetic::~StaticBTESynthetic()
{

    delete [] numCellLocalList;
    delete [] startCellList;



    for (int i = 0; i < numBandLocal; ++i) {
        for (int j = 0; j < numDirectionLocal; ++j) {
          delete []  energyDensity[i][j];
        }
        delete []  energyDensity[i];
    }
    delete []  energyDensity;





    delete [] energyDensityVertex;
    delete [] eboundLocal;
    delete [] ebound;
    delete [] temperature;
    delete [] temperatureOld;
    delete [] temperatureLocal;
    delete [] temperatureVertex;
    delete [] temperature1;

    delete [] totalEnergy;
    delete [] totalEnergyLocal;
    delete [] heatFluxXLocal;
    delete [] heatFluxYLocal;
    delete [] heatFluxZLocal;
    delete [] heatFluxXGlobal;
    delete [] heatFluxYGlobal;
    delete [] heatFluxZGlobal;
    delete [] heatFluxXOld;
    delete [] heatFluxYOld;
    delete [] heatFluxZOld;
    delete [] gradientX;
    delete [] gradientY;
    delete [] gradientZ;
    delete [] gradientFaceX;
    delete [] gradientFaceY;
    delete [] gradientFaceZ;
    delete [] limit;

    for (int i = 0; i < numBound; ++i) {
       delete [] boundaryTempMacro[i];
    }
    delete [] boundaryTempMacro;



    delete [] Re;
    delete [] ReMacro;
    delete [] ReMacroLocal;




    delete [] elementIndex;
    delete [] elementFaceSize;
    delete [] elementFaceNeighobr;
    delete [] elementFaceArea;
    delete [] elementFaceIndex;
    delete [] elementFaceNormX;
    delete [] elementFaceNormY;
    delete [] elementFaceNormZ;
    delete [] elementFaceCenterX;
    delete [] elementFaceCenterY;
    delete [] elementFaceCenterZ;
    delete [] elementFaceBound;
    delete [] elementFaceVertexesSize;
    delete [] elementFaceVertexes;
    delete [] elementVertexesSize;
    delete [] elementVertexes;
    delete [] elementHeatSource;
    delete [] elementVolume;

    delete [] elementCenterX;
    delete [] elementCenterY;
    delete [] elementCenterZ;

    delete [] boundaryType;
    delete [] boundaryThermal;
    delete [] boundaryTemp;
    delete [] boundaryConnect;

    for (int i = 0; i < numBound; ++i) {
       delete [] boundaryCell[i];
       delete [] boundaryFace[i];
    }

    delete [] boundaryCell;
    delete [] boundaryFace;

    for (int i = 0; i < numBound; ++i) {
        for (int j = 0; j < numBand; ++j) {
           delete [] boundaryTrans[i][j];
        }
        delete [] boundaryTrans[i];
    }

    delete [] boundaryTrans;

    delete [] boundaryFlux;

    for (int i = 0; i < numFace; ++i) {
       delete [] faceCell[i];
       delete [] faceFace[i];
    }

    delete [] faceBound;
    delete [] faceCell;
    delete [] faceFace;

    delete [] nodeX;
    delete [] nodeY;
    delete [] nodeZ;

    delete [] matter;
    delete [] kappaBulk;
    delete [] capacityBulk;

    for (int i = 0; i < numofMatter; ++i) {


        for (int j = 0; j < numBand; ++j) {
           delete [] groupVelocityX[i][j];
           delete [] groupVelocityY[i][j];
           delete [] groupVelocityZ[i][j];
           delete [] heatCapacity[i][j];
           delete [] relaxationTime[i][j];
           delete [] latticeRatio[i][j];
           delete [] heatRatio[i][j];
           delete [] modeWeight[i][j];
        }
        delete [] groupVelocityX[i];
        delete [] groupVelocityY[i];
        delete [] groupVelocityZ[i];
        delete [] heatCapacity[i];
        delete [] relaxationTime[i];
        delete [] latticeRatio[i];
        delete [] heatRatio[i];
        delete [] modeWeight[i];
    }

    delete [] groupVelocityX;
    delete [] groupVelocityY;
    delete [] groupVelocityZ;
    delete [] heatCapacity;
    delete [] relaxationTime;
    delete [] latticeRatio;
    delete [] heatRatio;
    delete [] modeWeight;

    delete [] directionX;
    delete [] directionY;
    delete [] directionZ;

    //cout << "~StaticBTESynthetic is activated !!" << endl;
#ifdef USE_GPU

#else

#endif

}
StaticBTESynthetic::StaticBTESynthetic(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank)
{
    this->numProc = num_proc;
    this->worldRank = world_rank;

    //this->mesh=mesh;
    //this->angles=angles;
    //this->bcs=bcs;
    //this->bands=bands;

    L_x=mesh->L_x;
    L_y=mesh->L_y;
    L_z=mesh->L_z;
    dimension=mesh->Dimension;
    numNode=mesh->Nodes.size();
    numBand = bands->nband;
    numDirection = angles->direction.size();
    numCell = mesh->Elements.size();
    numBound=mesh->Boundaries.size();
    numEdge=mesh->nedge;
    numFace=mesh->Element_Faces.size();
    numofMatter=bands->kappabulk.size();

    if (numDirection > num_proc)
    {
        numDirectionLocal = numDirection / num_proc;
        if (world_rank < numDirection % num_proc)
        {
            numDirectionLocal = numDirectionLocal + 1;
        }
        numBandLocal=numBand;
    }
    else
    {
        numDirectionLocal = 1;
        numBandLocal = numBand / (num_proc / numDirection);
    }
    numCellLocal=
            numCell / num_proc;
    if (world_rank < numCell % num_proc)
    {
        numCellLocal = numCellLocal + 1;
        startCell=numCellLocal*world_rank;
    } else
    {
        startCell=(numCell / num_proc+1)*(numCell % num_proc)+(numCell / num_proc)*(world_rank-numCell % num_proc);


    }
    numCellLocalList=new int [num_proc];
    for (int i = 0; i < num_proc; ++i) {
        if(i<numCell % num_proc)
        {
            numCellLocalList[i]=  numCell / num_proc+1;
            //if(world_rank==0)
            //cout<<numCellLocalList[i]<<endl;
        } else
        {
            numCellLocalList[i]=  numCell / num_proc;
            //if(world_rank==0)
            //cout<<numCellLocalList[i]<<endl;
        }
    }
    startCellList=new int [num_proc];
    for (int i = 0; i < num_proc; ++i) {
        if(i<numCell % num_proc)
        {
            startCellList[i]=  (numCell / num_proc+1)*i;
            //if(world_rank==0)
            //cout<<startCellList[i]<<endl;
        } else
        {
            startCellList[i]=(numCell / num_proc+1)*(numCell % num_proc)+(numCell / num_proc)*(i-numCell % num_proc);
            //if(world_rank==0)
            //cout<<startCellList[i]<<endl;
        }
    }

    energyDensity= new double ** [numBandLocal];
    for (int i = 0; i < numBandLocal; ++i) {
        energyDensity[i]=new double * [numDirectionLocal];
        for (int j = 0; j < numDirectionLocal; ++j) {
            energyDensity[i][j]=new double [numCell];
        }
    }
    energyDensityVertex=new double [numNode*numofMatter];
    ebound = new double[numBand * numDirection * numBound * 2];
    eboundLocal = new double[numBandLocal * mesh->Boundaries.size() * 2];
    temperature = new double[numCell];
    temperatureLocal = new double[numCell];
    temperatureOld = new double[numCell];
    temperatureVertex=new double [numNode*numofMatter];
    temperature1=new double [numCell];

    totalEnergyLocal = new double[numCell];
    totalEnergy = new double[numCell];
    heatFluxXOld = new double[numCell];
    heatFluxXLocal = new double[numCell];
    heatFluxXGlobal = new double[numCell];
    heatFluxYOld = new double[numCell];
    heatFluxYLocal = new double[numCell];
    heatFluxYGlobal = new double[numCell];
    heatFluxZOld = new double[numCell];
    heatFluxZLocal = new double[numCell];
    heatFluxZGlobal = new double[numCell];

    gradientX=new double [numCell];
    gradientY=new double [numCell];
    gradientZ=new double [numCell];
    limit=new double [numCell];

    gradientFaceX=new double [numFace*2];
    gradientFaceY=new double [numFace*2];
    gradientFaceZ=new double [numFace*2];



    boundaryTempMacro=new double *[numBound];
    for (int i = 0; i < numBound; ++i) {
        boundaryTempMacro[i]=new double [2];
    }

    Re= new double [numCell];
    ReMacro=new double [numCell];
    ReMacroLocal=new double [numCell];
    Re1.resize(numCell);
    sol1.resize(numCell);




    elementIndex=new int [numCell];
    for (int i = 0; i < numCell; ++i) {
        elementIndex[i]=mesh->Elements[i].index;
    }

    elementFaceSize=new int [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementFaceSize[i]=mesh->Elements[i].faces.size();
    }
    elementFaceNeighobr=new int [numCell * 6];
    elementFaceArea=new double [numCell * 6];
    elementFaceIndex=new int [numCell*6];
    for (int i = 0; i < numCell; ++i) {
        for (int j = 0; j < elementFaceSize[i] ; ++j) {
            elementFaceArea[j+i*6]=mesh->Elements[i].faces[j].area;
            elementFaceNeighobr[j+i*6]=mesh->Elements[i].faces[j].neighbor;
            elementFaceIndex[j+i*6]=mesh->Elements[i].faces[j].index;
        }
    }

    elementFaceNormX=new double [numCell * 6];
    elementFaceNormY=new double [numCell * 6];
    elementFaceNormZ=new double [numCell * 6];

    elementFaceCenterX=new double [numCell * 6];
    elementFaceCenterY=new double [numCell * 6];
    elementFaceCenterZ=new double [numCell * 6];
    elementFaceBound=new int [numCell * 6];
    for (int i = 0; i < numCell * 6 ; ++i) {
        elementFaceNormX[i]=0;
        elementFaceCenterX[i]=0;
        elementFaceNormY[i]=0;
        elementFaceCenterY[i]=0;
        elementFaceNormZ[i]=0;
        elementFaceCenterZ[i]=0;
    }
    for (int i = 0; i < numCell; ++i) {
        for (int j = 0; j < elementFaceSize[i] ; ++j) {
            elementFaceNormX[j+i*6]=mesh->Elements[i].faces[j].norm.x;
            elementFaceNormY[j+i*6]=mesh->Elements[i].faces[j].norm.y;
            elementFaceNormZ[j+i*6]=mesh->Elements[i].faces[j].norm.z;
            elementFaceCenterX[j+i*6]=mesh->Elements[i].faces[j].center.x;
            elementFaceCenterY[j+i*6]=mesh->Elements[i].faces[j].center.y;
            elementFaceCenterZ[j+i*6]=mesh->Elements[i].faces[j].center.z;
            elementFaceBound[j+i*6]=mesh->Elements[i].faces[j].bound;
        }
    }


    elementFaceVertexesSize=new int [numCell * 6];
    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j < elementFaceSize[i]  ; ++j) {
            elementFaceVertexesSize[j+i*6]=mesh->Elements[i].faces[j].vertexes.size();
        }
    }

    elementFaceVertexes=new int [numCell * 6 * 4];
    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j < elementFaceSize[i]  ; ++j) {
            for (int k = 0; k < elementFaceVertexesSize[j+i*6]; ++k) {
                elementFaceVertexes[k+j*4+i*6*4]=mesh->Elements[i].faces[j].vertexes[k];
            }
        }
    }

    elementVertexesSize=new int [numCell];
    for (int i = 0; i < numCell; ++i) {
        elementVertexesSize[i]=mesh->Elements[i].vertexes.size();
    }
    elementVertexes=new int [numCell * 8];
    for (int i = 0; i < numCell; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            elementVertexes[j + i * 8] = mesh->Elements[i].vertexes[j];
            //std::cout<<elementVertexes[j + i * 8]<<std::endl;
        }
    }

    elementVolume=new double [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementVolume[i]=mesh->Elements[i].volume;
    }
    elementHeatSource=new double [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementHeatSource[i]=mesh->Elements[i].heat_source;
    }


    elementCenterX=new double [numCell];
    elementCenterY=new double [numCell];
    elementCenterZ=new double [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementCenterX[i]=mesh->Elements[i].center.x;
        elementCenterY[i]=mesh->Elements[i].center.y;
        elementCenterZ[i]=mesh->Elements[i].center.z;
    }

    elementNeighborList.resize(numCell);
    for (int i = 0; i < numCell; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            int vertex=mesh->Elements[i].vertexes[j];
            for (int k = 0; k < mesh->Nodes[vertex].cells.size(); ++k) {
                int cell=mesh->Nodes[vertex].cells[k];
                if(!ishave(cell,elementNeighborList[i])&&mesh->Elements[cell].matter==mesh->Elements[i].matter&&cell!=i)
                {
                    elementNeighborList[i].push_back(cell);
                }
            }
        }
    }

    boundaryType=new int [numBound];
    boundaryThermal=new double [numBound];
    boundaryTemp=new double [numBound];
    boundaryConnect=new int [numBound];

    for (int i = 0; i <numBound ; ++i) {
        boundaryType[i]=mesh->Boundaries[i].type;
        boundaryThermal[i]=mesh->Boundaries[i].Temperature;

        boundaryConnect[i]=mesh->Boundaries[i].connection;
    }

    boundaryCell=new int *[numBound];
    boundaryFace=new int *[numBound];
    for (int i = 0; i < numBound; ++i) {
        boundaryCell[i]=new int [2];
        boundaryFace[i]=new int [2];
        for (int j = 0; j < 2; ++j) {
            boundaryCell[i][j]=-1;
            boundaryFace[i][j]=-1;
        }
        for (int j = 0; j < mesh->Boundaries[i].cellindex.size(); ++j) {
            boundaryCell[i][j]=mesh->Boundaries[i].cellindex[j];
            boundaryFace[i][j]=mesh->Boundaries[i].faceindex[j];
        }

    }

    boundaryNeighbors.resize(numBound);
    for (int i = 0; i < numBound; ++i) {
        boundaryNeighbors[i]=mesh->Boundaries[i].neighbors;
    }

    boundaryTrans=new double ** [numBound];
    for (int i = 0; i < numBound; ++i) {
        boundaryTrans[i]=new double * [numBand];
        for (int j = 0; j < numBand; ++j) {
            boundaryTrans[i][j]=new double [2];
        }
    }
    for (int i = 0; i < numBound; ++i) {
        for (int j = 0; j < numBand; ++j) {
            if (boundaryType[i]==4||boundaryType[i]==5)
            {
                for (int k = 0; k < 2; ++k) {
                    boundaryTrans[i][j][k]=mesh->Boundaries[i].trans[k][j];
                }

            } else
            {
                for (int k = 0; k < 2; ++k) {
                    boundaryTrans[i][j][k]=0;
                }
            }

        }
    }

    boundaryFlux=new double [numBound];

    faceBound=new int [numFace];
    for (int i = 0; i < numFace; ++i) {
        faceBound[i]=mesh->Element_Faces[i].Boundindex;
    }
    faceCell=new int *[numFace];
    faceFace=new int *[numFace];

    for (int i = 0; i < numFace; ++i) {
        faceCell[i]=new int [2];
        faceFace[i]=new int [2];
        for (int j = 0; j < 2; ++j) {
            faceCell[i][j]=-1;
            faceFace[i][j]=-1;
        }
        for (int j = 0; j < mesh->Element_Faces[i].cellindex.size(); ++j) {
            faceCell[i][j]=mesh->Element_Faces[i].cellindex[j];
            faceFace[i][j]=mesh->Element_Faces[i].faceindex[j];
        }
    }

    nodeX=new double [numNode];
    nodeY=new double [numNode];
    nodeZ=new double [numNode];
    for (int i = 0; i < numNode; ++i) {
        nodeX[i]=mesh->Nodes[i].x;
        nodeY[i]=mesh->Nodes[i].y;
        nodeZ[i]=mesh->Nodes[i].z;
    }


    matter=new int [numCell];
    for (int i = 0; i < numCell; ++i) {
        matter[i]=mesh->Elements[i].matter;
    }
    kappaBulk=new double [numofMatter];
    capacityBulk=new double [numofMatter];
    for (int i = 0; i < numofMatter; ++i) {
        kappaBulk[i]=bands->kappabulk[i];
        capacityBulk[i]=bands->capacitybulk[i];
        //cout<<kappaBulk[i]<<endl;
    }

    groupVelocityX=new double **[numofMatter];
    groupVelocityY=new double **[numofMatter];
    groupVelocityZ=new double **[numofMatter];
    heatCapacity=new double **[numofMatter];
    relaxationTime=new double **[numofMatter];
    latticeRatio=new double ** [numofMatter];
    heatRatio=new double ** [numofMatter];
    modeWeight=new double ** [numofMatter];
    for (int i = 0; i < numofMatter; ++i) {
        groupVelocityX[i]=new double *[numBand];
        groupVelocityY[i]=new double *[numBand];
        groupVelocityZ[i]=new double *[numBand];
        heatCapacity[i]=new double *[numBand];
        relaxationTime[i]=new double *[numBand];
        latticeRatio[i]=new double * [numBand];
        heatRatio[i]=new double * [numBand];
        modeWeight[i]=new double * [numBand];

        for (int j = 0; j < numBand; ++j) {
            groupVelocityX[i][j]=new double [numDirection];
            groupVelocityY[i][j]=new double [numDirection];
            groupVelocityZ[i][j]=new double [numDirection];
            heatCapacity[i][j]=new double [numDirection];
            relaxationTime[i][j]=new double [numDirection];
            latticeRatio[i][j]=new double [numDirection];
            heatRatio[i][j]=new double [numDirection];
            modeWeight[i][j]=new double [numDirection];
        }
    }


    for (int imatter = 0; imatter < numofMatter; ++imatter) {
        for (int i = 0; i < numBand; ++i) {
            for (int j = 0; j < numDirection ; ++j) {
                groupVelocityX[imatter][i][j]=bands->bands[i].group_velocity[imatter]*angles->direction[j].x;
                groupVelocityY[imatter][i][j]=bands->bands[i].group_velocity[imatter]*angles->direction[j].y;
                groupVelocityZ[imatter][i][j]=bands->bands[i].group_velocity[imatter]*angles->direction[j].z;
                heatCapacity[imatter][i][j]=bands->bands[i].heat_capacity[imatter];
                relaxationTime[imatter][i][j]=bands->bands[i].relaxation_time[imatter];
                latticeRatio[imatter][i][j]=bands->bands[i].lattice_ratio[imatter];
                heatRatio[imatter][i][j]=bands->bands[i].heat_ratio[imatter];
                modeWeight[imatter][i][j]=angles->weight[j]/ angles->SolidAngle ;
            }
        }
    }


    directionX=new double [numDirection];
    directionY=new double [numDirection];
    directionZ=new double [numDirection];
    for (int j = 0; j < numDirection ; ++j) {
        directionX[j]=angles->direction[j].x;
        directionY[j]=angles->direction[j].y;
        directionZ[j]=angles->direction[j].z;
    }


}

void StaticBTESynthetic::_set_face_matrix()
{
    FaceMatrix=new double **[numFace];
    for (int ifa = 0; ifa < numFace; ++ifa) {
        int i = faceCell[ifa][0];
        int j = faceFace[ifa][0];
        FaceMatrix[ifa]=new double *[3];
        for (int k = 0; k < 3; ++k) {
            FaceMatrix[ifa][k]=new double [elementFaceVertexesSize[j+i*6]];
        }
    }
}
void StaticBTESynthetic::_delete_face_matrix() const
{
    for (int ifa = 0; ifa < numFace; ++ifa) {
        for (int k = 0; k < 3; ++k) {
            delete []  FaceMatrix[ifa][k];
        }
        delete [] FaceMatrix[ifa];
    }
    delete [] FaceMatrix;
}

void StaticBTESynthetic::_set_cell_matrix()
{

    CellMatrix=new double  **[numCell];
    for (int i = 0; i < numCell; ++i) {
        CellMatrix[i]=new double *[3];
        int nump=0;
        for (int j = 0; j <elementFaceSize[i]; ++j)
        {
            if (elementFaceBound[i*6+j] == -1)
            {
                nump++;
            }
        }
        for (int j = 0; j < 3; ++j) {
            CellMatrix[i][j]=new double [nump];
        }
    }
}
void StaticBTESynthetic::_set_cell_matrix_larger()
{

    CellMatrix=new double  **[numCell];
    for (int i = 0; i < numCell; ++i) {
        CellMatrix[i]=new double *[3];
        for (int j = 0; j < 3; ++j) {
            CellMatrix[i][j]=new double [elementNeighborList[i].size()];
        }
    }
}
void StaticBTESynthetic::_delete_cell_matrix() const
{

    for (int i = 0; i < numCell; ++i) {

        for (int j = 0; j < 3; ++j) {
            delete [] CellMatrix[i][j];
        }
        delete []  CellMatrix[i];
    }
    delete [] CellMatrix;
}


void StaticBTESynthetic::_set_matrix(const std::string& Matrix_type) {

    if (Matrix_type=="LU")
    {
        stiffMatrix.resize(numCell, numCell);
        cg=new Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> * [numBandLocal];
        for (int i = 0; i < numBandLocal; ++i) {
            cg[i]=new Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> [numDirectionLocal];
        }

    }
    else if (Matrix_type=="Iterative")
    {
        stiffMatrix.resize(numCell, numCell);
        /*csrRowPtr=new unsigned int **[numBandLocal];
        csrColInd=new unsigned int **[numBandLocal];
        csrVal=new double **[numBandLocal];

        for (int i = 0; i < numBandLocal; ++i) {
            csrRowPtr[i]=new unsigned int *[numDirectionLocal];
            csrColInd[i]=new unsigned int *[numDirectionLocal];
            csrVal[i]=new double *[numDirectionLocal];
            for (int j = 0; j < numDirectionLocal; ++j) {
                csrRowPtr[i][j]= new unsigned int[numCell + 1];
                csrColInd[i][j] = new unsigned int[7 * numCell + 1];
                csrVal[i][j] = new double[7 * numCell + 1];
            }
        }*/

        //csrRowPtr_macro= new unsigned int[numCell + 1];
        //csrColInd_macro = new unsigned int[7 * numCell + 1];
        //csrVal_macro = new double[7 * numCell + 1];
    }
}
void StaticBTESynthetic::_delete_matrix(const string &Matrix_type) const {

    if (Matrix_type=="LU")
    {
        for (int i = 0; i < numBandLocal; ++i) {
           delete [] cg[i];
        }
        delete [] cg;
    }
    else if (Matrix_type=="Iterative")
    {
        /*for (int i = 0; i < numBandLocal; ++i) {

            for (int j = 0; j < numDirectionLocal; ++j) {
                delete []   csrRowPtr[i][j];
                delete []   csrColInd[i][j];
                delete []  csrVal[i][j];
            }
            delete []   csrRowPtr[i];
            delete []   csrColInd[i];
            delete []  csrVal[i];
        }
        delete []   csrRowPtr;
        delete []   csrColInd;
        delete []  csrVal;

        //delete []   csrRowPtr_macro;
       // delete []   csrColInd_macro;
       // delete []  csrVal_macro;
    }*/
    }


}


void StaticBTESynthetic::_set_initial(int Use_Backup) const
{
    for (int iband_local = 0; iband_local < numBandLocal; iband_local++)
    {
        for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++)
        {
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            for (int kk = 0; kk < numCell; kk++)
            {
                energyDensity[iband_local][inf_local][kk]= temperature[kk] * heatCapacity[matter[kk]][iband_local][inf_local];
            }
        }

    }


}

void StaticBTESynthetic::_get_coefficient(int iband_local, int inf_local)
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    std::vector<Tri> tripletList;

    vector<pair<int, double>> compressed_Ke;

    for (int ie = 0; ie < numCell; ++ie)
    {
        compressed_Ke.resize(0);
        double valueindia = 0;
        //int matter = mesh->Elements[ie].matter;
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface)
        {
            double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                    + groupVelocityY[matter[ie]][iband][inf]* elementFaceNormY[jface + ie * 6]
                    + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
            double temp = relaxationTime[matter[ie]][iband][inf] *
                          elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
            if (dotproduct >= 0)
            {
                valueindia += temp;
            }
            else
            {
                if (elementFaceBound[jface + ie * 6] == -1)
                {
                    int neiindex = elementFaceNeighobr[jface + ie * 6];

                    compressed_Ke.emplace_back(neiindex, temp);
                }
                else
                {
                    //int boundindex = elementFaceBound[jface+ie*6];
                    //if (boundaryType[boundindex] == 0)
                    //{
                    //    int neiindex =  elementFaceNeighobr[jface+ie*6];
                    //    compressed_Ke.emplace_back(neiindex, temp);
                    //}
                }
            }
        }
        valueindia += 1.0;
        compressed_Ke.emplace_back(ie, valueindia);
        sort(compressed_Ke.begin(), compressed_Ke.end());

        for (int j = 0; j < compressed_Ke.size(); ++j) {
            tripletList.push_back(Tri(ie, compressed_Ke[j].first, compressed_Ke[j].second));
              }
          }
    stiffMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    cg[iband_local][inf_local].analyzePattern(stiffMatrix);
    cg[iband_local][inf_local].factorize(stiffMatrix);
    cout<<"finish calculate coefficients: "<<iband_local<<" "<<inf_local<<endl;
};

void StaticBTESynthetic::_get_coefficient_macro()
{
    std::vector<Tri> tripletList;
    for (int ie = 0; ie < numCell; ++ie)
    {
        vector<pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(elementFaceSize[ie] + 1);
        double valueindia = 0;
        double ax,ay,az;
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface)
        {
            double temp = kappaBulk [matter[ie]] * elementFaceArea[jface + ie * 6] ;
            if (elementFaceBound[jface + ie * 6] == -1)
            {
                double ax= elementCenterX[elementFaceNeighobr[jface + ie * 6]] - elementCenterX[ie];
                double ay= elementCenterY[elementFaceNeighobr[jface + ie * 6]] - elementCenterY[ie];
                double az= elementCenterZ[elementFaceNeighobr[jface + ie * 6]] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[jface + ie * 6], temp);

            }
            else if (boundaryType[elementFaceBound[jface + ie * 6]] != 3)
            {
                double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
            }
        }
        compressed_Ke.emplace_back(ie, valueindia);
        sort(compressed_Ke.begin(), compressed_Ke.end());
        for (int j = 0; j < compressed_Ke.size(); ++j) {
            tripletList.push_back(Tri(ie, compressed_Ke[j].first, compressed_Ke[j].second));
        }
    }
    stiffMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    cg_macro.analyzePattern(stiffMatrix);
    cg_macro.factorize(stiffMatrix);
}

void StaticBTESynthetic::_get_CellMatrix() const
{
    if(dimension==2)
    {
        for (int i = 0; i < numCell; ++i)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i] ; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }

            }
            vector<double> d1(nump, 0);
            vector<vector<double>> J1;
            J1.resize(nump);
            for (int j = 0; j < nump; ++j) {
                J1[j].resize(2);
            }
            nump = 0;

            for (int j = 0; j <elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    J1[nump][0]=(1.0/L_x *(elementCenterX[elementFaceNeighobr[j+i*6]] - elementCenterX[i]));
                    J1[nump][1]=(1.0/L_x *(elementCenterY[elementFaceNeighobr[j+i*6]] - elementCenterY[i]));
                     //d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);
                    nump++;
                }
            }
            if(nump>2)
            {
                Get_inverse_2D(CellMatrix[i],J1);
            }

        }
    }
    else if (dimension == 3)
    {
        for (int i = 0; i < numCell; ++i)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i] ; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }

            }
            vector<double> d1(nump, 0);
            vector<vector<double>> J1;
            J1.resize(nump);
            for (int j = 0; j < nump; ++j) {
                J1[j].resize(3);
            }
            nump = 0;

            for (int j = 0; j <elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    J1[nump][0]=(1.0/L_x *(elementCenterX[elementFaceNeighobr[j+i*6]] - elementCenterX[i]));
                    J1[nump][1]=(1.0/L_x *(elementCenterY[elementFaceNeighobr[j+i*6]] - elementCenterY[i]));
                    J1[nump][2]=(1.0/L_x *(elementCenterZ[elementFaceNeighobr[j+i*6]] - elementCenterZ[i]));
                    //d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);
                    nump++;
                }
            }
            if(nump>2)
            {
                Get_inverse_3D(CellMatrix[i],J1);
            }

        }
    }
}

void StaticBTESynthetic::_get_CellMatrix_larger() const
{
    if(dimension==2)
    {
        for (int i = 0; i < numCell; ++i)
        {

            vector<double> d1(elementNeighborList[i].size(), 0);
            vector<vector<double>> J1;
            J1.resize(elementNeighborList[i].size());
            for (int j = 0; j < elementNeighborList[i].size(); ++j) {
                J1[j].resize(2);
            }

            for (int j = 0; j <elementNeighborList[i].size(); ++j)
            {

                J1[j][0]=(1.0/L_x *(elementCenterX[elementNeighborList[i][j]] - elementCenterX[i]));
                J1[j][1]=(1.0/L_x *(elementCenterY[elementNeighborList[i][j]] - elementCenterY[i]));
                //d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);

            }
            if(elementNeighborList[i].size()>2)
            {
                Get_inverse_2D(CellMatrix[i],J1);
            }

        }
    }
    else if (dimension == 3)
    {
        for (int i = 0; i < numCell; ++i)
        {

            vector<double> d1(elementNeighborList[i].size(), 0);
            vector<vector<double>> J1;
            J1.resize(elementNeighborList[i].size());
            for (int j = 0; j < elementNeighborList[i].size(); ++j) {
                J1[j].resize(3);
            }

            for (int j = 0; j <elementNeighborList[i].size(); ++j)
            {

                J1[j][0]=(1.0/L_x *(elementCenterX[elementNeighborList[i][j]] - elementCenterX[i]));
                J1[j][1]=(1.0/L_x *(elementCenterY[elementNeighborList[i][j]] - elementCenterY[i]));
                J1[j][2]=(1.0/L_x *(elementCenterZ[elementNeighborList[i][j]] - elementCenterZ[i]));

                //d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);

            }
            if(elementNeighborList[i].size()>3)
            {
                Get_inverse_3D(CellMatrix[i],J1);
            }

        }
    }
}

void StaticBTESynthetic::_get_gradient(int Use_limiter, int iband_local, int inf_local) const {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int i = 0; i < numCell; ++i)
    {
        gradientX[i] = 0;
        gradientY[i] = 0;
        gradientZ[i] = 0;
        limit[i] = 1;
    }
    if (dimension == 1)
    {
        for (int ie = 0; ie < numCell; ++ie)
        {
            if ( elementFaceBound[ie*6]==-1 && elementFaceBound[ie*6+1]==-1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie] - energyDensity[iband_local][inf_local][ie-1]) / (elementVolume[ie] / 2 + elementVolume[ie-1] / 2);
                double s2 = (energyDensity[iband_local][inf_local][ie+1] - energyDensity[iband_local][inf_local][ie]) / (elementVolume[ie]  / 2 + elementVolume[ie+1] / 2);
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
            else if (elementFaceBound[ie*6] != -1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie+1] - energyDensity[iband_local][inf_local][ie]) / (elementVolume[ie]  / 2 + elementVolume[ie+1] / 2);
                double s2 = s1;
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
            else if (elementFaceBound[ie*6+1] != -1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie] - energyDensity[iband_local][inf_local][ie-1]) / (elementVolume[ie]  / 2 + elementVolume[ie-1] / 2);
                double s2 = s1;
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
        }
    }
    else if (dimension == 2)
    {
        for (int i = 0; i < numCell; ++i)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i] ; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }
                //else if (boundaryType[elementFaceBound[i*6+j]] == 0)
                //{
                //    nump++;
                //}
            }
            if (nump < 3)
            {
                gradientX[i] = 0;
                gradientY[i] = 0;
                gradientZ[i] = 0;
            }
            else
            {
                vector<double> d1(nump, 0);
                //vector<vector<double>> J1(nump);
                nump = 0;
                for (int j = 0; j <elementFaceSize[i]; ++j)
                {
                    if (elementFaceBound[i*6+j] == -1)
                    {
                        d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);
                        nump++;
                    }
                }
                vec gradient;
                gradient.x = 0;
                gradient.y = 0;
                gradient.z = 0;
                for (int m = 0; m < d1.size(); ++m)
                {
                    gradient.x+= CellMatrix[i][0][m] * d1[m];
                    gradient.y+= CellMatrix[i][1][m] * d1[m];
                }
                if (isnan(gradient.x) || isnan(gradient.y) || isnan(gradient.z))
                {
                    gradientX[i] = 0;
                    gradientY[i] = 0;
                    gradientZ[i] = 0;
                } else
                {
                    gradientX[i] = gradient.x;
                    gradientY[i] = gradient.y;
                    gradientZ[i] = gradient.z;
                }

            }
        }
    }
    else if (dimension == 3)
    {
        for (int i = 0; i < numCell; ++i)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i] ; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }
                //else if (boundaryType[elementFaceBound[i*6+j]] == 0)
                //{
                //    nump++;
                //}
            }
            if (nump < 3)
            {
                gradientX[i] = 0;
                gradientY[i] = 0;
                gradientZ[i] = 0;
            }
            else
            {
                vector<double> d1(nump, 0);
                //vector<vector<double>> J1(nump);
                nump = 0;
                for (int j = 0; j <elementFaceSize[i]; ++j)
                {
                    if (elementFaceBound[i*6+j] == -1)
                    {
                        d1[nump] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] - energyDensity[iband_local][inf_local][i]);
                        nump++;
                    }
                }
                vec gradient;
                gradient.x = 0;
                gradient.y = 0;
                gradient.z = 0;
                for (int m = 0; m < d1.size(); ++m)
                {
                    gradient.x+= CellMatrix[i][0][m] * d1[m];
                    gradient.y+= CellMatrix[i][1][m] * d1[m];
                    gradient.z+= CellMatrix[i][2][m] * d1[m];
                }
                if (isnan(gradient.x) || isnan(gradient.y) || isnan(gradient.z))
                {
                    gradientX[i] = 0;
                    gradientY[i] = 0;
                    gradientZ[i] = 0;
                } else
                {
                    gradientX[i] = gradient.x;
                    gradientY[i] = gradient.y;
                    gradientZ[i] = gradient.z;
                }

            }
        }
    }
    if (Use_limiter == 1)
    {
        for (int i = 0; i < numCell; ++i)
        {
            double max = energyDensity[iband_local][inf_local][i];
            double min = energyDensity[iband_local][inf_local][i];
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    if (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] > max)
                    {
                        max = energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]];
                    }
                    if (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] < min)
                    {
                        min = energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]];
                    }
                }
                else if (boundaryType[elementFaceBound[i*6+j]] != 3 && groupVelocityX[matter[i]][iband][inf] * elementFaceNormX[j+i*6]+groupVelocityY[matter[i]][iband][inf]
                * elementFaceNormY[j+i*6]+groupVelocityZ[matter[i]][iband][inf] * elementFaceNormZ[j+i*6] < 0)
                {
                    if (ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0] > max)
                    {
                        max = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0];
                    }
                    if (ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0] < min)
                    {
                        min = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0];
                    }
                }
            }
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                double ax = elementFaceCenterX[i*6+j] - elementCenterX[i];
                double ay = elementFaceCenterY[i*6+j] - elementCenterY[i];
                double az = elementFaceCenterZ[i*6+j] - elementCenterZ[i];
                if ((gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az) * limit[i] + energyDensity[iband_local][inf_local][i] < min)
                {
                    double y = (min - energyDensity[iband_local][inf_local][i]) / (gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az);

                    limit[i] = (pow(y, 2) + 2 * y) / (pow(y, 2) + y + 2);
                }
                if ((gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az) * limit[i] + energyDensity[iband_local][inf_local][i] > max)
                {
                    double y = (max - energyDensity[iband_local][inf_local][i]) / (gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az);
                    limit[i] = (pow(y, 2) + 2 * y) / (pow(y, 2) + y + 2);
                }
            }
        }
    }

}

void StaticBTESynthetic::_get_gradient_larger(int Use_limiter, int iband_local, int inf_local) const {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int i = 0; i < numCell; ++i)
    {
        gradientX[i] = 0;
        gradientY[i] = 0;
        gradientZ[i] = 0;
        limit[i] = 1;
    }
    if (dimension == 1)
    {
        for (int ie = 0; ie < numCell; ++ie)
        {
            if ( elementFaceBound[ie*6]==-1 && elementFaceBound[ie*6+1]==-1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie] - energyDensity[iband_local][inf_local][ie-1]) / (elementVolume[ie] / 2 + elementVolume[ie-1] / 2);
                double s2 = (energyDensity[iband_local][inf_local][ie+1] - energyDensity[iband_local][inf_local][ie]) / (elementVolume[ie]  / 2 + elementVolume[ie+1] / 2);
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
            else if (elementFaceBound[ie*6] != -1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie+1] - energyDensity[iband_local][inf_local][ie]) / (elementVolume[ie]  / 2 + elementVolume[ie+1] / 2);
                double s2 = s1;
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
            else if (elementFaceBound[ie*6+1] != -1)
            {
                double s1 = (energyDensity[iband_local][inf_local][ie] - energyDensity[iband_local][inf_local][ie-1]) / (elementVolume[ie]  / 2 + elementVolume[ie-1] / 2);
                double s2 = s1;
                if ((abs(s1) + abs(s2)) != 0)
                    gradientX[ie] = (sgn(s1) + sgn(s2)) * abs(s1) * abs(s2) / (abs(s1) + abs(s2));
                else
                    gradientX[ie] = 0;
            }
        }
    }
    else if (dimension == 2)
    {
        for (int i = 0; i < numCell; ++i)
        {
            if (elementNeighborList[i].size() < 3)
            {
                gradientX[i] = 0;
                gradientY[i] = 0;
                gradientZ[i] = 0;
            }
            else
            {
                vector<double> d1(elementNeighborList[i].size() , 0);
                //vector<vector<double>> J1(nump);
                for (int j = 0; j <elementNeighborList[i].size() ; ++j)
                {
                    d1[j] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementNeighborList[i][j]] - energyDensity[iband_local][inf_local][i]);
                }
                vec gradient;
                gradient.x = 0;
                gradient.y = 0;
                gradient.z = 0;
                for (int m = 0; m < d1.size(); ++m)
                {
                    gradient.x+= CellMatrix[i][0][m] * d1[m];
                    gradient.y+= CellMatrix[i][1][m] * d1[m];
                    //cout<<CellMatrix[i][0][m]<<" "<<CellMatrix[i][1][m]<<endl;
                }
                if (isnan(gradient.x) || isnan(gradient.y) || isnan(gradient.z))
                {
                    gradientX[i] = 0;
                    gradientY[i] = 0;
                    gradientZ[i] = 0;
                } else
                {
                    gradientX[i] = gradient.x;
                    gradientY[i] = gradient.y;
                    gradientZ[i] = gradient.z;
                }

            }
        }
    }
    else if (dimension == 3)
    {
        for (int i = 0; i < numCell; ++i)
        {
            if (elementNeighborList[i].size() < 3)
            {
                gradientX[i] = 0;
                gradientY[i] = 0;
                gradientZ[i] = 0;
            }
            else
            {
                vector<double> d1(elementNeighborList[i].size() , 0);
                //vector<vector<double>> J1(nump);
                for (int j = 0; j <elementNeighborList[i].size() ; ++j)
                {
                    d1[j] = 1.0/L_x * (energyDensity[iband_local][inf_local][elementNeighborList[i][j]] - energyDensity[iband_local][inf_local][i]);
                }
                vec gradient;
                gradient.x = 0;
                gradient.y = 0;
                gradient.z = 0;
                for (int m = 0; m < d1.size(); ++m)
                {
                    gradient.x+= CellMatrix[i][0][m] * d1[m];
                    gradient.y+= CellMatrix[i][1][m] * d1[m];
                    gradient.z+= CellMatrix[i][2][m] * d1[m];
                }
                if (isnan(gradient.x) || isnan(gradient.y) || isnan(gradient.z))
                {
                    gradientX[i] = 0;
                    gradientY[i] = 0;
                    gradientZ[i] = 0;
                } else
                {
                    gradientX[i] = gradient.x;
                    gradientY[i] = gradient.y;
                    gradientZ[i] = gradient.z;
                }

            }
        }
    }
    if (Use_limiter == 1)
    {
        for (int i = 0; i < numCell; ++i)
        {
            double max = energyDensity[iband_local][inf_local][i];
            double min = energyDensity[iband_local][inf_local][i];
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    if (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] > max)
                    {
                        max = energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]];
                    }
                    if (energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]] < min)
                    {
                        min = energyDensity[iband_local][inf_local][elementFaceNeighobr[j+i*6]];
                    }
                }
                else if (boundaryType[elementFaceBound[i*6+j]] != 3 && groupVelocityX[matter[i]][iband][inf] * elementFaceNormX[j+i*6]+groupVelocityY[matter[i]][iband][inf]
                * elementFaceNormY[j+i*6]+groupVelocityZ[matter[i]][iband][inf] * elementFaceNormZ[j+i*6] < 0)
                {
                    if (ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0] > max)
                    {
                        max = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0];
                    }
                    if (ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0] < min)
                    {
                        min = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + elementFaceBound[i * 6 + j] * 2 + 0];
                    }
                }
            }
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                double ax = elementFaceCenterX[i*6+j] - elementCenterX[i];
                double ay = elementFaceCenterY[i*6+j] - elementCenterY[i];
                double az = elementFaceCenterZ[i*6+j] - elementCenterZ[i];
                if ((gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az) * limit[i] + energyDensity[iband_local][inf_local][i] < min)
                {
                    double y = (min - energyDensity[iband_local][inf_local][i]) / (gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az);

                    limit[i] = (pow(y, 2) + 2 * y) / (pow(y, 2) + y + 2);
                }
                if ((gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az) * limit[i] + energyDensity[iband_local][inf_local][i] > max)
                {
                    double y = (max - energyDensity[iband_local][inf_local][i]) / (gradientX[i]*ax+ gradientY[i]*ay + gradientZ[i]*az);
                    limit[i] = (pow(y, 2) + 2 * y) / (pow(y, 2) + y + 2);
                }
            }
        }
    }

}


void StaticBTESynthetic::_get_Re(int iband_local, int inf_local)
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    for (int i = 0; i < numCell; ++i) {
        Re[i]=0;
    }

    //internal
    for (int ie = 0; ie < numCell; ++ie) {
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface) {
            double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                    + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6]
                    + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
            double temp = relaxationTime[matter[ie]][iband][inf] *
                          elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
            if (dotproduct >= 0) {
                double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                double e = (energyDensity[iband_local][inf_local][ie] + (ax * gradientX[ie] + ay * gradientY[ie] + az * gradientZ[ie]) * limit[ie]);
                Re[ie] -= temp * e;
            } else {
                if (elementFaceBound[jface + ie * 6] == -1) {
                    int neiindex = elementFaceNeighobr[jface + ie * 6];
                    double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[neiindex];
                    double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[neiindex];
                    double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[neiindex];
                    double e = (energyDensity[iband_local][inf_local][neiindex] + (ax * gradientX[neiindex]+ay * gradientY[neiindex]+az * gradientZ[neiindex]) * limit[neiindex]);
                    Re[ie] -= temp * e;
                }
            }
        }
        Re[ie] += temperatureOld[ie] * heatCapacity[matter[ie]][iband][inf];
        Re[ie] += elementHeatSource[ie] * heatRatio[matter[ie]][iband][inf] *
                  relaxationTime[matter[ie]][iband][inf];
        Re[ie] -= energyDensity[iband_local][inf_local][ie];


    }

    for (int ib = 0; ib < numBound; ++ib) {
        for (int icell = 0; icell < 2; ++icell) {
            int ie = boundaryCell[ib][icell];
            int jface = boundaryFace[ib][icell];
            if(ie >= 0)
            {
                double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                        + groupVelocityY[matter[ie]][iband][inf]* elementFaceNormY[jface + ie * 6]
                        + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);

                if (dotproduct < 0) {
                    double temp = relaxationTime[matter[ie]][iband][inf] *
                                  elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
                    Re[ie] -= temp * ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell];
                    if(isnan(Re[ie]))
                    {
                        cout<<ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell]<<endl;

                        cout<<Re[ie]<<endl;
                         }

                }
            }
        }

    }

    for (int i = 0; i < numCell; ++i) {
        Re1(i)=Re[i];
    }

}

void StaticBTESynthetic::_get_bound_ee(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {
            int ie = boundaryCell[ib][icell];
            int jface = boundaryFace[ib][icell];
            if(ie>=0)
            {
                double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf]
                        * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                if (dotproduct>= 0)
                {
                    double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                    double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                    double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                    double e = (energyDensity[iband_local][inf_local][ie] + (ax * gradientX[ie] + ay * gradientY[ie] + az * gradientZ[ie]) * limit[ie]);
                    eboundLocal[iband_local * numBound * 2 + ib * 2 + icell] = e;
                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell] = e;
                    if(isnan(e))
                    {
                        cout<<e<<endl;
                    }
                }
            }
        }

    }
    /*MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);*/
}

void StaticBTESynthetic::_set_vertex_energydensity(int iband_local, int inf_local) const
{
    double *temp=new double [numNode*numofMatter];
    double *R_n=new double [numNode*numofMatter];

    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    for (int i = 0; i < numNode*numofMatter ; ++i) {
        temp[i]=0;
        R_n[i]=0;
    }

    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            double distance= sqrt(pow(nodeX[elementVertexes[j+i*8]]-elementCenterX[i],2)+
                    pow(nodeY[elementVertexes[j+i*8]]-elementCenterY[i],2)+
                    pow(nodeZ[elementVertexes[j+i*8]]-elementCenterZ[i],2));
            temp[elementVertexes[j+i*8]*numofMatter+matter[i]]+=energyDensity[iband_local][inf_local][i]/
                    distance;
            R_n[elementVertexes[j+i*8]*numofMatter+matter[i]]+=1.0/distance;
        }
    }

    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {
            int i = boundaryCell[ib][icell];
            int j = boundaryFace[ib][icell] ;
            if(i>=0)
            {
                for (int k = 0; k < elementFaceVertexesSize[j+i*6]; ++k)
                {

                    int index = elementFaceVertexes[k+j*4+i*6*4];
                    temp[index*numofMatter+matter[i]] = 0;
                    R_n[index*numofMatter+matter[i]] = 0;
                }
            }
        }


    }
    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {
            int i = boundaryCell[ib][icell];
            int j = boundaryFace[ib][icell] ;
            if(i>=0)
            {
                for (int k = 0; k < elementFaceVertexesSize[j+i*6]; ++k)
                {
                    int index = elementFaceVertexes[k+j*4+i*6*4];
                    double distance= sqrt(pow(nodeX[index]-elementFaceCenterX[j+i*6],2)+
                            pow(nodeY[index]-elementFaceCenterY[j+i*6],2)+
                            pow(nodeZ[index]-elementFaceCenterZ[j+i*6],2));
                    temp[index*numofMatter+matter[i]] += ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell] / distance;
                    R_n[index*numofMatter+matter[i]] += 1.0 / distance;
                }
            }
        }
    }
    for (int i = 0; i < numNode*numofMatter ; ++i) {
        energyDensityVertex[i]=temp[i]/R_n[i];
        if(isnan(energyDensityVertex[i]))
        {
            energyDensityVertex[i]=0;
        }
        //cout<<temp[i]<<" "<<R_n[i]<<endl;
    }

    delete[] temp;
    delete[] R_n;


}
void StaticBTESynthetic::_get_face_matrix() const
{
    if(dimension==3)
    {
        for (int ifa = 0; ifa < numFace; ++ifa)
        {

                int i = faceCell[ifa][0];
                int j = faceFace[ifa][0];

                std::vector<double> d1(elementFaceVertexesSize[j+i*6], 0);
                std::vector<std::vector<double>> J1(elementFaceVertexesSize[j+i*6]);
                for (int k = 0; k < elementFaceVertexesSize[j+i*6]; ++k) {
                    J1[k].resize(3);
                }
                J1[0][0]=elementFaceNormX[j+i*6];
                J1[0][1]=(elementFaceNormY[j+i*6]);
                J1[0][2]=(elementFaceNormZ[j+i*6]);
                //d1[0] = 0;

                for (int m = 0; m <  elementFaceVertexesSize[j+i*6]-1; m++)
                {
                    J1[m + 1][0]=(1.0/L_x * (nodeX[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeX[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1][1]=(1.0/L_x * (nodeY[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeY[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1][2]=(1.0/L_x * (nodeZ[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeZ[elementFaceVertexes[m+j*4+i*6*4]]));
                    //d1[m + 1] = 1.0/L_x * (energyDensityVertex[elementFaceVertexes[m+1+j*4+i*6*4]] - energyDensityVertex[elementFaceVertexes[m+j*4+i*6*4]] );
                }
                Get_inverse_3D(FaceMatrix[ifa],J1);
            }
            //vec a = Get_gradient_3D(J1, d1);
        }

}

void StaticBTESynthetic::_get_face_gradient(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ifa = 0; ifa < numFace; ++ifa)
    {
        double gradient_X=0;
        double gradient_Y=0;
        double gradient_Z=0;
        if (faceBound[ifa]==-1)
        {
            int i = faceCell[ifa][0];
            int j = faceFace[ifa][0];
            int i2 = faceCell[ifa][1];
            int j2 = faceFace[ifa][1];
            double ax=elementCenterX[i2]-elementCenterX[i];
            double ay=elementCenterY[i2]-elementCenterY[i];
            double az=elementCenterZ[i2]-elementCenterZ[i];
            double gradientn;
            gradientn=(energyDensity[iband_local][inf_local][i2]-energyDensity[iband_local][inf_local][i])*
                    (ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6])/
                    (ax*ax+ay*ay+az*az);
            if (dimension == 2)
            {
                double mX=nodeX[elementFaceVertexes[1+j*4+i*6*4]]-nodeX[elementFaceVertexes[0+j*4+i*6*4]];
                double mY=nodeY[elementFaceVertexes[1+j*4+i*6*4]]-nodeY[elementFaceVertexes[0+j*4+i*6*4]];
                gradient_X=(energyDensityVertex[elementFaceVertexes[1+j*4+i*6*4]*numofMatter+matter[i]]-energyDensityVertex[elementFaceVertexes[0+j*4+i*6*4]*numofMatter+matter[i]])/(mX*mX+mY*mY)*mX;
                gradient_Y=(energyDensityVertex[elementFaceVertexes[1+j*4+i*6*4]*numofMatter+matter[i]]-energyDensityVertex[elementFaceVertexes[0+j*4+i*6*4]*numofMatter+matter[i]])/(mX*mX+mY*mY)*mY;
            }
            else if (dimension == 3)
            {
                std::vector<double> d1(elementFaceVertexesSize[j+i*6], 0);
               d1[0] = 0;

                for (int m = 0; m <  elementFaceVertexesSize[j+i*6]-1; m++)
                {
                    d1[m + 1] = 1.0/L_x * (energyDensityVertex[elementFaceVertexes[m+1+j*4+i*6*4]*numofMatter+matter[i]] - energyDensityVertex[elementFaceVertexes[m+j*4+i*6*4]*numofMatter+matter[i]] );
                }

                gradient_X = 0;
                gradient_Y = 0;
                gradient_Z = 0;
                for (int m = 0; m < d1.size(); ++m)
                {
                    gradient_X+= FaceMatrix[ifa][0][m] * d1[m];
                    gradient_Y+= FaceMatrix[ifa][1][m] * d1[m];
                    gradient_Z+= FaceMatrix[ifa][2][m] * d1[m];
                }
            }
            if (isnan(gradient_X)||isnan(gradient_Y)||isnan(gradient_Z))
            {
                cout<<gradient_X<<" "<<gradient_Y<<" "<<gradient_Z<<endl;
            }
            gradient_X = gradientn * elementFaceNormX[j+i*6] + gradient_X;
            gradient_Y = gradientn * elementFaceNormY[j+i*6] + gradient_Y;
            gradient_Z = gradientn * elementFaceNormZ[j+i*6] + gradient_Z;
            gradientFaceX[ifa*2]=gradient_X;
            gradientFaceY[ifa*2]=gradient_Y;
            gradientFaceZ[ifa*2]=gradient_Z;
            gradientFaceX[ifa*2+1]=gradient_X;
            gradientFaceY[ifa*2+1]=gradient_Y;
            gradientFaceZ[ifa*2+1]=gradient_Z;
        }
        else if (boundaryType[faceBound[ifa]]!=3)
        {
            for (int icell = 0; icell < 2; ++icell) {
                int ie = faceCell[ifa][icell];
                int jface = faceFace[ifa][icell];
                int ibcell=0;
                for (int ik = 0; ik < 2; ++ik) {
                    if(boundaryCell[faceBound[ifa]][ik]==ie)
                    {
                        ibcell=ik;
                    }
                }
                if(ie >= 0)
                {
                    double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                    double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                    double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                    double gradientn;
                    double eface = 0;
                    eface = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + faceBound[ifa] * 2 + ibcell];
                    gradientn= (eface-energyDensity[iband_local][inf_local][ie]) *
                            (ax*elementFaceNormX[jface + ie * 6] + ay * elementFaceNormY[jface + ie * 6] + az * elementFaceNormZ[jface + ie * 6]) /
                            (ax*ax+ay*ay+az*az);
                    if (dimension == 2)
                    {
                        double mX= nodeX[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeX[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                        double mY= nodeY[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeY[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                        gradient_X= (energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - energyDensityVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]]) / (mX * mX + mY * mY) * mX;
                        gradient_Y= (energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - energyDensityVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]]) / (mX * mX + mY * mY) * mY;

                        if (isnan(gradient_X)||isnan(gradient_Y)||isnan(gradient_Z))
                        {
                            cout << elementFaceVertexes[1 + jface * 4 + ie * 6 * 4] << " " << elementFaceVertexes[jface * 4 + ie * 6 * 4] << endl;
                            cout << energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] << " " << energyDensityVertex[elementFaceVertexes[jface * 4 + ie * 6 * 4]] << endl;
                            cout<<gradient_X<<" "<<gradient_Y<<" "<<gradient_Z<<endl;
                        }
                    }
                    else if (dimension == 3)
                    {
                        std::vector<double> d1(elementFaceVertexesSize[jface + ie * 6], 0);
                        d1[0] = 0;

                        for (int m = 0; m < elementFaceVertexesSize[jface + ie * 6] - 1; m++)
                        {
                            d1[m + 1] = 1.0/L_x * (energyDensityVertex[elementFaceVertexes[m + 1 + jface * 4 + ie * 6 * 4]] - energyDensityVertex[elementFaceVertexes[m + jface * 4 + ie * 6 * 4]] );
                        }

                        gradient_X = 0;
                        gradient_Y = 0;
                        gradient_Z = 0;
                        for (int m = 0; m < d1.size(); ++m)
                        {
                            gradient_X+= FaceMatrix[ifa][0][m] * d1[m];
                            gradient_Y+= FaceMatrix[ifa][1][m] * d1[m];
                            gradient_Z+= FaceMatrix[ifa][2][m] * d1[m];
                        }
                    }


                    gradient_X = gradientn * elementFaceNormX[jface + ie * 6] + gradient_X;
                    gradient_Y = gradientn * elementFaceNormY[jface + ie * 6] + gradient_Y;
                    gradient_Z = gradientn * elementFaceNormZ[jface + ie * 6] + gradient_Z;

                    gradientFaceX[ifa*2+icell]=gradient_X;
                    gradientFaceY[ifa*2+icell]=gradient_Y;
                    gradientFaceZ[ifa*2+icell]=gradient_Z;
                }

            }
        }
        /*else if (boundaryType[faceBound[ifa]]!=3)
        {
            for (int icell = 0; icell < 2; ++icell) {
                int ie = faceCell[ifa][icell];
                int jface = faceFace[ifa][icell];
                int ibcell=0;
                for (int ik = 0; ik < 2; ++ik) {
                    if(boundaryCell[faceBound[ifa]][ik]==ie)
                    {
                        ibcell=ik;
                    }
                }
                if(ie >= 0)
                {
                    double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                    double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                    double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                    double gradientn;
                    double eface = 0;
                    eface = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + faceBound[ifa] * 2 + ibcell];
                    gradientn= (eface-energyDensity[iband_local][inf_local][ie]) *
                            (ax*elementFaceNormX[jface + ie * 6] + ay * elementFaceNormY[jface + ie * 6] + az * elementFaceNormZ[jface + ie * 6]) /
                            (ax*ax+ay*ay+az*az);
                    if (dimension == 2)
                    {
                        double mX= nodeX[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeX[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                        double mY= nodeY[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeY[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                        gradient_X= (energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] - energyDensityVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]]) / (mX * mX + mY * mY) * mX;
                        gradient_Y= (energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] - energyDensityVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]]) / (mX * mX + mY * mY) * mY;

                        if (isnan(gradient_X)||isnan(gradient_Y)||isnan(gradient_Z))
                        {
                            cout << elementFaceVertexes[1 + jface * 4 + ie * 6 * 4] << " " << elementFaceVertexes[jface * 4 + ie * 6 * 4] << endl;
                            cout << energyDensityVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] << " " << energyDensityVertex[elementFaceVertexes[jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] << endl;
                            cout<<gradient_X<<" "<<gradient_Y<<" "<<gradient_Z<<endl;
                        }
                    }
                    else if (dimension == 3)
                    {
                        std::vector<double> d1(boundaryNeighbors[faceBound[ifa]].size()+1, 0);
                        std::vector<std::vector<double>> J1(boundaryNeighbors[faceBound[ifa]].size()+1);
                        for (int k = 0; k < boundaryNeighbors[faceBound[ifa]].size()+1; ++k) {
                            J1[k].resize(3);
                        }
                        J1[0][0]=elementFaceNormX[jface+ie*6];
                        J1[0][1]=(elementFaceNormY[jface+ie*6]);
                        J1[0][2]=(elementFaceNormZ[jface+ie*6]);
                        d1[0] = 0;

                        for (int m = 0; m < boundaryNeighbors[faceBound[ifa]].size(); ++m) {
                            for (int icell1 = 0; icell1 < 2; ++icell1) {
                                int ie1=boundaryCell[boundaryNeighbors[faceBound[ifa]][m]][icell1];
                                int jface1=boundaryFace[boundaryNeighbors[faceBound[ifa]][m]][icell1];
                                if(ie1>=0)
                                {
                                    if(elementIndex[ie1]==elementIndex[ie])
                                    {
                                        J1[m + 1][0]=1.0/L_x *(elementFaceCenterX[jface1+ie1*6]-elementFaceCenterX[jface+ie*6]);
                                        J1[m + 1][1]=1.0/L_x *(elementFaceCenterY[jface1+ie1*6]-elementFaceCenterY[jface+ie*6]);
                                        J1[m + 1][2]=1.0/L_x *(elementFaceCenterZ[jface1+ie1*6]-elementFaceCenterZ[jface+ie*6]);
                                        d1[m + 1] = 1.0/L_x * (ebound[iband*numDirection*numBound*2+inf*numBound*2+boundaryNeighbors[faceBound[ifa]][m]*2+icell1]-eface);

                                    }
                                }

                            }

                        }
                        vec a= Get_gradient_3D(J1,d1);
                        gradient_X = 0;
                        gradient_Y = 0;
                        gradient_Z = 0;
                        for (int m = 0; m < d1.size(); ++m)
                        {
                            gradient_X+= FaceMatrix[ifa][0][m] * d1[m];
                            gradient_Y+= FaceMatrix[ifa][1][m] * d1[m];
                            gradient_Z+= FaceMatrix[ifa][2][m] * d1[m];
                        }
                        gradient_X = a.x;
                        gradient_Y = a.y;
                        gradient_Z = a.z;

                    }


                    gradient_X = gradientn * elementFaceNormX[jface + ie * 6] + gradient_X;
                    gradient_Y = gradientn * elementFaceNormY[jface + ie * 6] + gradient_Y;
                    gradient_Z = gradientn * elementFaceNormZ[jface + ie * 6] + gradient_Z;

                    gradientFaceX[ifa*2+icell]=gradient_X;
                    gradientFaceY[ifa*2+icell]=gradient_Y;
                    gradientFaceZ[ifa*2+icell]=gradient_Z;
                }

            }
        }*/
    }
}

void  StaticBTESynthetic::_get_non_fourier(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    for (int ifa = 0; ifa < numFace; ++ifa)
    {
        if (faceBound[ifa]==-1)
        {
            int ie = faceCell[ifa][0];
            int jface = faceFace[ifa][0];
            int ie2 = faceCell[ifa][1];
            int jface2 = faceFace[ifa][1];
            if(heatRatio[matter[ie]][iband][inf]>0)
            {
                ReMacroLocal[ie]-= modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * relaxationTime[matter[ie]][iband][inf] *
                        (groupVelocityX[matter[ie]][iband][inf] * gradientFaceX[ifa*2] + groupVelocityY[matter[ie]][iband][inf] * gradientFaceY[ifa*2] + groupVelocityZ[matter[ie]][iband][inf] * gradientFaceZ[ifa*2])
                        * (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                ReMacroLocal[ie]+= kappaBulk[matter[ie]] * modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * latticeRatio[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf]
                        * (elementFaceNormX[jface + ie * 6] * gradientFaceX[ifa*2] + elementFaceNormY[jface + ie * 6] * gradientFaceY[ifa*2] + elementFaceNormZ[jface + ie * 6] * gradientFaceZ[ifa*2]);

            }
            if(heatRatio[matter[ie2]][iband][inf]>0)
            {
                ReMacroLocal[ie2]-= modeWeight[matter[ie2]][iband][inf] * elementFaceArea[jface + ie * 6] * relaxationTime[matter[ie2]][iband][inf] *
                        (groupVelocityX[matter[ie2]][iband][inf] * gradientFaceX[ifa*2+1] + groupVelocityY[matter[ie2]][iband][inf] * gradientFaceY[ifa*2+1] + groupVelocityZ[matter[ie2]][iband][inf] * gradientFaceZ[ifa*2+1])
                        * (groupVelocityX[matter[ie2]][iband][inf] * elementFaceNormX[jface2 + ie2 * 6] + groupVelocityY[matter[ie2]][iband][inf] * elementFaceNormY[jface2 + ie2 * 6] + groupVelocityZ[matter[ie2]][iband][inf] * elementFaceNormZ[jface2 + ie2 * 6]);
                ReMacroLocal[ie2]+= kappaBulk[matter[ie2]] * modeWeight[matter[ie2]][iband][inf] * elementFaceArea[jface2 + ie2 * 6] * latticeRatio[matter[ie2]][iband][inf] / heatCapacity[matter[ie2]][iband][inf]
                        * (elementFaceNormX[jface2 + ie2 * 6] * gradientFaceX[ifa*2+1] + elementFaceNormY[jface2 + ie2 * 6] * gradientFaceY[ifa*2+1] + elementFaceNormZ[jface2 + ie2 * 6] * gradientFaceZ[ifa*2+1]);
            }
        }
        else if (boundaryType[faceBound[ifa]]!=3)
        {
            for (int icell = 0; icell < 2; ++icell) {
                int ie = faceCell[ifa][icell];
                int jface = faceFace[ifa][icell];
                if(ie>=0)
                {
                    if(heatRatio[matter[ie]][iband][inf]>0)
                    {
                        ReMacroLocal[ie]-= modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * relaxationTime[matter[ie]][iband][inf] *
                                (groupVelocityX[matter[ie]][iband][inf] * gradientFaceX[ifa*2+icell] + groupVelocityY[matter[ie]][iband][inf] * gradientFaceY[ifa*2+icell] + groupVelocityZ[matter[ie]][iband][inf] * gradientFaceZ[ifa*2+icell])
                                * (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                        ReMacroLocal[ie]+= kappaBulk[matter[ie]] * modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * latticeRatio[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf]
                                * (elementFaceNormX[jface + ie * 6] * gradientFaceX[ifa*2+icell] + elementFaceNormY[jface + ie * 6] * gradientFaceY[ifa*2+icell] + elementFaceNormZ[jface + ie * 6] * gradientFaceZ[ifa*2+icell]);

                    }

                }
            }

        }

    }
}

/*void  StaticBTESynthetic::_q_non_fourier_2D(int iband_local,int inf_local) const
{
    int inf=((inf_local)*numProc+worldRank)%numDirection;
    int iband=iband_local*(ceil(double(numProc)/double(numDirection)))+worldRank/numDirection;
    vector<double> ee_vertex(bands->kappabulk.size()*mesh->Nodes.size(),0);
    vector<double> tempall(bands->kappabulk.size()*mesh->Nodes.size(),0);
    vector<double> R_n(bands->kappabulk.size()*mesh->Nodes.size(),0);
    for (int i = 0; i <mesh->Nodes.size() ; ++i) {
        for (int j = 0; j <mesh->Nodes[i].cells.size() ; ++j) {
            tempall[i*bands->kappabulk.size()+mesh->Elements[mesh->Nodes[i].cells[j]].matter]+=energyDensity[iband_local][inf_local][mesh->Nodes[i].cells[j]]/ get_distance(mesh->Nodes[i],mesh->Elements[mesh->Nodes[i].cells[j]].center);
            R_n[i*bands->kappabulk.size()+mesh->Elements[mesh->Nodes[i].cells[j]].matter]+=1.0/get_distance(mesh->Nodes[i],mesh->Elements[mesh->Nodes[i].cells[j]].center);
        }

    }
    for (int i = 0; i < mesh->Boundaries.size() ; ++i) {
        for (int j = 0; j < mesh->Boundaries[i].cellindex.size(); ++j) {
            int ie=mesh->Boundaries[i].cellindex[j];
            int iface=mesh->Boundaries[i].faceindex[j];
            for (int k = 0; k <mesh->Elements[ie].faces[iface].vertexes.size() ; ++k) {
                int index=mesh->Elements[ie].faces[iface].vertexes[k];
                tempall[index*bands->kappabulk.size()+mesh->Elements[ie].matter]=0;
                R_n[index*bands->kappabulk.size()+mesh->Elements[ie].matter]=0;
            }
        }
    }
    for (int i = 0; i < mesh->Boundaries.size() ; ++i) {
        for (int j = 0; j < mesh->Boundaries[i].cellindex.size(); ++j) {
            int ie=mesh->Boundaries[i].cellindex[j];
            int iface=mesh->Boundaries[i].faceindex[j];
            for (int k = 0; k <mesh->Elements[ie].faces[iface].vertexes.size() ; ++k) {
                int index=mesh->Elements[ie].faces[iface].vertexes[k];
                tempall[index*bands->kappabulk.size()+mesh->Elements[ie].matter]+=ebound[iband*numDirection*mesh->Boundaries.size()*2+inf*mesh->Boundaries.size()*2+i*2+j]/ get_distance(mesh->Nodes[index],mesh->Elements[ie].faces[iface].center);
                R_n[index*bands->kappabulk.size()+mesh->Elements[ie].matter]+=1.0/ get_distance(mesh->Nodes[index],mesh->Elements[ie].faces[iface].center);
            }
        }
    }

    for (int i = 0; i <mesh->Nodes.size()*bands->kappabulk.size() ; ++i) {
        if  (R_n[i]!=0)
            ee_vertex[i]=tempall[i]/R_n[i];
    }


    vec gradient_face;



    for (int ie = 0; ie < numCell ; ++ie) {
        int imatter=mesh->Elements[ie].matter;
        for (int jface = 0; jface < mesh->Elements[ie].faces.size() ; ++jface) {
            if (mesh->Elements[ie].faces[jface].bound == -1)
            {
                vec gradient_face;
                gradient_face.x=0;
                gradient_face.y=0;
                gradient_face.z=0;
                vec a;
                a= mesh->Elements[mesh->Elements[ie].faces[jface].neighbor].center - mesh->Elements[ie].center;
                double gradientn;
                gradientn= (energyDensity[iband_local][inf_local][mesh->Elements[ie].faces[jface].neighbor] - energyDensity[iband_local][inf_local][ie]) * dot_product(mesh->Elements[ie].faces[jface].norm, a) / dot_product(a, a);

                vector<double> d1(mesh->Elements[ie].faces[jface].vertexes.size(), 0);
                vector<vector<double>> J1(mesh->Elements[ie].faces[jface].vertexes.size());
                J1[0].push_back(mesh->Elements[ie].faces[jface].norm.x);
                J1[0].push_back(mesh->Elements[ie].faces[jface].norm.y);
                d1[0]=0;
                for (int m = 0; m < mesh->Elements[ie].faces[jface].vertexes.size() - 1; m++) {
                    J1[m+1].push_back(mesh->Nodes[mesh->Elements[ie].faces[jface].vertexes[m + 1]].x - mesh->Nodes[mesh->Elements[ie].faces[jface].vertexes[m]].x);
                    J1[m+1].push_back(mesh->Nodes[mesh->Elements[ie].faces[jface].vertexes[m + 1]].y - mesh->Nodes[mesh->Elements[ie].faces[jface].vertexes[m]].y);
                    d1[m+1]= ee_vertex[mesh->Elements[ie].faces[jface].vertexes[m + 1] * bands->kappabulk.size() + mesh->Elements[ie].matter] - ee_vertex[mesh->Elements[ie].faces[jface].vertexes[m] * bands->kappabulk.size() + mesh->Elements[ie].matter];
                    //d1[m+1]=ee_vertex[mesh->Elements[i].faces[j].vertexes[m+1]]-ee_vertex[mesh->Elements[i].faces[j].vertexes[m]];

                }
                gradient_face=Get_gradient_2D(J1,d1);
                gradient_face.x=gradientn*mesh->Elements[ie].faces[jface].norm.x + gradient_face.x;
                gradient_face.y=gradientn*mesh->Elements[ie].faces[jface].norm.y + gradient_face.y;
                if(heatRatio[imatter][iband][inf]>0)
                {
                    ReMacroLocal[ie]-= modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * relaxationTime[matter[ie]][iband][inf] *
                            (groupVelocityX[matter[ie]][iband][inf] * gradient_face.x + groupVelocityY[matter[ie]][iband][inf] * gradient_face.y )
                            * (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                    ReMacroLocal[ie]+= kappaBulk[matter[ie]] * modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * latticeRatio[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf]
                            * (elementFaceNormX[jface + ie * 6] * gradient_face.x + elementFaceNormY[jface + ie * 6] * gradient_face.y );
               }
            }
            else if (mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].type != 3 && mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].type != 31 && mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].type != 32 && mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].type != 33)
            {

                vec gradient_face;
                vec a;
                a= mesh->Elements[ie].faces[jface].center - mesh->Elements[ie].center;
                double gradientn;
                double eface=0;
                for (int k = 0; k <mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].cellindex.size() ; ++k) {
                    if (mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].cellindex[k] == ie)
                    {
                        eface=ebound[iband*numDirection*mesh->Boundaries.size()*2 +inf*mesh->Boundaries.size()*2 + mesh->Elements[ie].faces[jface].bound * 2 + k];
                    }
                }
                gradientn= (eface-energyDensity[iband_local][inf_local][ie]) * dot_product(mesh->Elements[ie].faces[jface].norm, a) / dot_product(a, a);
                //
                vector<double> d1(mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors.size() + 1, 0);
                vector<vector<double>> J1(mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors.size() + 1);
                J1[0].push_back(mesh->Elements[ie].faces[jface].norm.x);
                J1[0].push_back(mesh->Elements[ie].faces[jface].norm.y);
                d1[0]=0;
                for (int m = 0; m < mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors.size(); m++) {
                    for (int k = 0; k < mesh->Boundaries[mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors[m]].cellindex.size() ; ++k) {
                        double ie1=mesh->Boundaries[mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors[m]].cellindex[k];
                        double iface1=mesh->Boundaries[mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors[m]].faceindex[k];
                        if (mesh->Elements[ie1].index==mesh->Elements[ie].index)
                        {
                            J1[m+1].push_back(mesh->Elements[ie1].faces[iface1].center.x-mesh->Elements[ie].faces[jface].center.x);
                            J1[m+1].push_back(mesh->Elements[ie1].faces[iface1].center.y-mesh->Elements[ie].faces[jface].center.y);
                            d1[m+1]= ebound[iband*numDirection*mesh->Boundaries.size()*2 +inf*mesh->Boundaries.size()*2 + mesh->Boundaries[mesh->Elements[ie].faces[jface].bound].neighbors[m] * 2 + k] - eface;
                        }
                    }

                }


                gradient_face=Get_gradient_2D(J1,d1);
                if(isnan(gradient_face.x)||isnan(gradient_face.y))
                {
                    cout<<1<<endl;
                }
                gradient_face.x=gradientn*mesh->Elements[ie].faces[jface].norm.x + gradient_face.x;
                gradient_face.y=gradientn*mesh->Elements[ie].faces[jface].norm.y + gradient_face.y;
                //gradient_face=gradient[i];
                //gradient_f2[i][j]=gradient_face;
                if(heatRatio[imatter][iband][inf]>0)
                {
                    ReMacroLocal[ie]-= modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * relaxationTime[matter[ie]][iband][inf] *
                            (groupVelocityX[matter[ie]][iband][inf] * gradient_face.x + groupVelocityY[matter[ie]][iband][inf] * gradient_face.y )
                            * (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                    ReMacroLocal[ie]+= kappaBulk[matter[ie]] * modeWeight[matter[ie]][iband][inf] * elementFaceArea[jface + ie * 6] * latticeRatio[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf]
                            * (elementFaceNormX[jface + ie * 6] * gradient_face.x + elementFaceNormY[jface + ie * 6] * gradient_face.y );
               }
            }
        }
    }
    cout<<"";
}*/

void StaticBTESynthetic::_set_bound_ee_1() const
{
    for (int inf_local = 0; inf_local < numDirectionLocal ; inf_local++)
    {
        int inf = ((inf_local) * numProc + worldRank) % numDirection;
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local)
        {
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            for (int ib = 0; ib < numBound; ++ib)
            {
                for (int icell = 0; icell < 2; ++icell) {
                    int ie =boundaryCell[ib][icell];
                    int jface =boundaryFace[ib][icell];
                    if(ie>=0)
                    {
                        if (heatRatio[matter[ie]][iband][inf] != 0)
                        {
                            double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf]
                                    * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                            if (dotproduct < 0)
                            {
                                if (boundaryType[ib]==1)
                                {
                                    double e = heatCapacity[matter[ie]][iband][inf] * boundaryThermal[ib];
                                    // mesh->Boundaries[ie].ee_bound[jface]->set(inf,iband,e);
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = e;
                                    //cout<<eboundLocal[iband_local * numBound * 2 + ib * 2 +icell]<<endl;
                                }
                                else if (boundaryType[ib]  == 2)
                                {
                                    double einsum1 = 0;
                                    double temp1 = 0;
                                    for (int nft = 0; nft < numDirection; ++nft)
                                    {
                                        double dotproduct1=(groupVelocityX[matter[ie]][iband][nft] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][nft]
                                                * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][nft] * elementFaceNormZ[jface + ie * 6]);
                                        if (dotproduct1 >= 0)
                                        {
                                            einsum1 += ebound[iband * numDirection * numBound * 2 + nft * numBound * 2 + ib * 2 + icell] *
                                                    (dotproduct1* modeWeight[matter[ie]][iband][nft]);
                                            temp1 += (dotproduct1*modeWeight[matter[ie]][iband][nft]);
                                        }
                                    }
                                    double e = einsum1 / temp1;
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = e;
                                }
                                else if (boundaryType[ib]  == 3)
                                {
                                    vec Reflectr;
                                    double dotproduct1=(directionX[inf] * elementFaceNormX[jface + ie * 6] + directionY[inf]
                                            * elementFaceNormY[jface + ie * 6] + directionZ[inf] * elementFaceNormZ[jface + ie * 6]);
                                    double ReflectrX=directionX[inf]- elementFaceNormX[jface + ie * 6] * dotproduct1 * 2;
                                    double ReflectrY=directionY[inf]- elementFaceNormY[jface + ie * 6] * dotproduct1 * 2;
                                    double ReflectrZ=directionZ[inf]- elementFaceNormZ[jface + ie * 6] * dotproduct1 * 2;
                                    //Reflectr = angles->direction[inf] - mesh->Elements[ie].faces[jface].norm * (angles->direction[inf] * mesh->Elements[ie].faces[jface].norm) * 2;
                                    double close = 1;
                                    int nf = -1;
                                    for (int k = 0; k < numDirection; ++k)
                                    {
                                        double length= sqrt(pow(ReflectrX-directionX[k],2)+pow(ReflectrY-directionY[k],2)+pow(ReflectrZ-directionZ[k],2));
                                        if (length < close)
                                        {
                                            nf = k;
                                            close = length;
                                        }
                                    }
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 + icell] = ebound[iband * numDirection * numBound * 2 + nf * numBound * 2 + ib * 2 +icell];

                                }
                                else if(boundaryType[ib]==5)
                                {
                                    double einsum1=0;
                                    double temp1=0;
                                    double einsum2=0;
                                    double temp2=0;
                                    int ie1=boundaryCell[ib][1-icell];
                                    int jface1=boundaryFace[ib][1-icell];
                                    for (int nft = 0; nft <numDirection ; ++nft) {
                                        double groupvelocity1= sqrt(pow(groupVelocityX[matter[ie]][iband][nft],2)+pow(groupVelocityY[matter[ie]][iband][nft],2)+pow(groupVelocityZ[matter[ie]][iband][nft],2));
                                        double groupvelocity2= sqrt(pow(groupVelocityX[matter[ie1]][iband][nft],2)+pow(groupVelocityY[matter[ie1]][iband][nft],2)+pow(groupVelocityZ[matter[ie1]][iband][nft],2));

                                        double dotproduct1=(directionX[nft] * elementFaceNormX[jface + ie * 6] + directionY[nft]
                                                * elementFaceNormY[jface + ie * 6] + directionZ[nft] * elementFaceNormZ[jface + ie * 6]);
                                        //double dotproduct1=(groupVelocityX[matter[ie]][iband][nft] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][nft]
                                         //* elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][nft] * elementFaceNormZ[jface + ie * 6]);
                                        if (dotproduct1 >= 0)
                                        {
                                            einsum1 +=groupvelocity1 * ebound[iband * numDirection * numBound * 2 + nft * numBound * 2 + ib * 2 + icell] *
                                                    (dotproduct1* modeWeight[matter[ie]][iband][nft]);
                                            //cout<<einsum1<<endl;
                                            temp1 += (dotproduct1*modeWeight[matter[ie]][iband][nft]);

                                        }
                                        else
                                        {
                                            einsum2+=groupvelocity2 * ebound[iband * numDirection * numBound * 2 + nft * numBound * 2 + ib * 2 + 1-icell] *
                                                    (dotproduct1* modeWeight[matter[ie]][iband][nft]);
                                            temp2+=(dotproduct1* modeWeight[matter[ie]][iband][nft]);
                                        }
                                    }
                                    double groupvelocity1= sqrt(pow(groupVelocityX[matter[ie]][iband][inf],2)+pow(groupVelocityY[matter[ie]][iband][inf],2)+pow(groupVelocityZ[matter[ie]][iband][inf],2));

                                    double e=einsum1/temp1*(1- boundaryTrans[ib][iband][icell])+einsum2/temp2*(1-boundaryTrans[ib][iband][icell]);
                                    //mesh->Boundaries[i].macro_Temp[j].Temp+=bands->bands[iband].lattice_ratio[mesh->Elements[ie].matter] * e / bands->bands[iband].heat_capacity[mesh->Elements[ie].matter] * angles->weight[inf] ;
                                    //mesh->Boundaries[i].ee_bound[j]->set(inf,iband,e);
                                    //cout<<e<<endl;
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 +icell]=e/groupvelocity1;
                                }
                                else if (boundaryType[ib] < 0)
                                {
                                    double e = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + boundaryConnect[ib] * 2 +icell];
                                    e = e + heatCapacity[matter[ie]][iband][inf] * (boundaryThermal[ib] - boundaryThermal[boundaryConnect[ib]]);
                                    // mesh->Boundaries[ie].ee_bound[jface]->set(inf,iband,e);
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = e;
                                }
                                else
                                {
                                    eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 +icell];
                                }
                            }
                            else
                            {
                                eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 +icell];
                            }
                        }
                        else
                        {
                            eboundLocal[iband_local * numBound * 2 + ib * 2 +icell] = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 +icell];
                        }
                    }

                }


            }
            //cout<<endl;
            for (int i = 0; i < numBound*2; ++i) {
               // cout<<eboundLocal[iband_local*numBound*2+i]<<endl;
            }
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                          numBound * 2,
                          MPI_DOUBLE,
                          (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                          numBound * 2,
                          MPI_DOUBLE,
                          MPI_COMM_WORLD);
        }
    }

}

void StaticBTESynthetic::_set_macro_boundary_condition() const
{
    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {
            boundaryTempMacro[ib][icell]=0;
            double sum = 0;
            int ie =boundaryCell[ib][icell];
            int jface =boundaryFace[ib][icell];
            if(ie>=0)
            {
                for (int iband = 0; iband < numBand; ++iband)
                {
                    for (int inf = 0; inf < numDirection; inf++)
                    {
                        double e = ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 +icell];
                        boundaryTempMacro [ib][icell] += latticeRatio[matter[ie]][iband][inf] * e / heatCapacity[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf];
                    }
                }
            }

        }

    }
    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {

            int ie =boundaryCell[ib][icell];
            int jface =boundaryFace[ib][icell];
            if(ie>=0)
            {
                if(boundaryType[ib]<0&&boundaryThermal[ib]>0)
                {
                    boundaryTempMacro[ib][icell]=boundaryThermal[ib]/2;
                    boundaryTempMacro[boundaryConnect[ib]][icell]=-boundaryThermal[ib]/2;
                }
            }
        }
    }
}

void StaticBTESynthetic::_set_vertex_temperature() const {


    double *temp=new double [numNode*numofMatter];
    double *R_n=new double [numNode*numofMatter];

    for (int i = 0; i < numNode*numofMatter ; ++i) {
        temp[i]=0;
        R_n[i]=0;
    }

    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            double distance= sqrt(pow(nodeX[elementVertexes[j+i*8]]-elementCenterX[i],2)+
                    pow(nodeY[elementVertexes[j+i*8]]-elementCenterY[i],2)+
                    pow(nodeZ[elementVertexes[j+i*8]]-elementCenterZ[i],2));
            temp[elementVertexes[j+i*8]*numofMatter+matter[i]]+=temperature[i]/
                    distance;
            R_n[elementVertexes[j+i*8]*numofMatter+matter[i]]+=1.0/distance;
        }
    }
    for (int i = 0; i < numNode*numofMatter ; ++i) {
        temperatureVertex[i]=temp[i]/R_n[i];
        //std::cout<<temperatureVertex[i]<<std::endl;
    }
    delete[] temp;
    delete[] R_n;


}

void  StaticBTESynthetic::_get_Re_macro()
{
    double *ReLocal=new double [numCellLocal];
    //cout<<"before while"<<endl;
    for (int icell_local = 0; icell_local < numCellLocal ; ++icell_local) {
        int ie= icell_local + startCell;
        //cout<<icell_local<<" "<<ie<<endl;
        ReLocal[icell_local]=ReMacro[ie];
        ReLocal[icell_local]-= elementHeatSource[ie] * elementVolume[ie];
        for (int jface = 0; jface < elementFaceSize[ie] ; ++jface) {
            double ax,ay,az;
            ax=0;ay=0;az=0;
            if (elementFaceBound[jface + ie * 6] >= 0)
            {
                if (boundaryType[elementFaceBound[jface + ie * 6]] != 3)
                {
                    for (int icell = 0; icell < 2; ++icell) {
                        if(boundaryCell[elementFaceBound[jface + ie * 6]][icell]==ie)
                        {
                            double temp = kappaBulk[matter[ie]] * elementFaceArea[jface + ie * 6] ;
                            ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                            ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                            az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                            temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                            ReLocal[icell_local] += -temp * boundaryTempMacro[elementFaceBound[jface + ie * 6]][icell];
                        }
                    }

                }
            } else
            {
                ax= elementCenterX[elementFaceNeighobr[jface + ie * 6]] - elementCenterX[ie];
                ay= elementCenterY[elementFaceNeighobr[jface + ie * 6]] - elementCenterY[ie];
                az= elementCenterZ[elementFaceNeighobr[jface + ie * 6]] - elementCenterZ[ie];
            }
            double a2x= elementFaceNormX[jface + ie * 6] - ax / (ax * elementFaceNormX[jface + ie * 6] + ay * elementFaceNormY[jface + ie * 6] + az * elementFaceNormZ[jface + ie * 6]);
            double a2y= elementFaceNormY[jface + ie * 6] - ay / (ax * elementFaceNormX[jface + ie * 6] + ay * elementFaceNormY[jface + ie * 6] + az * elementFaceNormZ[jface + ie * 6]);
            double a2z= elementFaceNormZ[jface + ie * 6] - az / (ax * elementFaceNormX[jface + ie * 6] + ay * elementFaceNormY[jface + ie * 6] + az * elementFaceNormZ[jface + ie * 6]);
            if (dimension == 2)
            {
                double gradient_X= nodeX[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeX[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                double gradient_Y= nodeY[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]] - nodeY[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]];
                ax= (temperatureVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] - temperatureVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]]) / (gradient_X * gradient_X + gradient_Y * gradient_Y) * gradient_X;
                ay= (temperatureVertex[elementFaceVertexes[1 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]] - temperatureVertex[elementFaceVertexes[0 + jface * 4 + ie * 6 * 4]*numofMatter+matter[ie]]) / (gradient_X * gradient_X + gradient_Y * gradient_Y) * gradient_Y;
            } else if (dimension == 3)
            {
                std::vector<double> d1(elementFaceVertexesSize[jface + ie * 6], 0);

                int ifa=elementFaceIndex[jface + ie * 6];

                int i1 = faceCell[ifa][0];
                int j1 = faceFace[ifa][0];

                    d1[0] = 0;
                    for (int m = 0; m <  elementFaceVertexesSize[j1+i1*6]-1; m++)
                    {
                        d1[m + 1] = 1.0/L_x * (temperatureVertex[elementFaceVertexes[m+1+j1*4+i1*6*4]*numofMatter+matter[i1]] - temperatureVertex[elementFaceVertexes[m+0+j1*4+i1*6*4]*numofMatter+matter[i1]] );
                    }

                    for (int m = 0; m < d1.size(); ++m)
                    {
                        ax+= FaceMatrix[ifa][0][m] * d1[m];
                        ay+= FaceMatrix[ifa][1][m] * d1[m];
                        az+= FaceMatrix[ifa][2][m] * d1[m];
                    }
                    //cout<<a.x<<" "<<a.y<<" "<<a.z<<endl;
                    //cout<<ax<<" "<<ay<<" "<<az<<endl;

                    //ax=a.x;
                    //ay=a.y;
                    //az=a.z;



            }
            if (elementFaceBound[jface + ie * 6] >= 0) {
                if (boundaryType[elementFaceBound[jface + ie * 6]] != 3) {
                    ReLocal[icell_local]-= kappaBulk[matter[ie]] * elementFaceArea[jface + ie * 6] * (a2x * ax + a2y * ay + a2z * az);
                }
            } else
            {
                ReLocal[icell_local] -= kappaBulk [matter[ie]]* elementFaceArea[jface + ie * 6] * (a2x * ax + a2y * ay + a2z * az);
            }

        }
    }
    MPI_Allgatherv(ReLocal,
                  numCellLocal,
                  MPI_DOUBLE,
                  Re,
                  numCellLocalList,
                  startCellList,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
    //cout<<"end gather"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < numCell; ++i) {
        Re1(i)=Re[i];
    }
    delete [] ReLocal;
}

/*void  StaticBTESynthetic::_get_Re_macro()
{
    //cout<<"before while"<<endl;
    for (int i= 0; i < numCell ; ++i) {

        //cout<<icell_local<<" "<<i<<endl;
        Re[i]=ReMacro[i];
        Re[i]-=elementHeatSource[i]*elementVolume[i];
        for (int j = 0; j <elementFaceSize[i] ; ++j) {
            double ax,ay,az;
            ax=0;ay=0;az=0;
            if (elementFaceBound[j+i*6]>=0)
            {
                if (boundaryType[elementFaceBound[j+i*6]]!=3)
                {
                    double temp = kappaBulk * elementFaceArea[j+i*6] ;
                    ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                    ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                    az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                    temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                    Re[i] += -temp * boundaryTempMacro[elementFaceBound[j + i * 6]];
                }
            } else
            {
                ax=elementCenterX[elementFaceNeighobr[j+i*6]]-elementCenterX[i];
                ay=elementCenterY[elementFaceNeighobr[j+i*6]]-elementCenterY[i];
                az=elementCenterZ[elementFaceNeighobr[j+i*6]]-elementCenterZ[i];
            }
            double a2x=elementFaceNormX[j+i*6]-ax/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            double a2y=elementFaceNormY[j+i*6]-ay/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            double a2z=elementFaceNormZ[j+i*6]-az/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            if (dimension == 2)
            {
                double gradient_X=nodeX[elementFaceVertexes[1+j*4+i*6*4]]-nodeX[elementFaceVertexes[0+j*4+i*6*4]];
                double gradient_Y=nodeY[elementFaceVertexes[1+j*4+i*6*4]]-nodeY[elementFaceVertexes[0+j*4+i*6*4]];
                ax=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradient_X*gradient_X+gradient_Y*gradient_Y)*gradient_X;
                ay=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradient_X*gradient_X+gradient_Y*gradient_Y)*gradient_Y;
            } else if (dimension == 3)
            {
                std::vector<double> d1(elementFaceVertexesSize[j+i*6], 0);

                int ifa=elementFaceIndex[j + i * 6];

                int i1 = faceCell[ifa];
                int j1 = faceFace[ifa];
                d1[0] = 0;
                for (int m = 0; m <  elementFaceVertexesSize[j1+i1*6]-1; m++)
                {
                    d1[m + 1] = 1e9 * (temperatureVertex[elementFaceVertexes[m+1+j1*4+i1*6*4]] - temperatureVertex[elementFaceVertexes[m+0+j1*4+i1*6*4]] );
                }

                for (int m = 0; m < d1.size(); ++m)
                {
                    ax+= FaceMatrix[ifa][0][m] * d1[m];
                    ay+= FaceMatrix[ifa][1][m] * d1[m];
                    az+= FaceMatrix[ifa][2][m] * d1[m];
                }

            }
            if (elementFaceBound[j+i*6]>=0) {
                if (boundaryType[elementFaceBound[j + i * 6]] != 3) {
                    Re[i]-= kappaBulk  * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
                }
            } else
            {
                Re[i] -= kappaBulk  * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
            }

        }
    }



    for (int i = 0; i < numCell; ++i) {
        Re1(i)=Re[i];
    }
}*/

void StaticBTESynthetic::_recover_temperature(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        temperatureLocal[ie] += latticeRatio[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie] * modeWeight[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf];
    }
}

void StaticBTESynthetic::_get_total_energy(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        totalEnergyLocal[ie] += energyDensity[iband_local][inf_local][ie] * modeWeight[matter[ie]][iband][inf] / capacityBulk[matter[ie]];
        
    }
}

void StaticBTESynthetic::_get_heat_flux(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        heatFluxXLocal[ie] += groupVelocityX[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie] ;
        heatFluxYLocal[ie] += groupVelocityY[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie];
        heatFluxZLocal[ie] += groupVelocityZ[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie];
    }

}

/*void StaticBTESynthetic::_get_bound_flux() const {
    for (int ib = 0; ib < numBound; ++ib) {
        boundaryFlux[ib]=0;
        int ie = boundaryCell[ib];
        int jface = boundaryFace[ib];
        for (int iband = 0; iband < numBand; ++iband) {
            for (int inf = 0; inf < numDirection; ++inf) {
                boundaryFlux[ib]+= (groupVelocityX[matter[ie]][iband][inf]*elementFaceNormX[boundaryCell[ib]*6+boundaryFace[ib]]+
                        groupVelocityY[matter[ie]][iband][inf]*elementFaceNormY[boundaryCell[ib]*6+boundaryFace[ib]]+
                        groupVelocityZ[matter[ie]][iband][inf]*elementFaceNormZ[boundaryCell[ib]*6+boundaryFace[ib]])
                                *modeWeight[matter[ie]][iband][inf] * ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 ] ;

            }
        }
        //cout<<boundaryFlux[ib]<<endl;
    }
}

void StaticBTESynthetic::_get_bound_temp() const {
    for (int ib = 0; ib < numBound; ++ib) {
        boundaryTemp[ib]=0;
        int ie = boundaryCell[ib];
        int jface = boundaryFace[ib];
        for (int iband = 0; iband < numBand; ++iband) {
            for (int inf = 0; inf < numDirection; ++inf) {
                boundaryTemp[ib]+= latticeRatio[matter[ie]][iband][inf]
                        *modeWeight[matter[ie]][iband][inf] * ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 ]/ heatCapacity[matter[ie]][iband][inf] ;

            }
        }
        //cout<<boundaryFlux[ib]<<endl;
    }
}
*/

void StaticBTESynthetic::_get_bound_flux() const {

}

void StaticBTESynthetic::_get_bound_temp() const {

}

bool StaticBTESynthetic::_get_magin_check_error(int nt, double error_temp_limit, double error_flux_limit)
{
    error = 0;
    double Tempmax = 0;
    for (int i = 0; i < numCell; ++i)
    {
        if (Tempmax < abs(temperature[i]))
        {
            Tempmax = abs(temperature[i]);
        }
    }
    for (int i = 0; i < numCell; ++i)
    {
        error += sqrt(pow(temperature[i] - temperatureOld[i], 2));
    }
    error = error / numCell / Tempmax;
    error_heat=0;
    double heat_max = 0;
    for (int i = 0; i < numCell; ++i)
    {
        double heat = sqrt(pow(heatFluxXGlobal[i], 2) + pow(heatFluxYGlobal[i], 2) + pow(heatFluxZGlobal[i], 2));
        double heat1 = sqrt(pow(heatFluxXOld[i], 2) + pow(heatFluxYOld[i], 2) + pow(heatFluxZOld[i], 2));
        if (heat_max < heat)
        {
            heat_max = heat;
        }
        error_heat += abs(heat - heat1);

    }
    error_heat = error_heat / numCell / heat_max;
    if (worldRank == 0)
    {
        cout << "----------------------------------------------------------------------------------" << endl;
        cout << "Iteration #" << nt << "\t Margin per band per cell: Temperature:" << error << " Heat Flux:" << error_heat << endl;
    }
    if ((error > errorOld&&error_heat>errorOld_heat)|| isnan(error)|| isnan(error_heat) || isinf(error_heat) || isinf(error))
    {
        errorIncreaseTime++;
    }
    else
    {
        errorIncreaseTime = 0;
    }
    errorOld = error;
    errorOld_heat=error_heat;
    return (error <= error_temp_limit && error_heat <= error_flux_limit);
}

void StaticBTESynthetic::_print_out() const
{
    ofstream output("TempLattice.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << temperature[i] << endl;
    }
    output.close();
    ofstream output1("Temperature.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output1 << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << totalEnergy[i] << endl;
    }
    output1.close();
    /*ofstream output1("Tempcell1.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output1 << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << temperature1[i] << endl;
    }
    output1.close();*/
    /*ofstream output2("HeatSource.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output2 << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << elementHeatSource[i] << endl;
    }
    output2.close();*/
    ofstream outputheat("HeatFlux.dat");
    for (int i = 0; i < numCell; ++i)
    {
        outputheat << setw(15) << elementCenterX[i]  << setw(15) << elementCenterY[i] << setw(15) << elementCenterZ[i] << setw(15) << heatFluxXGlobal[i] << setw(15) << heatFluxYGlobal[i] << setw(15) << heatFluxZGlobal[i] << endl;
    }
}

void StaticBTESynthetic::copy() const
{
    for (int i = 0; i < numCell; ++i)
    {
        temperatureOld[i] = temperature[i];
        heatFluxXOld[i] = heatFluxXGlobal[i];
        heatFluxYOld[i] = heatFluxYGlobal[i];
        heatFluxZOld[i] = heatFluxZGlobal[i];
    }
    for (int i = 0; i < numCell; ++i)
    {
        totalEnergyLocal[i]=0;
        temperatureLocal[i] = 0;
        heatFluxXLocal[i] = 0;
        heatFluxYLocal[i] = 0;
        heatFluxZLocal[i] = 0;
        temperature[i] = 0;
        heatFluxXGlobal[i] = 0;
        heatFluxYGlobal[i] = 0;
        heatFluxZGlobal[i] = 0;
        ReMacroLocal[i]=0;

    }
    for (int i = 0; i < numBandLocal * numBound * 2; ++i)
    {
        eboundLocal[i] = 0;
    }
}

void StaticBTESynthetic::_marco_solution()
{
    double internalerror=1000;
    intererrorIncreaseTime=0;
    while(internalerror>1e-6)
    {
        for (int i = 0; i < numCell; ++i)
        {
            temperature1[i]=temperature[i];
        }
        _set_vertex_temperature();

        _get_Re_macro();
        sol1=cg_macro.solve(Re1);
        //double *sol = _solve_matrix_hypre((int*)csrRowPtr_macro, (int*)csrColInd_macro, csrVal_macro, Re,numCell,worldRank,numProc);
        internalerror=0;
        for (int i = 0; i < numCell; ++i)
        {
            temperature[i]=sol1(i);
            internalerror+=pow(temperature1[i]-temperature[i],2);
        }
        internalerror= internalerror / numCell;
        if(internalerror>error_old)
        {
            intererrorIncreaseTime++;
        }
        error_old=internalerror;
        if(intererrorIncreaseTime>=10)
        {
            cout<<"cannot congerve, try limiter or 1 st order"<<endl;
            exit(1);
        }
        if(worldRank==0)
        {
            cout<<"internal error: "<<internalerror<<endl;
        }


        MPI_Barrier(MPI_COMM_WORLD);
        //delete [] sol;

    }
};

void StaticBTESynthetic::_get_coefficient_Iterative(int iband_local, int inf_local)
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    std::vector<Tri> tripletList;

    vector<pair<int, double>> compressed_Ke;

    for (int ie = 0; ie < numCell; ++ie)
    {
        compressed_Ke.resize(0);
        double valueindia = 0;
        //int matter = mesh->Elements[ie].matter;
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface)
        {
            double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                    + groupVelocityY[matter[ie]][iband][inf]* elementFaceNormY[jface + ie * 6]
                    + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
            double temp = relaxationTime[matter[ie]][iband][inf] *
                    elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
            if (dotproduct >= 0)
            {
                valueindia += temp;
            }
            else
            {
                if (elementFaceBound[jface + ie * 6] == -1)
                {
                    int neiindex = elementFaceNeighobr[jface + ie * 6];

                    compressed_Ke.emplace_back(neiindex, temp);
                }
                else
                {
                    //int boundindex = elementFaceBound[jface+ie*6];
                    //if (boundaryType[boundindex] == 0)
                    //{
                    //    int neiindex =  elementFaceNeighobr[jface+ie*6];
                    //    compressed_Ke.emplace_back(neiindex, temp);
                    //}
                }
            }
        }
        valueindia += 1.0;
        compressed_Ke.emplace_back(ie, valueindia);
        sort(compressed_Ke.begin(), compressed_Ke.end());

        for (int j = 0; j < compressed_Ke.size(); ++j) {
            tripletList.push_back(Tri(ie, compressed_Ke[j].first, compressed_Ke[j].second));
        }
    }
    stiffMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    //cg[iband_local][inf_local].analyzePattern(stiffMatrix);
    //cg[iband_local][inf_local].factorize(stiffMatrix);
    cout<<"finish calculate coefficients: "<<iband_local<<" "<<inf_local<<endl;
};

void StaticBTESynthetic::_get_coefficient_macro_Iterative() const {
    int csrRowPtr_iter = 0, csrColInd_iter = 0, csrVal_iter = 0;
    csrRowPtr_macro[0] = 0;
    csrRowPtr_iter++;

    for (int ie = 0; ie < numCell; ++ie)
    {
        vector<pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(elementFaceSize[ie] + 1);
        double valueindia = 0;
        double ax,ay,az;
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface)
        {
            double temp = kappaBulk [matter[ie]] * elementFaceArea[jface + ie * 6] ;
            if (elementFaceBound[jface + ie * 6] == -1)
            {
                double ax= elementCenterX[elementFaceNeighobr[jface + ie * 6]] - elementCenterX[ie];
                double ay= elementCenterY[elementFaceNeighobr[jface + ie * 6]] - elementCenterY[ie];
                double az= elementCenterZ[elementFaceNeighobr[jface + ie * 6]] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[jface + ie * 6], temp);

            }
            else if (boundaryType[elementFaceBound[jface + ie * 6]] != 3)
            {
                double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
            }
        }
        compressed_Ke.emplace_back(ie, valueindia);
        sort(compressed_Ke.begin(), compressed_Ke.end());
        for (int j = 0; j < compressed_Ke.size(); ++j) {
            csrColInd_macro[csrColInd_iter++] = compressed_Ke[j].first;
            //cout<< compressed_Ke[num].first<<" ";
            csrVal_macro[csrVal_iter++] = compressed_Ke[j].second;
        }
        csrRowPtr_macro[csrRowPtr_iter] = csrRowPtr_macro[csrRowPtr_iter - 1] + compressed_Ke.size();
        csrRowPtr_iter++;
    }
}

void StaticBTESynthetic::_get_Re_firstorder(int iband_local, int inf_local) {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    for (int i = 0; i < numCell; ++i) {
        Re[i]=0;
    }

    //internal
    for (int ie = 0; ie < numCell; ++ie) {

        Re[ie] += temperatureOld[ie] * heatCapacity[matter[ie]][iband][inf];
        Re[ie] += elementHeatSource[ie] * heatRatio[matter[ie]][iband][inf] *
                relaxationTime[matter[ie]][iband][inf];

    }

    for (int ib = 0; ib < numBound; ++ib) {
        for (int icell = 0; icell < 2; ++icell) {
            int ie = boundaryCell[ib][icell];
            int jface = boundaryFace[ib][icell];
            if(ie >= 0)
            {
                double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                        + groupVelocityY[matter[ie]][iband][inf]* elementFaceNormY[jface + ie * 6]
                        + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);

                if (dotproduct < 0) {
                    double temp = relaxationTime[matter[ie]][iband][inf] *
                            elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
                    Re[ie] -= temp * ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell];
                }
            }
        }

    }

    for (int i = 0; i < numCell; ++i) {
        Re1(i)=Re[i];
    }
}

void StaticBTESynthetic::_get_bound_ee_firstorder(int iband_local, int inf_local) const
{
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ib = 0; ib < numBound; ++ib)
    {
        for (int icell = 0; icell < 2; ++icell) {
            int ie = boundaryCell[ib][icell];
            int jface = boundaryFace[ib][icell];
            if(ie>=0)
            {
                double dotproduct=(groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6] + groupVelocityY[matter[ie]][iband][inf]
                        * elementFaceNormY[jface + ie * 6] + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                if (dotproduct>= 0)
                {
                    double e = (energyDensity[iband_local][inf_local][ie] );
                    eboundLocal[iband_local * numBound * 2 + ib * 2 + icell] = e;
                    ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell] = e;
                }
            }
        }

    }
}







