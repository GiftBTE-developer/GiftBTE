//
// Created by yuehu on 2023/5/21.
//

#include "TransientBTE/transient.h"
#include <algorithm>
#include <iomanip>
#include <chrono>

Transient::~Transient() {

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

    cout << "~StaticBTESynthetic is activated !!" << endl;
#ifdef USE_GPU

#else

#endif

}
Transient::Transient(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc,
                     int world_rank,double deltaT,double totalT) {
    this->deltaT=deltaT;
    this->totalT=totalT;
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
            if(world_rank==0)
                cout<<numCellLocalList[i]<<endl;
        } else
        {
            numCellLocalList[i]=  numCell / num_proc;
            if(world_rank==0)
                cout<<numCellLocalList[i]<<endl;
        }
    }
    startCellList=new int [num_proc];
    for (int i = 0; i < num_proc; ++i) {
        if(i<numCell % num_proc)
        {
            startCellList[i]=  (numCell / num_proc+1)*i;
            if(world_rank==0)
                cout<<startCellList[i]<<endl;
        } else
        {
            startCellList[i]=(numCell / num_proc+1)*(numCell % num_proc)+(numCell / num_proc)*(i-numCell % num_proc);
            if(world_rank==0)
                cout<<startCellList[i]<<endl;
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
        capacityBulk[i]=bands->kappabulk[i];
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


void Transient::_set_initial(int Use_Backup) const {
    for (int i = 0; i < numCell; ++i) {
        temperature[i]= cos(2*PI/L_x*elementCenterX[i]);
    }
    for (int iband_local = 0; iband_local < numBandLocal; iband_local++)
    {
        for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++)
        {
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            for (int kk = 0; kk < numCell; kk++)
            {
                energyDensity[iband_local][inf_local][kk]= temperature[kk] * heatCapacity[matter[kk]][iband][inf];
            }
        }

    }
}

void Transient::_get_bound_ee(int iband_local, int inf_local) const {
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
                }
            }
        }

    }
    MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
}

void Transient::_set_bound_ee_1() const {
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

void Transient::_recover_temperature(int iband_local, int inf_local) const {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        temperatureLocal[ie] += latticeRatio[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie] * modeWeight[matter[ie]][iband][inf] / heatCapacity[matter[ie]][iband][inf];
    }
}

void Transient::_get_total_energy(int iband_local, int inf_local) const {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        totalEnergyLocal[ie] += energyDensity[iband_local][inf_local][ie] * modeWeight[matter[ie]][iband][inf] / capacityBulk[matter[ie]];
    }
}

void Transient::_get_heat_flux(int iband_local, int inf_local) const {
    int inf = ((inf_local) * numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
    for (int ie = 0; ie < numCell; ++ie)
    {
        heatFluxXLocal[ie] += groupVelocityX[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie] ;
        heatFluxYLocal[ie] += groupVelocityY[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie];
        heatFluxZLocal[ie] += groupVelocityZ[matter[ie]][iband][inf] * modeWeight[matter[ie]][iband][inf] * energyDensity[iband_local][inf_local][ie];
    }
}

bool Transient::_get_magin_check_error(int nt, double error_temp_limit, double error_flux_limit) {
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

void Transient::_print_out() const {
    ofstream output("Tempcell.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << temperature[i] << endl;
    }
    output.close();
    ofstream output1("Tempcell1.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output1 << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << temperature1[i] << endl;
    }
    output1.close();
    ofstream output2("HeatSource.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output2 << elementCenterX[i] << " " << elementCenterY[i] << " " << elementCenterZ[i] << " " << elementHeatSource[i] << endl;
    }
    output2.close();
    ofstream outputheat("HeatFlux.dat");
    for (int i = 0; i < numCell; ++i)
    {
        outputheat << setw(15) << elementCenterX[i]  << setw(15) << elementCenterY[i] << setw(15) << elementCenterZ[i] << setw(15) << heatFluxXGlobal[i] << setw(15) << heatFluxYGlobal[i] << setw(15) << heatFluxZGlobal[i] << endl;
    }
}

void Transient::copy() const {
    for (int i = 0; i < numCell; ++i)
    {
        temperatureOld[i] = temperature[i];
        heatFluxXOld[i] = heatFluxXGlobal[i];
        heatFluxYOld[i] = heatFluxYGlobal[i];
        heatFluxZOld[i] = heatFluxZGlobal[i];
    }
    for (int i = 0; i < numCell; ++i)
    {
        temperatureLocal[i] = 0;
        heatFluxXLocal[i] = 0;
        heatFluxYLocal[i] = 0;
        heatFluxZLocal[i] = 0;
        temperature[i] = 0;
        heatFluxXGlobal[i] = 0;
        heatFluxYGlobal[i] = 0;
        heatFluxZGlobal[i] = 0;


    }
    for (int i = 0; i < numBandLocal * numBound * 2; ++i)
    {
        eboundLocal[i] = 0;
    }
}

void Transient::_delete_cell_matrix() const {

    for (int i = 0; i < numCell; ++i) {

        for (int j = 0; j < 3; ++j) {
            delete [] CellMatrix[i][j];
        }
        delete []  CellMatrix[i];
    }
    delete [] CellMatrix;
}

void Transient::_set_cell_matrix_larger() {
    CellMatrix=new double  **[numCell];
    for (int i = 0; i < numCell; ++i) {
        CellMatrix[i]=new double *[3];
        for (int j = 0; j < 3; ++j) {
            CellMatrix[i][j]=new double [elementNeighborList[i].size()];
        }
    }

}

void Transient::_get_CellMatrix_larger() const {
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

void Transient::_get_gradient_larger(int Use_limiter, int iband_local, int inf_local) const {
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


void Transient::_get_explicit_Re(int itime, int spatial_order, int Use_limiter,int iband_local, int inf_local,double deltaTime)
{
    int inf = ((inf_local)*numProc + worldRank) % numDirection;
    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;

    for (int i = 0; i < numCell; ++i) {
        Re[i] = 0;
    }
    //maxy要加在前面
    double max_y=0;
    for (int ie = 0; ie< numCell; ie++)
    {
        if (elementCenterY[ie] > max_y)
        {
            max_y = elementCenterY[ie];
        }
    }

    //internal
    for (int ie = 0; ie < numCell; ++ie) {
        for (int jface = 0; jface < elementFaceSize[ie]; ++jface) {
            double dotproduct = (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                    + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6]
                    + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
            double temp = elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;           //
            if (dotproduct >= 0) {
                double ax = elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                double ay = elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                double az = elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                double e = (energyDensity[iband_local][inf_local][ie] + (ax * gradientX[ie] + ay * gradientY[ie] + az * gradientZ[ie]) * limit[ie]);
                Re[ie] += temp * e;
                // if(ie==0)
                //cout<<ie<<"dot>0 "<<temp<<"  "<<e<<endl;
            }
            else {
                if (elementFaceBound[jface + ie * 6] == -1) {
                    int neiindex = elementFaceNeighobr[jface + ie * 6];
                    double ax = elementFaceCenterX[jface + ie * 6] - elementCenterX[neiindex];
                    double ay = elementFaceCenterY[jface + ie * 6] - elementCenterY[neiindex];
                    double az = elementFaceCenterZ[jface + ie * 6] - elementCenterZ[neiindex];
                    double e = (energyDensity[iband_local][inf_local][neiindex] + (ax * gradientX[neiindex] + ay * gradientY[neiindex] + az * gradientZ[neiindex]) * limit[neiindex]);
                    Re[ie] += temp * e;
                    //if(ie==0)
                    //cout<<ie<<"dot<0 "<<temp<<"  "<<e<<endl;
                }
            }
        }
        //equlibrium
        Re[ie] -= temperatureOld[ie] * heatCapacity[matter[ie]][iband][inf]/relaxationTime[matter[ie]][iband][inf];
        //if(ie==0)  cout<<ie<<"  req "<<temperatureOld[ie]* heatCapacity[matter[ie]][iband][inf]/relaxationTime[matter[ie]][iband][inf]<<endl;
        //heatsource
        //Re[ie] -= elementHeatSource[ie] * heatRatio[matter[ie]][iband][inf];

        //TDTR_heatsorce
        /*double rr = pow(elementCenterX[ie] ,2) + pow(elementCenterY[ie]  - 1e-5,2);
        double RR = 2e-6;
        double heatindex = 0;
        double heatratio = exp(-2*rr/pow(RR,2));
        // if (elementFaceCenterY[ie] == max_y) cout<<ie<<"  "<<elementFaceCenterX[ie]<<"  "<<elementFaceCenterY[ie]<< "  "<<max_y<<"   "<< pow(rr, 0.5)<<"   "<<2*RR <<"  "<<elementHeatSource[ie] <<endl;
        //cout<< pow(rr, 0.5)<<"  "<<2*RR<<endl;
        if (pow(rr, 0.5) <= 2*RR)
        {
            if (elementCenterY[ie] >9e-6)
            {
                rr=pow(elementCenterX[ie] ,2) + pow(elementCenterY[ie]- 1e-5,2);
                //cout<<"ie"<<ie<< "  "<<elementCenterY[ie]<<"   "<<max_y<<endl;
                int times=(1e-12)/(deltaTime);
                int new_itime=itime/times;

                int tt=12500;
                int numheat = (new_itime)/tt;
                int checkheat = (new_itime)-numheat*tt;
                double interg = 1.0;
                double TT=1.0/(tt*8*deltaTime);//2e8;
                if (checkheat>1)
                {
                    heatindex = 0;
                }
                else
                {
                    //double deltaTime=1e-14;
                    interg = interg+(double)(new_itime);
                    heatindex = sin(interg*deltaTime*2*PI*TT)+1;
                }
                //cout<<"sin   "<<sin(interg*deltaTime*2*PI*TT)<<"   "<<interg*deltaTime*2*PI*TT<<" deltaTime  "<<deltaTime   <<endl;
                //cout<<"ie   "<<ie<<"    itime  "<<itime<<"inf"<<inf<<"iband"<<iband<<"checkheat:   "<<checkheat<<"   numheat :"<<numheat<<" heatindex: "<<heatindex<<endl;
                //cout<<ie<<" heatindex  "<<heatindex<<"  heatratio  "<<heatratio<<" elementHeatSource "<<elementHeatSource[ie]<<"  heatRatio "<<heatRatio[matter[ie]][iband][inf]<<"modeWeight[matter[ie]][iband][inf] "<<modeWeight[matter[ie]][iband][inf] <<endl;
                //cout<<"inf"<<inf<<"iband"<<iband<<"re  "<<heatindex * heatratio  * heatRatio[matter[ie]][iband][inf] *modeWeight[matter[ie]][iband][inf] * elementHeatSource[ie]<<endl;
            }
        }
        double h =heatindex * heatratio  * heatRatio[matter[ie]][iband][inf]    * elementHeatSource[ie];//*  modeWeight[matter[ie]][iband][inf]
        //cout<<h<<"  heatratio  "<<heatratio<<endl;
        // heatsources[ie] = heatindex * heatratio * elementHeatSource[ie]  * elementVolume[ie];
        if (ie==309){
            //cout<<ie<<"   "<<heatindex<<"   "<<heatratio<<"  "<<elementHeatSource[ie]<<"   "<<heatRatio[matter[ie]][iband][inf]<<endl;
        }
        Re[ie] -= h;*/

    }
    for (int ib = 0; ib < numBound; ++ib) {
        for (int icell = 0; icell < 2; ++icell) {
            int ie = boundaryCell[ib][icell];
            int jface = boundaryFace[ib][icell];
            //cout<<"  ie "<<ie<<"  dot "<<endl;
            if (ie >= 0)
            {
                double dotproduct = (groupVelocityX[matter[ie]][iband][inf] * elementFaceNormX[jface + ie * 6]
                        + groupVelocityY[matter[ie]][iband][inf] * elementFaceNormY[jface + ie * 6]
                        + groupVelocityZ[matter[ie]][iband][inf] * elementFaceNormZ[jface + ie * 6]);
                // cout<<"  ie "<<ie<<"  dot "<<dotproduct<<"  groupVelocityX[matter[ie]][iband][inf]"<<groupVelocityX[matter[ie]][iband][inf]<<"iband local"<<iband<<"  "<<iband_local<<"inf local"<<inf<<"  "<<inf_local<<" elementFaceNormX[jface + ie * 6] "<<elementFaceNormX[jface + ie * 6]<<endl;
                if (dotproduct < 0) {
                    double temp = elementFaceArea[jface + ie * 6] / elementVolume[ie] * dotproduct;
                    Re[ie] += temp * ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell];
                    //   cout<<"eebound "<<ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell]<<endl;
                    // if(ie==0) cout<<ie<<" bound "<<temp<<"   "<<ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib * 2 + icell]<<endl;
                }
            }
        }

    }
}