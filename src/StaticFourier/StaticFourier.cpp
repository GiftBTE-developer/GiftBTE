//
// Created by huyue on 2023/1/17.
//

//
// Created by yuehu on 2022/12/12.
//

#include "StaticFourier/StaticFourier.h"
#include <algorithm>
#include <iomanip>
//#include <petscksp.h>
#include "utility/utility.h"
using namespace std;

StaticFourier::~StaticFourier()
{
   delete [] temperature;
   delete [] temperatureOld;
   delete [] temperatureVertex;
   delete [] kappa;

   delete [] elementFaceSize;

   delete [] elementFaceNeighobr;
   delete [] elementFaceArea;

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

   delete [] elementVolume;

   delete [] elementHeatSource;

   delete [] heatFluxXGlobal;
   delete [] heatFluxYGlobal;
   delete [] heatFluxZGlobal;

   delete [] elementCenterX;
   delete [] elementCenterY;
   delete [] elementCenterZ;

   delete []  boundaryType;
   delete [] boundaryTemp;
   delete [] boundaryThermal;
   delete []  boundaryFlux;
   delete [] boundaryCell;
   delete []  boundaryFace;
   delete [] boundaryResis;
   delete [] boundaryConnect;
   delete [] nodeX;
   delete [] nodeY;
   delete [] nodeZ;

   delete []   csrRowPtr;
   delete []   csrColInd;
   delete []  csrVal;
   delete [] Re;

   for (int i = 0; i < numCell; ++i) {

       for (int j = 0; j < 3; ++j) {
           delete [] CellMatrix[i][j];
       }
       delete []  CellMatrix[i];
   }
   delete [] CellMatrix;
}

StaticFourier::StaticFourier(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc, int world_rank)
{
    this->worldRank=world_rank;
    this->numProc=num_proc;
    L_x=mesh->L_x;
    L_y=mesh->L_y;
    L_z=mesh->L_z;
    dimension=mesh->Dimension;
    numCell=mesh->Elements.size();
    numnode=mesh->Nodes.size();
    temperature=new double [numCell];
    temperatureOld=new double [numCell];
    temperatureVertex=new double [mesh->Nodes.size()];
    kappa=new double [numCell];

    for (int i = 0; i < numCell ; ++i) {
        kappa[i]=bands->kappabulk[mesh->Elements[i].matter];
        //cout<<kappa[i]<<endl;
    }

    elementFaceSize=new int [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementFaceSize[i]=mesh->Elements[i].faces.size();
    }
    elementFaceNeighobr=new int [numCell * 6];
    elementFaceArea=new double [numCell * 6];
    for (int i = 0; i < numCell; ++i) {
        for (int j = 0; j < elementFaceSize[i] ; ++j) {
            elementFaceArea[j+i*6]=mesh->Elements[i].faces[j].area;
            elementFaceNeighobr[j+i*6]=mesh->Elements[i].faces[j].neighbor;
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


    heatFluxXGlobal = new double[numCell];

    heatFluxYGlobal = new double[numCell];

    heatFluxZGlobal = new double[numCell];

    elementCenterX=new double [numCell];
    elementCenterY=new double [numCell];
    elementCenterZ=new double [numCell];
    for (int i = 0; i < numCell ; ++i) {
        elementCenterX[i]=mesh->Elements[i].center.x;
        elementCenterY[i]=mesh->Elements[i].center.y;
        elementCenterZ[i]=mesh->Elements[i].center.z;
    }
    numBound=mesh->Boundaries.size();
    boundaryType=new int [numBound];
    boundaryTemp=new double [numBound];
    boundaryThermal=new double [numBound];
    boundaryFlux=new double [numBound];
    boundaryCell=new int [numBound];
    boundaryFace=new int [numBound];
    boundaryConnect=new int [numBound];

    for (int i = 0; i <numBound ; ++i) {
        boundaryConnect[i]=mesh->Boundaries[i].connection;
    }
    for (int i = 0; i <numBound ; ++i) {
        boundaryType[i]=mesh->Boundaries[i].type;
        boundaryThermal[i]=mesh->Boundaries[i].Temperature;
        boundaryCell[i]=mesh->Boundaries[i].cellindex[0];
        boundaryFace[i]=mesh->Boundaries[i].faceindex[0];

        boundaryFlux[i]=0;
    }
    boundaryResis=new double [numBound];
    for (int ib = 0; ib < numBound; ++ib) {
        int i=mesh->Boundaries[ib].cellindex[0];
        boundaryResis[ib]=0;
        if(boundaryType[ib]==5)
        {
            for (int iband = 0; iband < bands->nband; ++iband) {
                boundaryResis[ib]+=bands->bands[iband].heat_capacity[mesh->Elements[i].matter]
                        *bands->bands[iband].group_velocity[mesh->Elements[i].matter]*mesh->Boundaries[ib].trans[0][iband];
            }
        }
    }

    nodeX=new double [numnode];
    nodeY=new double [numnode];
    nodeZ=new double [numnode];
    for (int i = 0; i < numnode; ++i) {
        nodeX[i]=mesh->Nodes[i].x;
        nodeY[i]=mesh->Nodes[i].y;
        nodeZ[i]=mesh->Nodes[i].z;
    }

    csrRowPtr= new unsigned int[numCell + 1];
    csrColInd = new unsigned int[7 * numCell + 1];
    csrVal = new double[7 * numCell + 1];
    Re=new double [numCell];


    stiffMatrix.resize(numCell, numCell);
    Re1.resize(numCell);
    sol1.resize(numCell);
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
    CellMatrix=new double **[numCell];
    for (int i = 0; i < numCell; ++i) {
        CellMatrix[i]=new double *[3];

        for (int j = 0; j < 3; ++j) {
            CellMatrix[i][j]=new double [elementNeighborList[i].size()];
        }
    }
}

void  StaticFourier::_set_initial() {
    for (int i = 0; i < numCell ; ++i) {
        temperatureOld[i]=0;
        temperature[i]=0;
        Re[i]=0;
    }
}



void StaticFourier::_set_vertex_temperature() {


    double *temp=new double [numnode];
    double *R_n=new double [numnode];

    for (int i = 0; i <numnode ; ++i) {
        temp[i]=0;
        R_n[i]=0;
    }

    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            double distance= sqrt(pow(nodeX[elementVertexes[j+i*8]]-elementCenterX[i],2)+
                                  pow(nodeY[elementVertexes[j+i*8]]-elementCenterY[i],2)+
                                  pow(nodeZ[elementVertexes[j+i*8]]-elementCenterZ[i],2));
            temp[elementVertexes[j+i*8]]+=temperature[i]/
                    distance;
            R_n[elementVertexes[j+i*8]]+=1.0/distance;
        }
    }
    for (int i = 0; i <numnode ; ++i) {
        temperatureVertex[i]=temp[i]/R_n[i];
        //std::cout<<temperatureVertex[i]<<std::endl;
    }
    delete[] temp;
    delete[] R_n;


}
void StaticFourier::_get_CellMatrix() const
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

void  StaticFourier::_get_coefficient()
{
    int csrRowPtr_iter = 0, csrColInd_iter = 0, csrVal_iter = 0;
    csrRowPtr[0] = 0;
    csrRowPtr_iter++;
    for (int i = 0; i < numCell ; ++i) {
        std::vector<std::pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(6 + 1);
        double valueindia = 0;

        for (int j = 0; j <elementFaceSize[i] ; ++j) {


            if (elementFaceBound[j+i*6]==-1)
            {
                double temp = 2*(kappa[i]*kappa[elementFaceNeighobr[j+i*6]])/(kappa[i]+kappa[elementFaceNeighobr[j+i*6]]) * elementFaceArea[j+i*6] ;
                double ax=elementCenterX[elementFaceNeighobr[j+i*6]]-elementCenterX[i];
                double ay=elementCenterY[elementFaceNeighobr[j+i*6]]-elementCenterY[i];
                double az=elementCenterZ[elementFaceNeighobr[j+i*6]]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[j+i*6], temp);

            }
            else if (boundaryType[elementFaceBound[j+i*6]]==5)
            {
                double temp = 2*(kappa[i]*kappa[elementFaceNeighobr[j+i*6]])/(kappa[i]+kappa[elementFaceNeighobr[j+i*6]]) * elementFaceArea[j+i*6] ;
                double ax=elementCenterX[elementFaceNeighobr[j+i*6]]-elementCenterX[i];
                double ay=elementCenterY[elementFaceNeighobr[j+i*6]]-elementCenterY[i];
                double az=elementCenterZ[elementFaceNeighobr[j+i*6]]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6] * az);
                temp=1.0/(1.0/temp+1.0/(boundaryResis[elementFaceBound[j+i*6]]*elementFaceArea[j+i*6]));
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[j+i*6], temp);
            }
            else if (boundaryType[elementFaceBound[j+i*6]]==1)
            {
                double temp = kappa[i] * elementFaceArea[j+i*6] ;
                double ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                double ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                double az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
            }
            else if(boundaryType[elementFaceBound[j+i*6]]<0)
            {
                if(boundaryThermal[elementFaceBound[j+i*6]]!=0||boundaryThermal[boundaryConnect[elementFaceBound[j+i*6]]]!=0)
                {
                    double temp = kappa[i] * elementFaceArea[j+i*6] ;
                    double ax= elementFaceCenterX[j+i*6] - elementCenterX[i];
                    double ay= elementFaceCenterY[j+i*6] - elementCenterY[i];
                    double az= elementFaceCenterZ[j+i*6] - elementCenterZ[i];
                    temp=temp/(elementFaceNormX[j+i*6] * ax + elementFaceNormY[j+i*6] * ay + elementFaceNormZ[j+i*6] * az);
                    //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                    valueindia += -temp;
                }
            }

        }
        compressed_Ke.emplace_back(i,valueindia);
        sort(compressed_Ke.begin(), compressed_Ke.end());
        for (int num=0;num<compressed_Ke.size();num++) {
            csrColInd[csrColInd_iter++] = compressed_Ke[num].first;
            //cout<< compressed_Ke[num].first<<" ";
            csrVal[csrVal_iter++] = compressed_Ke[num].second;
            if(isnan(compressed_Ke[num].first)||isnan(compressed_Ke[num].second)||isinf(compressed_Ke[num].first)||isinf(compressed_Ke[num].second))
                cout<<compressed_Ke[num].first<<" "<<compressed_Ke[num].second<<endl;
            //cout<<compressed_Ke[num].second<<endl;
        }

        csrRowPtr[csrRowPtr_iter] = csrRowPtr[csrRowPtr_iter - 1] + compressed_Ke.size();

        csrRowPtr_iter++;
    }
}

void StaticFourier::_get_coefficient_macro()
{
    std::vector<Tri> tripletList;
    for (int ie = 0; ie < numCell; ++ie)
    {
        vector<pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(elementFaceSize[ie] + 1);
        double valueindia = 0;

        for (int jface = 0; jface < elementFaceSize[ie]; ++jface)
        {
           if (elementFaceBound[jface + ie * 6] == -1)
            {
                double temp = 2*(kappa[ie]*kappa[elementFaceNeighobr[jface+ie*6]])/(kappa[ie]+kappa[elementFaceNeighobr[jface+ie*6]]) * elementFaceArea[jface+ie*6] ;

                double ax= elementCenterX[elementFaceNeighobr[jface + ie * 6]] - elementCenterX[ie];
                double ay= elementCenterY[elementFaceNeighobr[jface + ie * 6]] - elementCenterY[ie];
                double az= elementCenterZ[elementFaceNeighobr[jface + ie * 6]] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[jface + ie * 6], temp);

            }
           else if (boundaryType[elementFaceBound[jface+ie*6]]==5)
           {
               double temp = 2*(kappa[ie]*kappa[elementFaceNeighobr[jface+ie*6]])/(kappa[ie]+kappa[elementFaceNeighobr[jface+ie*6]]) * elementFaceArea[jface+ie*6] ;
               double ax=elementCenterX[elementFaceNeighobr[jface+ie*6]]-elementCenterX[ie];
               double ay=elementCenterY[elementFaceNeighobr[jface+ie*6]]-elementCenterY[ie];
               double az=elementCenterZ[elementFaceNeighobr[jface+ie*6]]-elementCenterZ[ie];
               temp=temp/(elementFaceNormX[jface+ie*6] * ax+elementFaceNormY[jface+ie*6]* ay+elementFaceNormZ[jface+ie*6] * az);
               temp=1.0/(1.0/temp+1.0/(boundaryResis[elementFaceBound[jface+ie*6]]*elementFaceArea[jface+ie*6]));
               //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
               valueindia += -temp;
               compressed_Ke.emplace_back(elementFaceNeighobr[jface+ie*6], temp);
           }
           else if (boundaryType[elementFaceBound[jface+ie*6]]==1)
           {
               double temp = kappa[ie] * elementFaceArea[jface+ie*6] ;
                double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
            }
           else if(boundaryType[elementFaceBound[jface+ie*6]]<0)
           {
               if(boundaryThermal[elementFaceBound[jface+ie*6]]!=0||boundaryThermal[boundaryConnect[elementFaceBound[jface+ie*6]]]!=0)
               {
                   double temp = kappa[ie] * elementFaceArea[jface+ie*6] ;
                   double ax= elementFaceCenterX[jface + ie * 6] - elementCenterX[ie];
                   double ay= elementFaceCenterY[jface + ie * 6] - elementCenterY[ie];
                   double az= elementFaceCenterZ[jface + ie * 6] - elementCenterZ[ie];
                   temp=temp/(elementFaceNormX[jface + ie * 6] * ax + elementFaceNormY[jface + ie * 6] * ay + elementFaceNormZ[jface + ie * 6] * az);
                   //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                   valueindia += -temp;
               }
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

void  StaticFourier::_get_Re()
{
    for (int i = 0; i < numCell ; ++i) {
        Re[i]=0;
    }
    for (int i = 0; i < numCell ; ++i) {
        Re[i]-=elementHeatSource[i]*elementVolume[i];
    }
    for (int i = 0; i < numCell ; ++i) {
        for (int j = 0; j <elementFaceSize[i] ; ++j) {
            double ax,ay,az;
            if (elementFaceBound[j+i*6]>=0)
            {
                if (boundaryType[elementFaceBound[j+i*6]]==1)
                {
                    double temp = kappa[i] * elementFaceArea[j+i*6] ;
                    ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                    ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                    az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                    temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                    Re[i] += -temp*boundaryTemp[elementFaceBound[j+i*6]];
                }
                else if (boundaryType[elementFaceBound[j+i*6]]==5)
                {
                    //Re[i] += (temperature[elementFaceNeighobr[j+i*6]]-temperature[i])*boundaryResis[elementFaceBound[j+i*6]];
                }
                else if(boundaryType[elementFaceBound[j+i*6]]<0)
                {
                    if(boundaryThermal[elementFaceBound[j+i*6]]!=0||boundaryThermal[boundaryConnect[elementFaceBound[j+i*6]]]!=0)
                    {double temp = kappa[i] * elementFaceArea[j+i*6] ;
                    ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                    ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                    az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                    temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                    Re[i] += -temp*boundaryTemp[elementFaceBound[j+i*6]];}
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
                double gradientX=nodeX[elementFaceVertexes[1+j*4+i*6*4]]-nodeX[elementFaceVertexes[0+j*4+i*6*4]];
                double gradientY=nodeY[elementFaceVertexes[1+j*4+i*6*4]]-nodeY[elementFaceVertexes[0+j*4+i*6*4]];
                ax=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradientX*gradientX+gradientY*gradientY)*gradientX;
                ay=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradientX*gradientX+gradientY*gradientY)*gradientY;
            } else if (dimension == 3)
            {
                std::vector<double> d1(elementFaceVertexesSize[j+i*6], 0);
                std::vector<std::vector<double>> J1(elementFaceVertexesSize[j+i*6]);
                J1[0].push_back(elementFaceNormX[j+i*6]);
                J1[0].push_back(elementFaceNormY[j+i*6]);
                J1[0].push_back(elementFaceNormZ[j+i*6]);
                d1[0] = 0;

                for (int m = 0; m <  elementFaceVertexesSize[j+i*6]-1; m++)
                {
                    J1[m + 1].push_back(1.0/L_x * (nodeX[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeX[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1].push_back(1.0/L_x * (nodeY[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeY[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1].push_back(1.0/L_x * (nodeZ[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeZ[elementFaceVertexes[m+j*4+i*6*4]]));
                    d1[m + 1] = 1.0/L_x * (temperatureVertex[elementFaceVertexes[m+1+j*4+i*6*4]] - temperatureVertex[elementFaceVertexes[m+0+j*4+i*6*4]] );
                }
                vec a = Get_gradient_3D(J1, d1);
                ax=a.x;
                ay=a.y;
                az=a.z;
            }
            if(elementFaceBound[j+i*6]==-1)
            {
                Re[i] -= 2*(kappa[i]*kappa[elementFaceNeighobr[j+i*6]])/(kappa[i]+kappa[elementFaceNeighobr[j+i*6]]) * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
            }
            else {
                if (boundaryType[elementFaceBound[j + i * 6]] == 1||boundaryType[elementFaceBound[j + i * 6]] == 4) {
                    Re[i] -= kappa[i] * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
                } else if(boundaryType[elementFaceBound[j + i * 6]] <0)
                {
                    if(boundaryThermal[elementFaceBound[j+i*6]]!=0||boundaryThermal[boundaryConnect[elementFaceBound[j+i*6]]]!=0)
                    {
                        Re[i] -= kappa[i] * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
                    }
                }
            }

        }
    }
    for (int i = 0; i < numCell; ++i) {
        Re1(i)=Re[i];
    }
}

void StaticFourier::_get_heat_flux() const
{
    if(dimension==1)
    {
        for (int ie = 0; ie < numCell; ++ie)
        {
            if ( elementFaceBound[ie*6]==-1 && elementFaceBound[ie*6+1]==-1)
            {
                heatFluxXGlobal[ie]=-2*(kappa[ie]*kappa[ie+1])/(kappa[ie]+kappa[ie+1])* (temperature[ie+1] - temperature[ie]) / (elementVolume[ie+1] / 2 + elementVolume[ie] / 2);
            }
            else if (elementFaceBound[ie*6] != -1)
            {
                heatFluxXGlobal[ie]=-2.0/(kappa[ie]+kappa[ie+1])*(kappa[ie]*kappa[ie+1])* (temperature[ie+1] - temperature[ie]) / (elementVolume[ie+1] / 2 + elementVolume[ie] / 2);

            }
            else if (elementFaceBound[ie*6+1] != -1)
            {
                heatFluxXGlobal[ie]=-2.0/(kappa[ie]+kappa[ie-1])*(kappa[ie]*kappa[ie-1])* (temperature[ie] - temperature[ie-1]) / (elementVolume[ie] / 2 + elementVolume[ie-1] / 2);

            }
        }
    }
    else if(dimension==2)
    {
        for (int i = 0; i < numCell ; ++i)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }
            }
            if (nump < 2)
            {
                heatFluxXGlobal[i] = 0;
                heatFluxYGlobal[i] = 0;
            }
            else
            {

                vector<double> d1(elementNeighborList[i].size() , 0);
                //vector<vector<double>> J1(nump);
                for (int j = 0; j <elementNeighborList[i].size() ; ++j)
                {
                    d1[j] = 1.0/L_x * (temperature[elementNeighborList[i][j]] - temperature[i]);
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

                heatFluxXGlobal[i]=-kappa[i]*gradient.x;
                heatFluxYGlobal[i]=-kappa[i]*gradient.y;
                //cout<<endl;

            }
        }
    }
    else if(dimension==3)
    {
        for (int i = 0; i < numCell ; ++i) {

            vector<double> d1(elementNeighborList[i].size() , 0);
            //vector<vector<double>> J1(nump);
            for (int j = 0; j <elementNeighborList[i].size() ; ++j)
            {
                d1[j] = 1.0/L_x * (temperature[elementNeighborList[i][j]] - temperature[i]);
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

            heatFluxXGlobal[i]=-kappa[i]*gradient.x;
            heatFluxYGlobal[i]=-kappa[i]*gradient.y;
            heatFluxZGlobal[i]=-kappa[i]*gradient.z;
        }
    }
}

void StaticFourier::_get_bound_temp() const
{
    for (int ib = 0; ib < numBound; ++ib) {
        int i = boundaryCell[ib];
        int j = boundaryFace[ib];
        boundaryTemp[ib]=temperature[i]+(-heatFluxXGlobal[ib]/kappa[i]*(elementFaceCenterX[i*6+j]-elementCenterX[i])+
                -heatFluxYGlobal[ib]/kappa[i]*(elementFaceCenterY[i*6+j]-elementCenterY[i])+
                -heatFluxZGlobal[ib]/kappa[i]*(elementFaceCenterZ[i*6+j]-elementCenterZ[i]));
        if(boundaryType[ib]==1)
        {
            boundaryTemp[ib]=boundaryThermal[ib];
        } else if(boundaryType[ib]==4)
        {
        } else if(boundaryType[ib]<0&&boundaryThermal[ib]!=0)
        {
            boundaryTemp[ib]=boundaryThermal[ib]/2;
            boundaryTemp[boundaryConnect[ib]]=-boundaryThermal[ib]/2;
        }
    }
}

void StaticFourier::_get_bound_flux()
{
    for (int ib = 0; ib < numBound; ++ib) {
        int i = boundaryCell[ib];
        int ij=boundaryFace[ib];
        double gradientX=0;
        double gradientY=0;
        double gradientZ=0;
        if(dimension==1)
        {
            if (elementFaceBound[i*6] != -1)
            {
                gradientX= (heatFluxXGlobal[i+1] - heatFluxXGlobal[i]) / (elementVolume[i+1] / 2 + elementVolume[i] / 2);

            }
            else if (elementFaceBound[i*6+1] != -1)
            {
                gradientX=(heatFluxXGlobal[i] - heatFluxXGlobal[i-1]) / (elementVolume[i] / 2 + elementVolume[i-1] / 2);
            }
        }
        else if(dimension==2)
        {
            int nump = 0;
            for (int j = 0; j < elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    nump++;
                }
            }
            if (nump < 2)
            {
                gradientX = 0;
                gradientY = 0;
            }
            else
            {
                vector<double> d1(nump, 0);
                vector<vector<double>> J1(nump);
                nump = 0;
                for (int j = 0; j < elementFaceSize[i]; ++j)
                {
                    if (elementFaceBound[i*6+j]== -1)
                    {
                        d1[nump] =(heatFluxXGlobal[elementFaceNeighobr[j+i*6]] - heatFluxXGlobal[i])*elementFaceNormX[ij+i*6]
                                +(heatFluxYGlobal[elementFaceNeighobr[j+i*6]] - heatFluxYGlobal[i])*elementFaceNormY[ij+i*6];
                        nump++;
                    }

                }

                for (int m = 0; m < d1.size(); ++m)
                {
                    gradientX+= CellMatrix[i][0][m] * d1[m];
                    gradientY+= CellMatrix[i][1][m] * d1[m];
                }

            }
        }
        else if(dimension==3)
        {
            int nump = 0;
            for (int j = 0; j <elementFaceSize[i] ; ++j) {
                if (elementFaceNeighobr[j+i*6] >= 0)
                {
                    nump++;
                }

            }
            vector<double> d1(nump, 0);
            nump = 0;
            for (int j = 0; j <elementFaceSize[i]; ++j)
            {
                if (elementFaceBound[i*6+j] == -1)
                {
                    d1[nump] = 1.0/L_x * ((heatFluxXGlobal[elementFaceNeighobr[j+i*6]] - heatFluxXGlobal[i])*elementFaceNormX[ij+i*6]
                            +(heatFluxYGlobal[elementFaceNeighobr[j+i*6]] - heatFluxYGlobal[i])*elementFaceNormY[ij+i*6]
                            +(heatFluxZGlobal[elementFaceNeighobr[j+i*6]] - heatFluxZGlobal[i])*elementFaceNormZ[ij+i*6]);
                    nump++;
                }
            }


            for (int m = 0; m < d1.size(); ++m)
            {
                gradientX+= CellMatrix[i][0][m] * d1[m];
                gradientY+= CellMatrix[i][1][m] * d1[m];
                gradientZ+= CellMatrix[i][2][m] * d1[m];
            }
        }

        boundaryFlux[ib]= (heatFluxXGlobal[i]*elementFaceNormX[ij+i*6]
                +heatFluxYGlobal[i]*elementFaceNormY[ij+i*6]
                +heatFluxZGlobal[i]*elementFaceNormZ[ij+i*6])
                        +(gradientX*(elementFaceCenterX[i*6+ij]-elementCenterX[i])+
                        gradientY*(elementFaceCenterY[i*6+ij]-elementCenterY[i])+
                        gradientZ*(elementFaceCenterZ[i*6+ij]-elementCenterZ[i]));
        //cout<<boundaryFlux[ib]<<endl;
    }

}



void StaticFourier::_print_out() {
    std::ofstream output("Temperature.dat");
    for (int i = 0; i < numCell ; ++i) {
        output<<elementCenterX[i]<<" "<<elementCenterY[i]<<" "<<elementCenterZ[i]<<" "<<temperature[i]<<std::endl;
    }
    output.close();

    /*std::ofstream output1("TempcellVertex_macro.dat");
    for (int i = 0; i <numnode ; ++i) {
        output1<<nodeX[i]<<" "<<nodeY[i]<<" "<<nodeZ[i]<<" "<<temperatureVertex[i]<<std::endl;
    }
    output1.close();*/

    std::ofstream output2("HeatFlux.dat");
    for (int i = 0; i < numCell; ++i)
    {
        output2 << setw(15) << elementCenterX[i]  << setw(15) << elementCenterY[i] << setw(15) << elementCenterZ[i] << setw(15) << heatFluxXGlobal[i] << setw(15) << heatFluxYGlobal[i] << setw(15) << heatFluxZGlobal[i] << endl;
    }
    output2.close();

    ofstream outputVTK("Result_macro.vtk");
    outputVTK<<"# vtk DataFile Version 2.0"<<endl;
    outputVTK<<"RESULT"<<endl;
    outputVTK<<"ASCII"<<endl;
    outputVTK<<"DATASET UNSTRUCTURED_GRID"<<endl;
    outputVTK<<"POINTS "<<numnode<<" float"<<endl;
    for (int i = 0; i < numnode; ++i)
    {
        outputVTK<<nodeX[i]/L_x <<" "
        << nodeY[i]/L_y <<" " << nodeZ[i]/L_z  << endl;
    }
    int totalnum=0;
    for (int i = 0; i < numCell; ++i) {
        totalnum+=elementVertexesSize[i]+1;
    }
    outputVTK<<endl;
    outputVTK<<"CELLS "<<numCell<<" "<<totalnum<<endl;
    for (int i = 0; i < numCell; ++i) {
        outputVTK<<elementVertexesSize[i];
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            outputVTK<<" "<<elementVertexes[j+i*8];
        }
        outputVTK<<endl;
    }
    outputVTK<<endl;
    outputVTK<<"CELL_TYPES "<<numCell<<endl;
    for (int i = 0; i < numCell; ++i) {
        if( dimension==2&&elementVertexesSize[i]==4)
        {
            outputVTK<<8<<endl;
        } else if( dimension==2&&elementVertexesSize[i]==3)
        {
            outputVTK<<5<<endl;
        }else if( dimension==3&&elementVertexesSize[i]==4)
        {
            outputVTK<<10<<endl;
        }
        else if( dimension==3&&elementVertexesSize[i]==8)
        {
            outputVTK<<11<<endl;
        }
        else if( dimension==3&&elementVertexesSize[i]==6)
        {
            outputVTK<<13<<endl;
        }
    }
    outputVTK<<endl;
    outputVTK<<"CELL_DATA"<<" "<<numCell<<endl;
    outputVTK<<"SCALARS"<<" Temperature"<<" float 1"<<endl;
    outputVTK<<"LOOKUP_TABLE default"<<endl;
    for (int i = 0; i < numCell; ++i) {
        outputVTK<<temperature[i]<<endl;
    }
    outputVTK<<"VECTORS"<<" HeatFLux"<<" float"<<endl;
    //outputVTK<<"LOOKUP_TABLE default"<<endl;
    for (int i = 0; i < numCell; ++i) {
        outputVTK<<heatFluxXGlobal[i]<<" "<<heatFluxYGlobal[i]<<" "<<heatFluxZGlobal[i]<<endl;
    }
}

void  StaticFourier::_solve() {
    //_set_initial();
    _get_coefficient_macro();
    _get_CellMatrix();
    double error = 1000;

    while (error > 1e-5)
    {
        for (int i = 0; i < numCell ; ++i) {
            temperatureOld[i]=temperature[i];
        }
        _set_vertex_temperature();
        if(worldRank==0)
            _print_out();
        _get_Re();
        for (int i = 0; i < numCell; ++i) {

            //cout<<Re[i]<<endl;
        }
        //double *sol = _solve_matrix((int*)csrRowPtr, (int*)csrColInd, csrVal, Re, numCell);
        sol1=cg_macro.solve(Re1);
        for (int i = 0; i < numCell ; ++i) {
            temperature[i]=sol1(i);
            //cout<< temperature[i]<<endl;
        }
        _get_heat_flux();
        //_get_bound_temp();
        //_get_bound_flux();

        error=0;
        for (int i = 0; i < numCell; ++i) {
            error+= pow(temperatureOld[i]-temperature[i],2);
        }
        error= error / numCell;
        if(worldRank==0)
        {_print_out();}

        if(worldRank==0)
        std::cout<<"iterations error  "<<error<<std::endl;
    }


}