//
// Created by huyue on 2023/1/17.
//
//
// Created by yuehu on 2022/12/12.
//

#include "TransientFourier/TransientFourier.h"
#include <algorithm>
#include <iomanip>
#include <petscksp.h>
#include "mpi.h"
#include "utility/utility.h"

TransientFourier::TransientFourier(BTEMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc,
                                   int world_rank) {
    deltaT=10;
    totaltime=500;
    Dimension=mesh->Dimension;
    numcell=mesh->Elements.size();
    numnode=mesh->Nodes.size();
    temperaturen=new double [numcell];
    temperaturen_1=new double [numcell];
    temperatureTemp=new double [numcell];
    temperatureVertex=new double [mesh->Nodes.size()];
    kappa=new double [numcell];
    capacity=new double [numcell];
    for (int i = 0; i <numcell ; ++i) {
        kappa[i]=bands->kappabulk[mesh->Elements[i].matter];
        capacity[i]=bands->capacitybulk[mesh->Elements[i].matter];
    }

    elementFaceSize=new int [numcell];
    for (int i = 0; i < numcell ; ++i) {
        elementFaceSize[i]=mesh->Elements[i].faces.size();
    }
    elementFaceNeighobr=new int [numcell*6];
    elementFaceArea=new double [numcell*6];
    for (int i = 0; i < numcell; ++i) {
        for (int j = 0; j < elementFaceSize[i] ; ++j) {
            elementFaceArea[j+i*6]=mesh->Elements[i].faces[j].area;
            elementFaceNeighobr[j+i*6]=mesh->Elements[i].faces[j].neighbor;
        }
    }

    elementFaceNormX=new double [numcell*6];
    elementFaceNormY=new double [numcell*6];
    elementFaceNormZ=new double [numcell*6];

    elementFaceCenterX=new double [numcell*6];
    elementFaceCenterY=new double [numcell*6];
    elementFaceCenterZ=new double [numcell*6];
    elementFaceBound=new int [numcell*6];
    for (int i = 0; i < numcell*6 ; ++i) {
        elementFaceNormX[i]=0;
        elementFaceCenterX[i]=0;
        elementFaceNormY[i]=0;
        elementFaceCenterY[i]=0;
        elementFaceNormZ[i]=0;
        elementFaceCenterZ[i]=0;
    }
    for (int i = 0; i < numcell; ++i) {
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


    elementFaceVertexesSize=new int [numcell*6];
    for (int i = 0; i <numcell ; ++i) {
        for (int j = 0; j < elementFaceSize[i]  ; ++j) {
            elementFaceVertexesSize[j+i*6]=mesh->Elements[i].faces[j].vertexes.size();
        }
    }

    elementFaceVertexes=new int [numcell*6*4];
    for (int i = 0; i <numcell ; ++i) {
        for (int j = 0; j < elementFaceSize[i]  ; ++j) {
            for (int k = 0; k < elementFaceVertexesSize[j+i*6]; ++k) {
                elementFaceVertexes[k+j*4+i*6*4]=mesh->Elements[i].faces[j].vertexes[k];
            }
        }
    }

    elementVertexesSize=new int [numcell];
    for (int i = 0; i < numcell; ++i) {
        elementVertexesSize[i]=mesh->Elements[i].vertexes.size();
    }
    elementVertexes=new int [numcell*8];
    for (int i = 0; i < numcell; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            elementVertexes[j + i * 8] = mesh->Elements[i].vertexes[j];
            //std::cout<<elementVertexes[j + i * 8]<<std::endl;
        }
    }

    elementVolume=new double [numcell];
    for (int i = 0; i < numcell ; ++i) {
        elementVolume[i]=mesh->Elements[i].volume;
    }
    elementHeatSource=new double [numcell];
    for (int i = 0; i < numcell ; ++i) {
        elementHeatSource[i]=mesh->Elements[i].heat_source;
    }


    elementCenterX=new double [numcell];
    elementCenterY=new double [numcell];
    elementCenterZ=new double [numcell];
    for (int i = 0; i < numcell ; ++i) {
        elementCenterX[i]=mesh->Elements[i].center.x;
        elementCenterY[i]=mesh->Elements[i].center.y;
        elementCenterZ[i]=mesh->Elements[i].center.z;
    }
    numBound=mesh->Boundaries.size();
    boundaryType=new int [numBound];
    boundaryTemp=new double [numBound];
    for (int i = 0; i <numBound ; ++i) {
        boundaryType[i]=mesh->Boundaries[i].type;
        boundaryTemp[i]=mesh->Boundaries[i].Temperature;
    }

    nodeX=new double [numnode];
    nodeY=new double [numnode];
    nodeZ=new double [numnode];
    for (int i = 0; i < numnode; ++i) {
        nodeX[i]=mesh->Nodes[i].x;
        nodeY[i]=mesh->Nodes[i].y;
        nodeZ[i]=mesh->Nodes[i].z;
    }

    csrRowPtr= new unsigned int[numcell + 1];
    csrColInd = new unsigned int[7 * numcell + 1];
    csrVal = new double[7 * numcell + 1];
    Re=new double [numcell];

}

void  TransientFourier::_set_initial() {
    for (int i = 0; i <numcell ; ++i) {
        temperaturen[i]=0;
        temperaturen_1[i]=0;
        temperatureTemp[i]=0;
        Re[i]=0;
    }
}


void TransientFourier::_set_vertex_temperature() {


    double *temp=new double [numnode];
    double *R_n=new double [numnode];

    for (int i = 0; i <numnode ; ++i) {
        temp[i]=0;
        R_n[i]=0;
    }

    for (int i = 0; i <numcell ; ++i) {
        for (int j = 0; j < elementVertexesSize[i]; ++j) {
            double distance= sqrt(pow(nodeX[elementVertexes[j+i*8]]-elementCenterX[i],2)+
                                  pow(nodeY[elementVertexes[j+i*8]]-elementCenterY[i],2)+
                                  pow(nodeZ[elementVertexes[j+i*8]]-elementCenterZ[i],2));
            temp[elementVertexes[j+i*8]]+=(temperaturen[i]+temperaturen_1[i])/2/
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

void  TransientFourier::_get_coefficient()
{
    int csrRowPtr_iter = 0, csrColInd_iter = 0, csrVal_iter = 0;
    csrRowPtr[0] = 0;
    csrRowPtr_iter++;
    for (int i = 0; i < numcell ; ++i) {
        std::vector<std::pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(6 + 1);
        double valueindia = 0;

        for (int j = 0; j <elementFaceSize[i] ; ++j) {
            double temp = kappa[i] * elementFaceArea[j+i*6]/2 ;

            if (elementFaceNeighobr[j+i*6] >= 0)
            {
                double ax=elementCenterX[elementFaceNeighobr[j+i*6]]-elementCenterX[i];
                double ay=elementCenterY[elementFaceNeighobr[j+i*6]]-elementCenterY[i];
                double az=elementCenterZ[elementFaceNeighobr[j+i*6]]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6] * az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
                compressed_Ke.emplace_back(elementFaceNeighobr[j+i*6], temp);

            }
            else if (boundaryType[elementFaceBound[j+i*6]]==1)
            {
                double ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                double ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                double az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                //temp=temp/ sqrt(ax*ax+ay*ay+az*az);
                valueindia += -temp;
            }
        }
        compressed_Ke.emplace_back(i,valueindia-1/deltaT*elementVolume[i]*capacity[i]);
        sort(compressed_Ke.begin(), compressed_Ke.end());
        for (int num=0;num<compressed_Ke.size();num++) {
            csrColInd[csrColInd_iter++] = compressed_Ke[num].first;
            //cout<< compressed_Ke[num].first<<" ";
            csrVal[csrVal_iter++] = compressed_Ke[num].second;
            //cout<<compressed_Ke[num].second<<endl;
        }

        csrRowPtr[csrRowPtr_iter] = csrRowPtr[csrRowPtr_iter - 1] + compressed_Ke.size();

        csrRowPtr_iter++;
    }
}

void  TransientFourier::_get_Re()
{
    for (int i = 0; i <numcell ; ++i) {
        Re[i]=0;
    }
    for (int i = 0; i <numcell ; ++i) {
        Re[i]-=elementHeatSource[i]*elementVolume[i];
        Re[i]-=temperaturen_1[i]/deltaT*elementVolume[i]*capacity[i];
    }
    for (int i = 0; i < numcell ; ++i) {
        for (int j = 0; j <elementFaceSize[i] ; ++j) {
            double ax,ay,az;
            if (elementFaceBound[j+i*6]>=0)
            {
                if (boundaryType[elementFaceBound[j+i*6]]==1)
                {
                    double temp = kappa[i] * elementFaceArea[j+i*6]/2;
                    ax=elementFaceCenterX[j+i*6]-elementCenterX[i];
                    ay=elementFaceCenterY[j+i*6]-elementCenterY[i];
                    az=elementFaceCenterZ[j+i*6]-elementCenterZ[i];
                    temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                    Re[i] += -temp*boundaryTemp[elementFaceBound[j+i*6]];

                    Re[i] += -temp*boundaryTemp[elementFaceBound[j+i*6]];

                    Re[i] += +temp*temperaturen_1[i];

                }
            } else
            {
                double temp = kappa[i] * elementFaceArea[j+i*6]/2;
                ax=elementCenterX[elementFaceNeighobr[j+i*6]]-elementCenterX[i];
                ay=elementCenterY[elementFaceNeighobr[j+i*6]]-elementCenterY[i];
                az=elementCenterZ[elementFaceNeighobr[j+i*6]]-elementCenterZ[i];
                temp=temp/(elementFaceNormX[j+i*6] * ax+elementFaceNormY[j+i*6]* ay+elementFaceNormZ[j+i*6]* az);
                Re[i] += -temp*temperaturen_1[elementFaceNeighobr[j+i*6]];
                Re[i] += +temp*temperaturen_1[i];

            }
            double a2x=elementFaceNormX[j+i*6]-ax/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            double a2y=elementFaceNormY[j+i*6]-ay/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            double a2z=elementFaceNormZ[j+i*6]-az/(ax*elementFaceNormX[j+i*6]+ay*elementFaceNormY[j+i*6]+az*elementFaceNormZ[j+i*6]);
            if (Dimension==2)
            {
                double gradientX=nodeX[elementFaceVertexes[1+j*4+i*6*4]]-nodeX[elementFaceVertexes[0+j*4+i*6*4]];
                double gradientY=nodeY[elementFaceVertexes[1+j*4+i*6*4]]-nodeY[elementFaceVertexes[0+j*4+i*6*4]];
                ax=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradientX*gradientX+gradientY*gradientY)*gradientX;
                ay=(temperatureVertex[elementFaceVertexes[1+j*4+i*6*4]]-temperatureVertex[elementFaceVertexes[0+j*4+i*6*4]])/(gradientX*gradientX+gradientY*gradientY)*gradientY;
            } else if (Dimension==3)
            {
                std::vector<double> d1(elementFaceVertexesSize[j+i*6], 0);
                std::vector<std::vector<double>> J1(elementFaceVertexesSize[j+i*6]);
                J1[0].push_back(elementFaceNormX[j+i*6]);
                J1[0].push_back(elementFaceNormY[j+i*6]);
                J1[0].push_back(elementFaceNormZ[j+i*6]);
                d1[0] = 0;

                for (int m = 0; m <  elementFaceVertexesSize[j+i*6]-1; m++)
                {
                    J1[m + 1].push_back(1e9 * (nodeX[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeX[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1].push_back(1e9 * (nodeY[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeY[elementFaceVertexes[m+j*4+i*6*4]]));
                    J1[m + 1].push_back(1e9 * (nodeZ[elementFaceVertexes[m+1+j*4+i*6*4]] - nodeZ[elementFaceVertexes[m+j*4+i*6*4]]));
                    d1[m + 1] = 1e9 * (temperatureVertex[elementFaceVertexes[m+1+j*4+i*6*4]] - temperatureVertex[elementFaceVertexes[m+0+j*4+i*6*4]] );
                }
                vec a = Get_gradient_3D(J1, d1);
                ax=a.x;
                ay=a.y;
                az=a.z;
            }
            if (elementFaceBound[j+i*6]>=0) {
                if (boundaryType[elementFaceBound[j + i * 6]] == 1) {
                    Re[i] -= kappa[i] * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
                }
            } else
            {
                Re[i] -= kappa[i] * elementFaceArea[j + i * 6] * (a2x * ax + a2y * ay + a2z * az);
            }

        }
    }

}

void TransientFourier::_print_out() {
    std::ofstream output("Tempcell_old.dat");
    for (int i = 0; i <numcell ; ++i) {
        output<<elementCenterX[i]<<" "<<elementCenterY[i]<<" "<<elementCenterZ[i]<<" "<<temperaturen[i]<<std::endl;
    }
    output.close();
}

void  TransientFourier::_solve()
{
    _set_initial();
    _get_coefficient();

    double difference=1000;
    for (int it = 0; it < totaltime/deltaT ; ++it) {
        for (int i = 0; i <numcell ; ++i) {
            temperaturen_1[i]=temperaturen[i];
        }
        double error = 1000;
        while (error > 1e-6)
        {
            for (int i = 0; i <numcell ; ++i) {
                temperatureTemp[i]=temperaturen[i];
            }
            _set_vertex_temperature();
            _print_out();
            _get_Re();

            double *sol = _solve_matrix((int*)csrRowPtr, (int*)csrColInd, csrVal, Re,numcell);
            for (int i = 0; i <numcell ; ++i) {
                temperaturen[i]=sol[i];
            }
            error=0;
            for (int i = 0; i < numcell; ++i) {
                error+= pow(temperaturen[i]-temperatureTemp[i],2);
            }
            error=error/numcell;

            std::cout<<"internalerror  "<<error<<std::endl;
        }
        difference=0;
        for (int i = 0; i < numcell; ++i) {
            difference+= pow(temperaturen[i]-temperaturen_1[i],2);
        }
        std::cout<<"#####################################################"<<std::endl;
        std::cout<<"currenttime "<<(it+1)*deltaT<<", difference  "<<difference<<std::endl;
        _print_out();

    }


}