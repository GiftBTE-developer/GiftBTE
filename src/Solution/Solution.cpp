//
// Created by yuehu on 2023/3/6.
//
#include <iomanip>
#include "Solution/Solution.h"

SolutionAll::~SolutionAll() {
   delete [] solution.Temperature;
   delete [] solution.heatFluxX;
   delete [] solution.heatFluxY;
   delete [] solution.heatFluxZ;

   for (int j = 0; j < solution.numBandLocal; ++j) {

       for (int k = 0; k < solution.numDirectionLocal; ++k) {
          delete [] solution.energyDensity[j][k];
       }
       delete [] solution.energyDensity[j];
   }
   delete [] solution.energyDensity;
   delete [] solution.ebound;
   delete [] solution.boundaryFlux;
   delete [] solution.boundaryTemp;
   delete [] solution.Temperature_old;
   delete []  solution.heatFluxX_old;
   delete [] solution.heatFluxY_old;
   delete []  solution.heatFluxZ_old;

}

SolutionAll::SolutionAll(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles,int num_proc, int world_rank) {
    solution.solverType="BTE";
    solution.numCell=mesh->volumeElementsFourier.size();
    solution.Temperature=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxX=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxY=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxZ=new double [mesh->volumeElementsFourier.size()];
    int numDirectionLocal;
    int numBandLocal;
    int numDirection=angles->direction.size();
    int numBand=bands->nband;
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
    solution.numDirectionLocal=numDirectionLocal;
    solution.numBandLocal=numBandLocal;
    solution.numDirection=numDirection;
    solution.numBand=numBand;
    solution.energyDensity= new double ** [numBandLocal];
    for (int j = 0; j < numBandLocal; ++j) {
        solution.energyDensity[j]=new double * [numDirectionLocal];
        for (int k = 0; k < numDirectionLocal; ++k) {
            solution.energyDensity[j][k]=new double [mesh->volumeElementsFourier.size()];
        }
    }
    solution.numBound=mesh->FourierMeshes.Boundaries.size();
    solution.ebound = new double[numBand * numDirection * solution.numBound * 2];
    solution.boundaryFlux = new double [solution.numBound];
    solution.boundaryTemp = new double [solution.numBound];
    solution.Temperature_old=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxX_old=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxY_old=new double [mesh->volumeElementsFourier.size()];
    solution.heatFluxZ_old=new double [mesh->volumeElementsFourier.size()];

}

void SolutionAll::_set_initial(DistributeMesh *mesh,BTEBand *bands, BTEAngle *angles) {

    for (int j = 0; j < solution.numCell; ++j) {
        //solution.Temperature[j]=0;
        solution.Temperature[j]=mesh->FourierMeshes.Elements[j].initial_temperature;//jiaxuan: for both steady and transient state
                solution.heatFluxX[j]=0;
        solution.heatFluxY[j]=0;
        solution.heatFluxZ[j]=0;
    }
    int *boundaryType=new int [solution.numBound];
    auto *boundaryTemp=new double [solution.numBound];
    int *boundaryCell=new int [solution.numBound];
    int *boundaryFace=new int [solution.numBound];
    int *boundaryConnect=new int [solution.numBound];

    for (int j = 0; j <solution.numBound ; ++j) {
        boundaryType[j]=mesh->FourierMeshes.Boundaries[j].type;
        boundaryTemp[j]=mesh->FourierMeshes.Boundaries[j].Temperature;
        boundaryCell[j]=mesh->FourierMeshes.Boundaries[j].cellindex[0];
        boundaryFace[j]=mesh->FourierMeshes.Boundaries[j].faceindex[0];
        boundaryConnect[j]=mesh->FourierMeshes.Boundaries[j].connection;
    }
    for (int ib = 0; ib < solution.numBound; ++ib)
    {
        if (boundaryType[ib] == 1)
        {
            solution.boundaryTemp[ib]=boundaryTemp[ib];
        } else if(boundaryType[ib]==2||boundaryType[ib]==3)
        {
            solution.boundaryFlux[ib]=0;
        } else if(boundaryType[ib]==0)
        {
            solution.boundaryFlux[ib]=0;
        }
        else if(boundaryType[ib]<0&&boundaryTemp[ib]!=0)
        {
            solution.boundaryTemp[ib]=boundaryTemp[ib]/2;
            solution.boundaryTemp[boundaryConnect[ib]]=-boundaryTemp[ib]/2;
        }
    }
    delete [] boundaryType;
    delete [] boundaryTemp;
    delete [] boundaryCell;
    delete [] boundaryFace;
    delete [] boundaryConnect;

}



void SolutionAll::_print_out(DistributeMesh *mesh)
{
    ofstream output("Temperature.dat");
    for (int i = 0; i < mesh->volumeElementsFourier.size(); ++i) {
        output<<mesh->FourierMeshes.Elements[i].center.x<<" "<<mesh->FourierMeshes.Elements[i].center.y
        <<" "<<mesh->FourierMeshes.Elements[i].center.z<<" "<<solution.Temperature[i]<<endl;
    }
    output.close();

    ofstream outputheat("HeatFlux.dat");

    for (int i = 0; i < mesh->volumeElementsFourier.size(); ++i)
    {
        outputheat<<setw(15) <<mesh->FourierMeshes.Elements[i].center.x << setw(15)
        << mesh->FourierMeshes.Elements[i].center.y <<setw(15) << mesh->FourierMeshes.Elements[i].center.z <<setw(15)
        << solution.heatFluxX[i] <<setw(15) << solution.heatFluxY[i]<<setw(15) <<solution.heatFluxZ[i] << endl;
    }

    outputheat.close();

    double sumHx=0;
    double sumHy=0;
    double sumHz=0;
    double sumA=0;
    for (int i = 0; i < mesh->volumeElementsFourier.size(); ++i) {
        sumHx+=solution.heatFluxX[i]*mesh->FourierMeshes.Elements[i].volume;
        sumHy+=solution.heatFluxY[i]*mesh->FourierMeshes.Elements[i].volume;
        sumHz+=solution.heatFluxZ[i]*mesh->FourierMeshes.Elements[i].volume;
        sumA+=mesh->FourierMeshes.Elements[i].volume;
    }
    //cout<<"======================================heat flux ave====================="<<endl;
    //ofstream outave("aveHeat");
    //outave<<sumHx/sumA<<"   "<<sumHy/sumA<<"   "<<sumHz/sumA<<endl;
    //cout<<"========================================================================"<<endl;

    ofstream outputVTK("Result.vtk");
    outputVTK<<"# vtk DataFile Version 2.0"<<endl;
    outputVTK<<"RESULT"<<endl;
    outputVTK<<"ASCII"<<endl;
    outputVTK<<"DATASET UNSTRUCTURED_GRID"<<endl;
    outputVTK<<"POINTS "<<mesh->FourierMeshes.Nodes.size()<<" float"<<endl;
    for (int i = 0; i < mesh->FourierMeshes.Nodes.size(); ++i)
    {
        outputVTK<<mesh->FourierMeshes.Nodes[i].x/mesh->FourierMeshes.L_x <<" "
        << mesh->FourierMeshes.Nodes[i].y/mesh->FourierMeshes.L_y <<" " << mesh->FourierMeshes.Nodes[i].z/mesh->FourierMeshes.L_z  << endl;
    }
    int totalnum=0;
    for (int i = 0; i < mesh->FourierMeshes.Elements.size(); ++i) {
        totalnum+=mesh->FourierMeshes.Elements[i].vertexes.size()+1;
    }
    outputVTK<<endl;
    outputVTK<<"CELLS "<<mesh->FourierMeshes.Elements.size()<<" "<<totalnum<<endl;
    for (int i = 0; i < mesh->FourierMeshes.Elements.size(); ++i) {
        outputVTK<<mesh->FourierMeshes.Elements[i].vertexes.size();
        for (int j = 0; j < mesh->FourierMeshes.Elements[i].vertexes.size(); ++j) {
            outputVTK<<" "<<mesh->FourierMeshes.Elements[i].vertexes[j];
        }
        outputVTK<<endl;
    }
    outputVTK<<endl;
    outputVTK<<"CELL_TYPES "<<mesh->FourierMeshes.Elements.size()<<endl;
    for (int i = 0; i < mesh->FourierMeshes.Elements.size(); ++i) {
        if( mesh->FourierMeshes.Dimension==2&&mesh->FourierMeshes.Elements[i].vertexes.size()==4)
        {
            outputVTK<<8<<endl;
        } else if(mesh->FourierMeshes.Dimension==2&&mesh->FourierMeshes.Elements[i].vertexes.size()==3)
        {
            outputVTK<<5<<endl;
        }else if(mesh->FourierMeshes.Dimension==3&&mesh->FourierMeshes.Elements[i].vertexes.size()==4)
        {
            outputVTK<<10<<endl;
        }
        else if(mesh->FourierMeshes.Dimension==3&&mesh->FourierMeshes.Elements[i].vertexes.size()==8)
        {
            outputVTK<<11<<endl;
        }
        else if(mesh->FourierMeshes.Dimension==3&&mesh->FourierMeshes.Elements[i].vertexes.size()==6)
        {
            outputVTK<<13<<endl;
        }
    }
    outputVTK<<endl;
    outputVTK<<"CELL_DATA"<<" "<<mesh->FourierMeshes.Elements.size()<<endl;
    outputVTK<<"SCALARS"<<" Temperature"<<" float 1"<<endl;
    outputVTK<<"LOOKUP_TABLE default"<<endl;
    for (int i = 0; i < mesh->FourierMeshes.Elements.size(); ++i) {
        outputVTK<<solution.Temperature[i]<<endl;
    }
    outputVTK<<"VECTORS"<<" HeatFLux"<<" float"<<endl;
    //outputVTK<<"LOOKUP_TABLE default"<<endl;
    for (int i = 0; i < mesh->FourierMeshes.Elements.size(); ++i) {
        outputVTK<<solution.heatFluxX[i]<<" "<<solution.heatFluxY[i]<<" "<<solution.heatFluxZ[i]<<endl;
    }
    /*ofstream outputboundtemp("Boundary_Temp_Flux.dat");

    for (int ib = 0; ib < mesh->FourierMeshes.Boundaries.size(); ++ib) {
        int i=mesh->FourierMeshes.Boundaries[ib].cellindex[0];
        int j=mesh->FourierMeshes.Boundaries[ib].faceindex[0];
        outputboundtemp<<setw(15)<<mesh->FourierMeshes.Boundaries[ib].index<<setw(15)
        <<mesh->FourierMeshes.Elements[i].faces[j].center.x<<setw(15)<<mesh->FourierMeshes.Elements[i].faces[j].center.y
        <<setw(15)<<mesh->FourierMeshes.Elements[i].faces[j].center.z<<setw(15)<<solution.boundaryTemp[ib]
        <<setw(15)<<solution.boundaryFlux[ib]<<endl;
    }

    outputboundtemp.close();*/
}


void SolutionAll::_Fourier_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, int num_proc, int world_rank) const
{
    StaticFourier solver(&mesh->FourierMeshes, bcs, bands,num_proc, world_rank);
    for (int j = 0; j < solution.numCell; ++j) {
        solver.temperature[j]=solution.Temperature[j];
        solver.heatFluxXGlobal[j]=solution.heatFluxX[j];
        solver.heatFluxYGlobal[j]=solution.heatFluxY[j];
        solver.heatFluxZGlobal[j]=solution.heatFluxZ[j];
    }
    for (int j = 0; j < solution.numBound; ++j) {
        solver.boundaryTemp[j]=solution.boundaryTemp[j];
        solver.boundaryFlux[j]=solution.boundaryFlux[j];
    }
    solver._solve();
    for (int j = 0; j < solution.numCell; ++j) {
        solution.Temperature[j] = solver.temperature[j];
        solution.heatFluxX[j] = solver.heatFluxXGlobal[j];
        solution.heatFluxY[j] = solver.heatFluxYGlobal[j];
        solution.heatFluxZ[j] = solver.heatFluxZGlobal[j];
    }
    for (int j = 0; j < solution.numBound; ++j) {
        solution.boundaryTemp[j]=solver.boundaryTemp[j];
        solution.boundaryFlux[j]=solver.boundaryFlux[j];
    }
}

void SolutionAll::_BTE_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank
                                  ,int Use_Backup, int Num_Max_Iter,string Order, string Method, string Matrix_solver, double error_temp_limit, double error_flux_limit) const
{
    MPI_Barrier(MPI_COMM_WORLD);
    StaticBTESynthetic solver(&mesh->BTEMeshes, bcs, bands,angles,num_proc, world_rank);

    for (int i = 0; i < mesh->volumeElementsBTE.size(); ++i) {
        solver.temperature[i]=solution.Temperature[mesh->elementProjection[i]];
        solver.heatFluxXGlobal[i]=solution.heatFluxX[mesh->elementProjection[i]];
        solver.heatFluxYGlobal[i]=solution.heatFluxY[mesh->elementProjection[i]];
        solver.heatFluxZGlobal[i]=solution.heatFluxZ[mesh->elementProjection[i]];
    }

    for (int ib = 0; ib < mesh->BTEMeshes.Boundaries.size(); ++ib) {
        int ie = mesh->BTEMeshes.Boundaries[ib].cellindex[0];
        int jface = mesh->BTEMeshes.Boundaries[ib].faceindex[0];
        if(mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].faces[jface].bound==-1)
        {
            double distance1= get_distance(mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].center,mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].faces[jface].center);
            double distance2= get_distance(mesh->FourierMeshes.Elements[mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].faces[jface].neighbor].center,mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].faces[jface].center);
            solver.boundaryThermal[ib]=solution.Temperature[mesh->elementProjection[ie]]/distance1+ solution.Temperature[mesh->FourierMeshes.Elements[mesh->elementProjection[ie]].faces[jface].neighbor]/distance2;
            solver.boundaryThermal[ib]=solver.boundaryThermal[ib]/(1.0/distance1+1.0/distance2);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(Order=="1")
    {
        if(Method=="Synthetic")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve_firstorder(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_firstorder_Iterative(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else if(Method=="DOM")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve_DOM_firstorder(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);
            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_DOM_firstorder_Iterative(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);
            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else
        {
            cout<<"no this method"<<endl;
            exit(1);
        }
    } else if (Order=="2")
    {
        if(Method=="Synthetic")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_Iterative(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else if(Method=="DOM")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve_DOM(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_DOM_Iterative(Use_Backup,Num_Max_Iter,0,error_temp_limit,error_flux_limit);

            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else
        {
            cout<<"no this method"<<endl;
            exit(1);
        }
    } else if(Order=="2limit")
    {
        if(Method=="Synthetic")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve(Use_Backup,Num_Max_Iter,1,error_temp_limit,error_flux_limit);

            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_Iterative(Use_Backup,Num_Max_Iter,1,error_temp_limit,error_flux_limit);

            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else if(Method=="DOM")
        {
            if(Matrix_solver=="LU")
            {
                solver.solve_DOM(Use_Backup,Num_Max_Iter,1,error_temp_limit,error_flux_limit);

            }
            else if(Matrix_solver=="Iterative")
            {
                solver.solve_DOM_Iterative(Use_Backup,Num_Max_Iter,1,error_temp_limit,error_flux_limit);

            } else
            {
                cout<<"no this matrix solver"<<endl;
                exit(1);
            }
        } else
        {
            cout<<"no this method"<<endl;
            exit(1);
        }
    } else
    {
        cout<<"no this matrix solver"<<endl;
        exit(1);
    }

    for (int i = 0; i < mesh->volumeElementsBTE.size(); ++i) {
        solution.Temperature[mesh->elementProjection[i]]=solver.totalEnergy[i];
        solution.heatFluxX[mesh->elementProjection[i]]=solver.heatFluxXGlobal[i];
        solution.heatFluxY[mesh->elementProjection[i]]=solver.heatFluxYGlobal[i];
        solution.heatFluxZ[mesh->elementProjection[i]]=solver.heatFluxZGlobal[i];
    }
  }


void SolutionAll::_get_solver(int i, DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank
                              ,int Use_Backup, int Num_Max_Iter, int Use_Sythetic, int Use_Limiter, double error_temp_limit, double error_flux_limit ){

    if(solutions[i].solverType=="BTE")
    {
        StaticBTESynthetic solver(&mesh->mesh[i], bcs, bands, angles, num_proc, world_rank);
        for (int j = 0; j < solutions[i].numBandLocal; ++j) {
            for (int k = 0; k < solutions[i].numDirectionLocal; ++k) {
                for (int l = 0; l < solutions[i].numCell; ++l) {
                    solver.energyDensity[j][k][l]=solutions[i].energyDensity[j][k][l];
                }
            }
        }
        for (int j = 0; j < solutions[i].numCell; ++j) {
            solver.temperature[j]=solutions[i].Temperature[j];
            solver.heatFluxXGlobal[j]=solutions[i].heatFluxX[j];
            solver.heatFluxYGlobal[j]=solutions[i].heatFluxY[j];
            solver.heatFluxZGlobal[j]=solutions[i].heatFluxZ[j];
        }
        for (int j = 0; j < solutions[i].numBand * solutions[i].numDirection * solutions[i].numBound * 2; ++j) {
            solver.ebound[j] = solutions[i].ebound[j];
        }
        for (int j = 0; j < solutions[i].numBandLocal; ++j) {
            for (int k = 0; k < solutions[i].numDirectionLocal; ++k) {
                for (int l = 0; l < solutions[i].numCell; ++l) {
                    solutions[i].energyDensity[j][k][l]=solver.energyDensity[j][k][l];
                }
            }
        }
        for (int j = 0; j < solutions[i].numCell; ++j) {
            solutions[i].Temperature[j]=solver.temperature[j];
            solutions[i].heatFluxX[j]=solver.heatFluxXGlobal[j];
            solutions[i].heatFluxY[j]=solver.heatFluxYGlobal[j];
            solutions[i].heatFluxZ[j]=solver.heatFluxZGlobal[j];
        }
        for (int j = 0; j < solutions[i].numBand * solutions[i].numDirection * solutions[i].numBound * 2; ++j) {
            solutions[i].ebound[j] = solver.ebound[j];
        }
        for (int j = 0; j < solutions[i].numBound; ++j) {
            solutions[i].boundaryTemp[j]= solver.boundaryTemp[j];
            solutions[i].boundaryFlux[j]=solver.boundaryFlux[j];
        }

    }
    else if (solutions[i].solverType=="Fourier")
    {
        StaticFourier solver(&mesh->mesh[i], bcs, bands,num_proc, world_rank);
        for (int j = 0; j < solutions[i].numCell; ++j) {
            solver.temperature[j]=solutions[i].Temperature[j];
            solver.heatFluxXGlobal[j]=solutions[i].heatFluxX[j];
            solver.heatFluxYGlobal[j]=solutions[i].heatFluxY[j];
            solver.heatFluxZGlobal[j]=solutions[i].heatFluxZ[j];
        }
        for (int j = 0; j < solutions[i].numBound; ++j) {
            solver.boundaryTemp[j]=solutions[i].boundaryTemp[j];
            solver.boundaryFlux[j]=solutions[i].boundaryFlux[j];
        }
        solver._solve();
        for (int j = 0; j < solutions[i].numCell; ++j) {
            solutions[i].Temperature[j] = solver.temperature[j];
            solutions[i].heatFluxX[j] = solver.heatFluxXGlobal[j];
            solutions[i].heatFluxY[j] = solver.heatFluxYGlobal[j];
            solutions[i].heatFluxZ[j] = solver.heatFluxZGlobal[j];
        }
        for (int j = 0; j < solutions[i].numBound; ++j) {
            solutions[i].boundaryTemp[j]=solver.boundaryTemp[j];
            solutions[i].boundaryFlux[j]=solver.boundaryFlux[j];
        }
    }

}

void SolutionAll::_solve(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles, int num_proc, int world_rank,int Use_Backup, int Num_Max_Iter, int Use_Sythetic, int Use_Limiter, double error_temp_limit, double error_flux_limit)
{
if(world_rank==0)
    _print_out(mesh);

    exit(0);


}

void
SolutionAll::_Transient_BTE_Solver(DistributeMesh *mesh, BTEBoundaryCondition *bcs, BTEBand *bands, BTEAngle *angles,
                                   int num_proc, int world_rank, int Use_Backup, string Order, double error_temp_limit,
                                   double error_flux_limit, double deltaT, double totalT,int use_TDTR,double pulse_time,double repetition_frequency,double modulation_frequency,double xy_r) {
    Transient solver(&mesh->FourierMeshes, bcs, bands, angles, num_proc,
              world_rank,deltaT,totalT, use_TDTR, pulse_time, repetition_frequency, modulation_frequency, xy_r);
    for (int j = 0; j < solution.numCell; ++j) {
        solver.temperature[j]=solution.Temperature[j];
        solver.heatFluxXGlobal[j]=solution.heatFluxX[j];
        solver.heatFluxYGlobal[j]=solution.heatFluxY[j];
        solver.heatFluxZGlobal[j]=solution.heatFluxZ[j];
    }
    if(Order=="2")
    {
        solver.solve(Use_Backup,error_temp_limit,error_flux_limit,deltaT,totalT,use_TDTR, pulse_time, repetition_frequency, modulation_frequency, xy_r);
    }

    for (int j = 0; j < solution.numCell; ++j) {
        solution.Temperature[j] = solver.temperature[j];
        solution.heatFluxX[j] = solver.heatFluxXGlobal[j];
        solution.heatFluxY[j] = solver.heatFluxYGlobal[j];
        solution.heatFluxZ[j] = solver.heatFluxZGlobal[j];
    }
}

