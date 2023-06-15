//
// Created by yuehu on 2021/9/8.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "StaticBTESynthetic/StaticBTESynthetic.h"
#include "StaticFourier/StaticFourier.h"
#include "TransientFourier/TransientFourier.h"
#include "BTEMesh/Distribute.h"
#include "Solution/Solution.h"
using namespace std;

int main(int argc, char **argv)
{

#ifndef USE_GPU
    //PetscInitialize(nullptr, nullptr, (char *)nullptr, nullptr);
    MPI_Init(&argc,&argv);
    int num_proc, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    cout << "the process" << world_rank << " of" << num_proc << endl;
#else
    int num_proc, world_rank;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int device_count;
    int device_id;
    cudaGetDeviceCount(&device_count);
    device_id = world_rank % device_count;
    cudaSetDevice(device_id);
    if (world_rank == 0)
    {
        cout << num_proc << " process(es),   " << device_count << " device(s)" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Bind solver rank " << world_rank << " to device " << device_id << "." << endl;
#endif

    string Name_Mesh_Type; //necessary
    string Name_Heat_Type; //yufei adding not necessary
    string Name_Mesh_File; //necessary
    string Name_Band_File; //necessary
    string Name_Bc_File; //necessary
    string Name_Heat_File; //not necessary
    string Name_multiscale_File; //not necessary
    string Name_InitialTemp_File; //jiaxuan adding not necessary
    double L_x = 0; //necessary
    double L_y = 0; //necessary
    double L_z = 0; //necessary
    int Dimension_Geometry=-1; //necessary
    string Order = "2"; //not necessary
    string Method = "DOM"; //not necessary
    string Matrix_solver = "LU"; //not necessary
    int Num_Max_Iter=10000; //not necessary
    int Angle_method=2; //not necessary
    int Dimension_Material=3; //not necessary
    int Num_Phi=4; //not necessary
    int Num_Theta=4; //not necessary
    double error_temp_limit=1e-5;  //not necessary
    double error_flux_limit=1e-3;  //not necessary
    bool Use_Backup=false;
    double Uniform_heat=0; //not necessary
    double DeltaT=-1; //necessary for transient
    double TotalT=-1; //necessary for transient
    string State="Steady"; //not necessary
    double Tref=300;
    //for TDTR case
    int use_TDTR=0;  //not necessary
    double repetition_frequency=0;//necessary for transient TDTR
    double modulation_frequency=0;//necessary for transient TDTR
    double xy_r=0;//necessary for transient TDTR
    double pulse_time=0;//necessary for transient TDTR
    string ant="!";

    ifstream fin_const1("input/CONTROL");
    if (!fin_const1.is_open())
    {
        cout << "Error: need CONTROL" << endl;
        exit(0);
    }
    else
    {
        string str;
        char new_line;
        while(getline(fin_const1,str))
        {
            if(str.find("State")>=0 && str.find("State") < str.length() && str.find(ant) == str.npos)
            {
                    fin_const1 >> State;
            }
            if(str.find("Order") >= 0 && str.find("Order") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Order ;
            }
            if(str.find("IterativeScheme") >= 0 && str.find("IterativeScheme") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Method ;
            }
            if(str.find("MatrixSolver") >= 0 && str.find("MatrixSolver") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Matrix_solver;
            }
            if(str.find("Tref") >= 0 && str.find("Tref")< str.length() && str.find(ant) == str.npos)
            {
                fin_const1>>Tref;
            }
            if(str.find("MaxNumIter") >= 0 && str.find("MaxNumIter") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Num_Max_Iter;
            }
            if(str.find("ResidualTemp") >= 0 && str.find("ResidualTemp") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> error_temp_limit ;
            }
            if(str.find("ResidualFlux") >= 0 && str.find("ResidualFlux") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> error_flux_limit ;
            }
            if(str.find("Multiscale") >= 0 && str.find("Multiscale") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Name_multiscale_File ;
            }
            if(str.find("TimeStep") >= 0 && str.find("TimeStep") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> DeltaT ;
            }
            if(str.find("TotalTime") >= 0 && str.find("TotalTime") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> TotalT ;
            }
            if(str.find("InitialTemperatureFile") >= 0 && str.find("InitialTemperatureFile") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> Name_InitialTemp_File ;//jiaxuan
            }
            //TDTR
            if(str.find("IsTDTR") >= 0 && str.find("IsTDTR") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> use_TDTR; cout<<"use_TDTR: "<<use_TDTR<<endl;
            }
            if(str.find("PulseTime") >= 0 && str.find("PulseTime") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> pulse_time;cout<<"pulse_time: "<<pulse_time<<endl;
            }
            if(str.find("RepetitionFrequency") >= 0 && str.find("RepetitionFrequency") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> repetition_frequency;cout<<"repetition_frequency: "<<repetition_frequency<<endl;
            }
            if(str.find("ModulationFrequency") >= 0 && str.find("ModulationFrequency") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> modulation_frequency;cout<<"modulation_frequency: "<<modulation_frequency<<endl;
            }
            if(str.find("RadiusProbe") >= 0 && str.find("RadiusProbe") < str.length() && str.find(ant) == str.npos)
            {
                fin_const1 >> xy_r;cout<<"RadiusProbe: "<<xy_r<<endl;
            }
        }

    }
    fin_const1.close();
    if(State=="Transient"&&(TotalT==-1||DeltaT==-1))
    {
        cout<<"Error: Please provide DeltaTime and TotalTime for transient solver"<<endl;
        exit(0);
    }
    ifstream fin_const2("input/GEOMETRY");
    if (!fin_const2.is_open())
    {
        cout << "Error: need GEOMETRY" << endl;
        exit(0);
    }
    else
    {
        string str;
        char new_line;
        while (getline(fin_const2,str))
        {
            if(str.find("GeometryDimension") >= 0 && str.find("GeometryDimension") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Dimension_Geometry;
            }
            if(str.find("ScaleX") >= 0 && str.find("ScaleX") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> L_x ;
            }
            if(str.find("ScaleY") >= 0 && str.find("ScaleY") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> L_y;
            }
            if(str.find("ScaleZ") >= 0 && str.find("ScaleZ") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> L_z ;
            }
            if(str.find("UniformHeat") >= 0 && str.find("UniformHeat") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Uniform_heat;
            }
            if(str.find("BCFile") >= 0 && str.find("BCFile") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Name_Bc_File;
            }
            if(str.find("MeshType") >= 0 && str.find("MeshType") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Name_Mesh_Type;
            }
            if(str.find("MeshFile") >= 0 && str.find("MeshFile") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Name_Mesh_File;
            }
            if(str.find("HeatFileType") >= 0 && str.find("HeatFileType") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Name_Heat_Type;
            } // yufei adding
            if(str.find("HeatFile") >= 0 && str.find("HeatFile") < str.length() && str.find(ant) == str.npos)
            {
                fin_const2 >> Name_Heat_File;
            }
        }
    }
    fin_const2.close();

    if(Dimension_Geometry==-1)
    {
        cout<<"Error: Please provide DimensionGeometry"<<endl;
        exit(0);
    } else if(Dimension_Geometry==1&&L_x==0)
    {
        cout<<"Error: Please provide Lx for one dimensional problem"<<endl;
        exit(0);
    }
    else if(Dimension_Geometry==2&&(L_x==0||L_y==0))
    {
        cout<<"Error: Please provide Lx and Ly for two dimensional problem"<<endl;
        exit(0);
    }
    else if(Dimension_Geometry==3&&(L_x==0||L_y==0||L_z==0))
    {
        cout<<"Error: Please provide Lx, Ly and Lz for two dimensional problem"<<endl;
        exit(0);
    }


    ifstream inputangle("input/PHONON");
    string str;
    while(getline(inputangle, str))
    {
        if(str.find("MaterialDimension") >= 0 && str.find("MaterialDimension") < str.length() && str.find(ant) == str.npos)
        {
            inputangle >> Dimension_Material;
        }
        if(str.find("Angle_method") >= 0 && str.find("Angle_method") < str.length() && str.find(ant) == str.npos)
        {
            inputangle >> Angle_method;
        }
        if(str.find("Ntheta") >= 0 && str.find("Ntheta") < str.length() && str.find(ant) == str.npos)
        {
            inputangle >> Num_Theta;
        }
        if(str.find("Nphi") >= 0 && str.find("Nphi") < str.length() && str.find(ant) == str.npos)
        {
            inputangle >> Num_Phi;
        }
    }
    inputangle.close();
    if (world_rank == 0)
    {
        cout << "mesh " << Name_Mesh_File  << " bcfile " << Name_Bc_File << " L_x " << L_x << " ntheta " << Num_Theta << " nphi " << Num_Phi << endl;
    }

    ifstream bandFile("input/PHONON");
    auto bands = new BTEBand(bandFile, Dimension_Material);
    bandFile.close();

    ifstream bcFile(Name_Bc_File);
    auto bcs = new BTEBoundaryCondition(bcFile, bands);
    bcFile.close();

    BTEMesh *mesh;

    ifstream geofile(Name_Mesh_File);
    ifstream heatfile(Name_Heat_File);
    ifstream initialtempfile(Name_InitialTemp_File);//jiaxuan
    //mesh = new BTEMesh(Dimension_Geometry, geofile, L_x, L_y, L_z, bands, bcs, "default");
    //mesh->setMeshParams(bands);
    //mesh->setMeshParams1(bcs);
    //mesh->BTEMesh_heatin(heatfile, Uniform_heat);
    ifstream geofile1(Name_Mesh_File);
    DistributeMesh *distributeMesh;
    distributeMesh=new DistributeMesh(Dimension_Geometry, geofile1, L_x, L_y, L_z,
                                  Name_Mesh_Type,bands,bcs,heatfile, Uniform_heat, Name_Heat_Type,initialtempfile,Tref);  //yufei jiaxuan
    geofile.close();
    heatfile.close();
    initialtempfile.close();//jiaxuan



    BTEAngle *angles;
    angles = new BTEAngle(Num_Theta, Num_Phi, Dimension_Material, Angle_method);

#ifndef USE_GPU

    SolutionAll solutionAll(distributeMesh,bcs,bands,angles,num_proc,world_rank);
    solutionAll._set_initial(distributeMesh,bands,angles);
    if(State != "Transient")
    {
        solutionAll._Fourier_Solver(distributeMesh,bcs,bands,num_proc,world_rank);
        MPI_Barrier(MPI_COMM_WORLD);
        distributeMesh->_build_BTEMesh(Dimension_Geometry, L_x, L_y, L_z,bands,bcs,heatfile, Uniform_heat,Name_multiscale_File);
        MPI_Barrier(MPI_COMM_WORLD);
        solutionAll._BTE_Solver(distributeMesh,bcs,bands,angles,num_proc,world_rank,Use_Backup,
                                Num_Max_Iter,Order,Method,Matrix_solver,error_temp_limit,error_flux_limit);
        solutionAll._print_out(distributeMesh);
    } else
    {
        solutionAll._Transient_BTE_Solver(distributeMesh,bcs,bands,angles,num_proc,world_rank,Use_Backup,Order,error_temp_limit,error_flux_limit,DeltaT,TotalT
                ,use_TDTR,pulse_time,repetition_frequency,modulation_frequency,xy_r);
    }


//StaticFourier fourierStatic(mesh, bcs, bands,num_proc, world_rank);
    //fourierStatic._solve();

    //TransientFourier fourierTransient(mesh, bcs, bands,num_proc, world_rank);
    //fourierTransient._solve();
#else
    StaticBTESolver solver(mesh, bcs, bands, angles, num_proc, world_rank, Num_Theta, Num_Phi);
    solver.solve(Use_Backup, Num_Max_Iter, Use_Sythetic, Use_Limiter, error_temp_limit, error_flux_limit);
#endif
    if (world_rank == 0)
    {
        cout << "******************************" << endl << "Calculation Finished !" << endl << "******************************" << endl;
    }

    delete bands;
    delete bcs;
   // delete mesh;
    delete angles;

    MPI_Finalize();
}
