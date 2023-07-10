// 2023/07
// This script is used for generating mesh file for 1D lines, 2D square/rectangle or 3D cubic geometries
// The mesh_control file specifies the dimension, length scale and the number of discrete meshes
// Please put the generated mesh file (inputmesh.txt) in the input folder

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
int main() {
    int Dimension;
    int LX=0;
    int LY=0;
    int LZ=0;
    int MeshNumX=0;
    int MeshNumY=0;
    int MeshNumZ=0;

    ifstream fin_const1("mesh_control");
    if (!fin_const1.is_open())
    {
        cout << "Error: need mesh_control" << endl;
        exit(0);
    }
    else
    {
        string str;
        char new_line;
        while(getline(fin_const1,str))
        {
            if(str.find("GeometryDimension")>=0 && str.find("GeometryDimension") < str.length())
            {
                fin_const1 >> Dimension;
            }
            if(str.find("X_Length")>=0 && str.find("X_Length") < str.length())
            {
                fin_const1 >> LX;
            }
            if(str.find("Y_Length")>=0 && str.find("Y_Length") < str.length())
            {
                fin_const1 >> LY;
            }
            if(str.find("Z_Length")>=0 && str.find("Z_Length") < str.length())
            {
                fin_const1 >> LZ;
            }
            if(str.find("X_MeshNumber")>=0 && str.find("X_MeshNumber") < str.length())
            {
                fin_const1 >> MeshNumX;
            }
            if(str.find("Y_MeshNumber")>=0 && str.find("Y_MeshNumber") < str.length())
            {
                fin_const1 >> MeshNumY;
            }
            if(str.find("Z_MeshNumber")>=0 && str.find("Z_MeshNumber") < str.length())
            {
                fin_const1 >> MeshNumZ;
            }
        }
    }
    fin_const1.close();

    //1D Line geometry
    if (Dimension==1)
    {
        ofstream output("inputmesh.txt");
        output<<"# COMSOL"<<endl<<"1 # sdim"<<endl<<MeshNumX+1<<" # number of mesh vertices"<<endl;
        output<<"0 # lowest mesh vertex index"<<endl<<endl;
        output<<"# Mesh vertex coordinates"<<endl;
        for(int i=0; i<MeshNumX+1; i++)
        {
            double temp;
            temp=i*(1.0/MeshNumX);
            output<<temp<<endl;
        }
        output<<endl<<"2 # number of element types"<<endl<<endl;

        //
        output<<"# Type #0"<<endl<<endl;
        output<<"1 # number of vertices per element"<<endl;
        output<<"2 # number of elements"<<endl;
        output<<"# Elements"<<endl;
        output<<0<<endl<<MeshNumX<<endl;

        output<<endl<<"2 # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        output<<0<<endl<<1<<endl;

        //
        output<<endl<<"# Type #1"<<endl<<endl;
        output<<"2 # number of vertices per element"<<endl;
        output<<MeshNumX<<" # number of elements"<<endl;
        output<<"# Elements"<<endl;
        for (int i=0; i<MeshNumX; i++)
        {
            output<<i<<" "<<i+1<<endl;
        }
        output<<endl;

        output<<MeshNumX<<" # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        for (int i=0; i<MeshNumX; i++)
        {
            output<<1<<endl;
        }
        output.close();

        //
        cout<<"Dimension = "<<Dimension<<endl;
        cout<<"Total number of discrete meshes = "<<MeshNumX<<endl;

    }
    else if (Dimension==2)
    {
        ofstream output("inputmesh.txt");
        output<<"# COMSOL"<<endl<<"2 # sdim"<<endl<<(MeshNumX+1)*(MeshNumY+1)<<" # number of mesh vertices"<<endl;
        output<<"0 # lowest mesh vertex index"<<endl<<endl;
        output<<"# Mesh vertex coordinates"<<endl;
        for (int i=0; i<MeshNumX+1; i++)
        {
            for(int j=0; j<MeshNumY+1; j++)
            {
                double temp1=i*(1.0/MeshNumX)*LX;
                double temp2=j*(1.0/MeshNumY)*LY;
                output<<temp1<<" "<<temp2<<endl;
            }
        }
        output<<endl<<"3 # number of element types"<<endl<<endl;

        //
        output<<"# Type #0"<<endl<<endl;
        output<<"# Type #1"<<endl<<endl;

        output<<"2 # number of vertices per element"<<endl;
        output<<MeshNumX*2+MeshNumY*2<<" # number of elements"<<endl;
        output<<"# Elements"<<endl;
        for(int i=0; i<MeshNumY; i++)
        {
            output<<i<<" "<<i+1<<endl; // nodes in boundary 1
        }
        for(int i=0; i<MeshNumX; i++)
        {
            output<<i*(MeshNumY+1)<<" "<<(i+1)*(MeshNumY+1)<<endl; //  nodes in boundary 2
        }
        for(int i=0; i<MeshNumX; i++)
        {
            output<<(i+1)*(MeshNumY+1)-1<<" "<<(i+2)*(MeshNumY+1)-1<<endl; // nodes in boundary 3
        }
        for(int i=0; i<MeshNumY; i++)
        {
            output<<MeshNumX*(MeshNumY+1)+i<<" "<<MeshNumX*(MeshNumY+1)+i+1<<endl; // nodes in boundary 4
        }

        output<<endl<<MeshNumX*2+MeshNumY*2<<" # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        for(int i=0; i<MeshNumY; i++)
        {
            output<<1-1<<endl;
        }
        for(int i=0; i<MeshNumX; i++)
        {
            output<<2-1<<endl;
        }
        for(int i=0; i<MeshNumX; i++)
        {
            output<<3-1<<endl;
        }
        for(int i=0; i<MeshNumY; i++)
        {
            output<<4-1<<endl;
        }

        //
        output<<endl<<"# Type #2"<<endl<<endl;

        output<<"4 # number of vertices per element"<<endl;
        output<<MeshNumX*MeshNumY<<" # number of elements"<<endl;
        output<<"# Elements"<<endl;
        for(int i=0; i<MeshNumX; i++)
        {
            for(int j=0; j<MeshNumY; j++)
            {
                output<<i*(MeshNumY+1)+j<<" "<<i*(MeshNumY+1)+j+1<<" "<<(i+1)*(MeshNumY+1)+j<<" "<<(i+1)*(MeshNumY+1)+j+1<<endl;
            }
        }
        output<<endl;

        output<<MeshNumX*MeshNumY<<" # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        for(int i=0; i<MeshNumX*MeshNumY; i++)
        {
            output<<1<<endl;
        }
        output.close();

        //
        cout<<"Dimension = "<<Dimension<<endl;
        cout<<"Total number of discrete meshes = "<<MeshNumX<<" x "<<MeshNumY<<" = "<<MeshNumY*MeshNumX<<endl;
    }
    else if (Dimension==3)
    {
        ofstream output("inputmesh.txt");
        output<<"# COMSOL"<<endl<<"3 # sdim"<<endl<<(MeshNumX+1)*(MeshNumY+1)*(MeshNumZ+1)<<" # number of mesh vertices"<<endl;
        output<<"0 # lowest mesh vertex index"<<endl<<endl;
        output<<"# Mesh vertex coordinates"<<endl;
        for(int i=1; i<MeshNumX+2; i++)
        {
            for(int j=1; j<MeshNumY+2; j++)
            {
                for(int k=1; k<MeshNumZ+2; k++)
                {
                    double temp1=(i-1)*(1.0/MeshNumX)*LX;
                    double temp2=(j-1)*(1.0/MeshNumY)*LY;
                    double temp3=(k-1)*(1.0/MeshNumZ)*LZ;
                    output<<temp1<<" "<<temp2<<" "<<temp3<<endl;
                }
            }
        }
        output<<endl<<"4 # number of element types"<<endl<<endl;

        //
        output<<"# Type #0"<<endl<<endl;
        output<<"# Type #1"<<endl<<endl;
        output<<"# Type #2"<<endl<<endl;

        output<<"4 quad # type name"<<endl;
        output<<MeshNumX*MeshNumY*2+MeshNumX*MeshNumZ*2+MeshNumZ*MeshNumY*2<<" # number of elements"<<endl;
        output<<"# Elements"<<endl;
        for (int i=1; i<MeshNumY+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)+j-1<<" "<<(i-1)*(MeshNumZ+1)+j<<" "<<i*(MeshNumZ+1)+j-1<<" "<<i*(MeshNumZ+1)+j<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j-1<<" "<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j-1<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumY+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)<<" "<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumY+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+MeshNumZ<<" "<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+MeshNumZ<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+MeshNumZ<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+MeshNumZ<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j-1+(MeshNumZ+1)*MeshNumY<<" "<<(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j+(MeshNumZ+1)*MeshNumY<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j-1+(MeshNumZ+1)*MeshNumY<<" "<<i*(MeshNumZ+1)*(MeshNumY+1)+j+(MeshNumZ+1)*MeshNumY<<endl;
            }
        }
        for (int i=1; i<MeshNumY+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<(i-1)*(MeshNumZ+1)+j-1+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX<<" "<<(i-1)*(MeshNumZ+1)+j+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX<<" "<<i*(MeshNumZ+1)+j-1+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX<<" "<<i*(MeshNumZ+1)+j+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX<<endl;
            }
        }

        output<<endl<<MeshNumX*MeshNumY*2+MeshNumX*MeshNumZ*2+MeshNumZ*MeshNumY*2<<" # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        for (int i=1; i<MeshNumY+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<1-1<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<2-1<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumY+1; j++)
            {
                output<<3-1<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumY+1; j++)
            {
                output<<4-1<<endl;
            }
        }
        for (int i=1; i<MeshNumX+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<5-1<<endl;
            }
        }
        for (int i=1; i<MeshNumY+1; i++)
        {
            for (int j=1; j<MeshNumZ+1; j++)
            {
                output<<6-1<<endl;
            }
        }

        //
        output<<endl<<"# Type #3"<<endl<<endl;

        output<<"3 hex # type name"<<endl;
        output<<"8 # number of vertices per element"<<endl;
        output<<MeshNumX*MeshNumY*MeshNumZ<<" # number of elements"<<endl;
        output<<"# Elements"<<endl;
        for(int i=1; i<MeshNumX+1; i++)
        {
            for(int j=1; j<MeshNumY+1; j++)
            {
                for(int k=1; k<MeshNumZ+1; k++)
                {
                    double temp1=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k-1;
                    double temp2=i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k-1;
                    double temp3=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k-1;
                    double temp4=i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k-1;
                    double temp5=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k;
                    double temp6=i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k;
                    double temp7=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k;
                    double temp8=i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k;
                    output<<temp1<<" "<<temp2<<" "<<temp3<<" "<<temp4<<" "<<temp5<<" "<<temp6<<" "<<temp7<<" "<<temp8<<endl;
                }
            }
        }

        output<<endl<<MeshNumX*MeshNumY*MeshNumZ<<" # number of geometric entity indices"<<endl;
        output<<"# Geometric entity indices"<<endl;
        for(int i=1; i<MeshNumX*MeshNumY*MeshNumZ+1; i++)
        {
            output<<1<<endl;
        }
        output.close();

        //
        cout<<"Dimension = "<<Dimension<<endl;
        cout<<"Total number of discrete meshes = "<<MeshNumX<<" x "<<MeshNumY<<" x "<<MeshNumZ<<" = "<<MeshNumY*MeshNumX*MeshNumZ<<endl;
    }

    cout<<endl<<"***********************************************************"<<endl;
    cout<<"***      Finish Generating Mesh File: inputmesh.txt     ***"<<endl;
    cout<<"***     Please put inputmesh.txt in the input folder    ***"<<endl;
    cout<<"***********************************************************"<<endl;



    return 0;
}
