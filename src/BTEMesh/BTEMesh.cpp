//
// Created by yuehu on 2021/9/8.
//

// #include "StaticBTESynthetic/BTEMesh.h"
// #include "StaticBTESynthetic/BTEBand.h"
// #include "StaticBTESynthetic/BTEBoundaryCondition.h"
#include "BTEMesh/BTEMesh.h"
#include "BTEBand/BTEBand.h"
#include "BTEBoundaryCondition/BTEBoundaryCondition.h"
#include <sstream>
#include <algorithm> //yufei:for heat source
using namespace std;

BTEMesh::BTEMesh(int Dimension_Geometry,std::ifstream& inFile, double L_x, double L_y, double L_z, BTEBand *bands, BTEBoundaryCondition *bcs,std::string mesh_type) { 
    this->Dimension=Dimension_Geometry;
    if (Dimension_Geometry==1)
    {
        if (!inFile.is_open()) {
            std::cout << "DEBUG: Mesh file not open" << std::endl;
            exit(1);

        }
        string line;
        this->L_x = L_x;
        this->L_y = 0;
        this->L_z = 0;
        getline(inFile,line);
        int num;
        inFile>>num>>line;
        this->Elements.resize(num);
        this->Nodes.resize(num+1);
        for (int i = 0; i <num+1 ; ++i) {

            Nodes[i].x=L_x/num*i;

        }

        for (int i = 0; i <num ; ++i) {
            Elements[i].vertexes.push_back(i);
            Elements[i].vertexes.push_back(i+1);
            Elements[i].center.x=1.0/2*(Nodes[Elements[i].vertexes[0]].x+Nodes[Elements[i].vertexes[1]].x);
            Elements[i].faces.resize(2);
            Elements[i].volume=L_x/num;
            Elements[i].faces[0].norm.x=-1;
            Elements[i].faces[1].norm.x=1;
            Elements[i].faces[0].area=1;
            Elements[i].faces[1].area=1;
            Elements[i].faces[0].vertexes.push_back(i);
            Elements[i].faces[1].vertexes.push_back(i+1);
            Elements[i].faces[0].center.x=Nodes[Elements[i].vertexes[0]].x;
            Elements[i].faces[1].center.x=Nodes[Elements[i].vertexes[1]].x;
            Elements[i].faces[0].direction.y=1;
            Elements[i].faces[1].direction.y=-1;
        }
        getline(inFile,line);
        int num_geo;
        inFile>>num_geo>>line;
        getline(inFile,line);
        vector<double> distribution(num_geo,0);
        for (int i = 0; i <num_geo ; ++i) {
            getline(inFile,line);
            std::stringstream ss(line);
            int temp;
            ss>>temp>>distribution[i];
        }
        for (int i = 0; i < num ; ++i) {
            if (Elements[i].center.x<distribution[0]*L_x)
            {
                Elements[i].index=0;
                Elements[i].matter=bands->geo_matter_index[0];
            }
            for (int j = 1; j < num_geo ; ++j) {
                if (Elements[i].center.x<distribution[j]*L_x&&Elements[i].center.x>distribution[j-1]*L_x)
                {
                    Elements[i].index=j;
                    if(j>=bands->geo_matter_index.size())
                    {
                        cout<<"Error: some region does not have material, check PHONON and GEOMETRY"<<endl;
                        exit(0);
                    }
                    Elements[i].matter=bands->geo_matter_index[j];
                }
            }
        }

        vector<double> boundx;
        boundx.push_back(0);
        for (int i = 0; i <num_geo ; ++i) {
            boundx.push_back(distribution[i]*L_x);
        }
        for (int i = 0; i < num_geo+1; ++i) {
            Boundary bound;
            bound.index=i;
            for (int j = 0; j <num ; ++j) {
                for (int k = 0; k <2 ; ++k) {
                    if (abs(Nodes[Elements[j].faces[k].vertexes[0]].x-boundx[i])<L_x/1e8)
                    {
                        bound.cellindex.push_back(j);
                        bound.faceindex.push_back(k);
                    }
                }
            }
            bound.macro_Temp.resize(bound.cellindex.size());
            bound.emit_temp.resize(bound.cellindex.size());

            for (int j = 0; j <bound.cellindex.size() ; ++j) {
                bound.macro_Temp[j].matter=Elements[bound.cellindex[j]].matter;
                bound.macro_Temp[j].Temp=0;
            }
            Boundaries.push_back(bound);
        }
        for (int i = 0; i <num ; ++i) {
            Nodes[Elements[i].vertexes[0]].cells.push_back(i);
            Nodes[Elements[i].vertexes[1]].cells.push_back(i);
        }
        for (int i = 0; i <num ; ++i) {
            for (int j = 0; j <2 ; ++j) {
                Elements[i].faces[j].neighbor=-1;
                for (int k = 0; k < Nodes[Elements[i].faces[j].vertexes[0]].cells.size() ; ++k) {
                    int index=Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index!=i){
                        if (Elements[i].faces[j].vertexes[0]==Elements[index].vertexes[0]||Elements[i].faces[j].vertexes[0]==Elements[index].vertexes[1])
                        {
                            Elements[i].faces[j].neighbor=index;
                        }
                    }
                }
                Elements[i].faces[j].bound=-1;

            }
        }

    }
    else if (Dimension_Geometry==2)
    {
        nedge=4;
        vector<vector<int>> Boundaryindex;
        if (!inFile.is_open()) {
            std::cout << "DEBUG: Mesh file not open" << std::endl;
            exit(1);

        }
        string line;
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = 0;
        string strmesh;
        int numnode;
        while(getline(inFile, strmesh)){
            if(strmesh.find("# number of mesh")>0
            && strmesh.find("# number of mesh")<strmesh.length()){
                numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                break;
            }
        }
        Nodes.resize(numnode);
        for (int i = 0; i < 3; i++)
            getline(inFile, strmesh);
        for (int i = 0; i < numnode; i++) {
            getline(inFile, strmesh);
            stringstream sss;
            sss.str(strmesh);
            string coord_x, coord_y;
            sss >> coord_x >> coord_y;
            Nodes[i].x = strtod(coord_x.c_str(), NULL) * L_x;
            Nodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;
            Nodes[i].z=0;
        }
        int numoftype;
        while(getline(inFile, strmesh)){
            if(strmesh.find("# number of element types")>0
            && strmesh.find("# number of element types")<strmesh.length()){
                numoftype = stoi(strmesh.substr(0, strmesh.find(" ")));
                break;
            }
        }
        while (numoftype>0)
        {
            numoftype--;
            int numvert;
            while(getline(inFile, strmesh)){
                if(strmesh.find("# number of vertices per element")>0
                && strmesh.find("# number of vertices per element")<strmesh.length()||strmesh.find("# number of nodes per element")>0
                && strmesh.find("# number of nodes per element")<strmesh.length()){
                    numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            if (numvert==2)
            {
                while(getline(inFile, strmesh)){
                    if(strmesh.find("# number of elements")>0
                    && strmesh.find("# number of elements")<strmesh.length()){
                        numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                Boundaries.resize(numvert);
                getline(inFile, strmesh);
                Boundaryindex.resize(numvert);
                for (int i=0;i<numvert;i++)
                {
                    Boundaryindex[i].resize(2);
                    inFile>>Boundaryindex[i][0]>>Boundaryindex[i][1];

                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i=0;i<numvert;i++)
                {
                    inFile>>Boundaries[i].index;
                }
            }
            else if (numvert==3)
            {
                while(getline(inFile, strmesh)){
                    if(strmesh.find("# number of elements")>0
                    && strmesh.find("# number of elements")<strmesh.length()){
                        numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i <numvert ; ++i) {
                    Cell c1;
                    c1.vertexes.resize(3);
                    inFile>>c1.vertexes[0]>>c1.vertexes[1]>>c1.vertexes[2];
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i=Elements.size()-numvert;i<Elements.size();i++)
                {
                    inFile>>Elements[i].index;
                }

            }
            else if (numvert==4)
            {
                while(getline(inFile, strmesh)){
                    if(strmesh.find("# number of elements")>0
                    && strmesh.find("# number of elements")<strmesh.length()){
                        numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i <numvert ; ++i) {
                    Cell c1;
                    c1.vertexes.resize(4);
                    inFile>>c1.vertexes[0]>>c1.vertexes[1]>>c1.vertexes[2]>>c1.vertexes[3];
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i=Elements.size()-numvert;i<Elements.size();i++)
                {
                    inFile>>Elements[i].index;
                }

            }

        }
        //set_index and matter
        for (int i = 0; i <Elements.size()  ; ++i) {
            Elements[i].index=Elements[i].index-1;
            if(Elements[i].index>=bands->geo_matter_index.size())
            {
                cout<<"Error: some region does not have material, check PHONON and GEOMETRY"<<endl;
                exit(0);
            }
            Elements[i].matter=bands->geo_matter_index[Elements[i].index];
        }

        //get_volume
        for (int i = 0; i <Elements.size()  ; ++i) {
            if (Elements[i].vertexes.size()==3)
            {
                vector<Point> tri;
                tri.push_back(Nodes[Elements[i].vertexes[0]]);
                tri.push_back(Nodes[Elements[i].vertexes[1]]);
                tri.push_back(Nodes[Elements[i].vertexes[2]]);
                Elements[i].volume=get_tri_area(tri);
                Elements[i].center=((Nodes[Elements[i].vertexes[0]]+Nodes[Elements[i].vertexes[1]])+Nodes[Elements[i].vertexes[2]])/3;
            } else if (Elements[i].vertexes.size()==4)
            {
                vector<Point> tri (3);
                tri[0]=Nodes[Elements[i].vertexes[0]];
                tri[1]=Nodes[Elements[i].vertexes[1]];
                tri[2]=Nodes[Elements[i].vertexes[2]];
                Elements[i].volume=get_tri_area(tri);
                tri[0]=Nodes[Elements[i].vertexes[1]];
                tri[1]=Nodes[Elements[i].vertexes[2]];
                tri[2]=Nodes[Elements[i].vertexes[3]];
                Elements[i].volume+=get_tri_area(tri);
                Elements[i].center=((Nodes[Elements[i].vertexes[0]]+Nodes[Elements[i].vertexes[1]])+Nodes[Elements[i].vertexes[2]]+Nodes[Elements[i].vertexes[3]])/4;
            }
        }
        //get_faces
        for (int i = 0; i <Elements.size()  ; ++i) {
            if (Elements[i].vertexes.size()==4)
            {
                Elements[i].faces.resize(4);
                for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                    Elements[i].faces[j].vertexes.resize(2);

                }
                Elements[i].faces[0].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1]=Elements[i].vertexes[1];

                Elements[i].faces[1].vertexes[0]=Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1]=Elements[i].vertexes[3];

                Elements[i].faces[2].vertexes[0]=Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[1]=Elements[i].vertexes[2];

                Elements[i].faces[3].vertexes[0]=Elements[i].vertexes[2];
                Elements[i].faces[3].vertexes[1]=Elements[i].vertexes[0];

                for (int j = 0; j <4 ; ++j) {
                    Elements[i].faces[j].area=get_distance(Nodes[Elements[i].faces[j].vertexes[0]],Nodes[Elements[i].faces[j].vertexes[1]]);
                    Elements[i].faces[j].center=(Nodes[Elements[i].faces[j].vertexes[0]]+Nodes[Elements[i].faces[j].vertexes[1]])/2;
                    vec a=(Nodes[Elements[i].faces[j].vertexes[0]]-Nodes[Elements[i].faces[j].vertexes[1]])/Elements[i].faces[j].area;
                    Elements[i].faces[j].norm.x=a.y;
                    Elements[i].faces[j].norm.y=-a.x;
                    Elements[i].faces[j].norm.z=0;
                    vec a1=Elements[i].center-Elements[i].faces[j].center;
                    if (a1*Elements[i].faces[j].norm>0)
                    {
                        Elements[i].faces[j].norm.x=-Elements[i].faces[j].norm.x;
                        Elements[i].faces[j].norm.y=-Elements[i].faces[j].norm.y;
                    }
                }
            } else if (Elements[i].vertexes.size()==3)
            {
                Elements[i].faces.resize(3);
                for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                    Elements[i].faces[j].vertexes.resize(2);
                    Elements[i].faces[j].vertexes[0]=Elements[i].vertexes[(j+1)%3];
                    Elements[i].faces[j].vertexes[1]=Elements[i].vertexes[(j+2)%3];
                    Elements[i].faces[j].area=get_distance(Nodes[Elements[i].faces[j].vertexes[0]],Nodes[Elements[i].faces[j].vertexes[1]]);
                    Elements[i].faces[j].center=(Nodes[Elements[i].faces[j].vertexes[0]]+Nodes[Elements[i].faces[j].vertexes[1]])/2;
                    vec a=(Nodes[Elements[i].faces[j].vertexes[0]]-Nodes[Elements[i].faces[j].vertexes[1]])/Elements[i].faces[j].area;
                    Elements[i].faces[j].norm.x=a.y;
                    Elements[i].faces[j].norm.y=-a.x;
                    Elements[i].faces[j].norm.z=0;
                    vec a1=Elements[i].center-Elements[i].faces[j].center;
                    if (a1*Elements[i].faces[j].norm>0)
                    {
                        Elements[i].faces[j].norm.x=-Elements[i].faces[j].norm.x;
                        Elements[i].faces[j].norm.y=-Elements[i].faces[j].norm.y;
                    }
                }


            }
        }
        for (int i = 0; i < Elements.size() ; ++i) {
            for (int j = 0; j <Elements[i].vertexes.size() ; ++j) {
                Nodes[Elements[i].vertexes[j]].cells.push_back(i);
            }
        }

        //get_neighbor
        for (int i = 0; i <Elements.size() ; ++i) {
            for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                Elements[i].faces[j].neighbor=-1;
                for (int k = 0; k < Nodes[Elements[i].faces[j].vertexes[0]].cells.size() ; ++k) {
                    int index=Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index!=i){
                        if (ishave(Elements[i].faces[j].vertexes[0],Elements[index].vertexes)==1&&
                        ishave(Elements[i].faces[j].vertexes[1],Elements[index].vertexes)==1 )
                        {
                            Elements[i].faces[j].neighbor=index;
                        }
                    }
                }
                Elements[i].faces[j].bound=-1;

            }
        }
        for (int i = 0; i <Boundaryindex.size() ; ++i) {
            for (int j = 0; j <Nodes[Boundaryindex[i][0]].cells.size() ; ++j) {
                int index=Nodes[Boundaryindex[i][0]].cells[j];
                for (int k = 0; k <Elements[index].faces.size() ; ++k) {
                    if (Boundaryindex[i][1] ==  Elements[index].faces[k].vertexes[0] ||
                    Boundaryindex[i][1] == Elements[index].faces[k].vertexes[1] )
                    {
                        if (Boundaryindex[i][0] ==  Elements[index].faces[k].vertexes[0] ||
                        Boundaryindex[i][0] == Elements[index].faces[k].vertexes[1])
                        {
                            Boundaries[i].cellindex.push_back(index);
                            Boundaries[i].faceindex.push_back(k);
                        }

                    }
                }
            }
            Boundaries[i].macro_Temp.resize(Boundaries[i].cellindex.size());
            Boundaries[i].emit_temp.resize(Boundaries[i].cellindex.size());
            for (int j = 0; j <Boundaries[i].cellindex.size() ; ++j) {
                Boundaries[i].macro_Temp[j].matter=Elements[Boundaries[i].cellindex[j]].matter;
                Boundaries[i].macro_Temp[j].Temp=0;
            }
        }

    }
    else if (Dimension_Geometry==3)
    {
        nedge=6;
        vector<vector<int>> Boundaryindex;
        if (!inFile.is_open()) {
            std::cout << "DEBUG: Mesh file not open" << std::endl;
            exit(1);

        }
        string line;
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
        string strmesh;
        int numnode;
        while(getline(inFile, strmesh)){
            if(strmesh.find("# number of mesh")>0
            && strmesh.find("# number of mesh")<strmesh.length()){
                numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                break;
            }
        }
        Nodes.resize(numnode);
        for (int i = 0; i < 3; i++)
            getline(inFile, strmesh);
        for (int i = 0; i < numnode; i++) {
            getline(inFile, strmesh);
            stringstream sss;
            sss.str(strmesh);
            string coord_x, coord_y, coord_z;
            sss >> coord_x >> coord_y >> coord_z;
            Nodes[i].x = strtod(coord_x.c_str(), NULL) * L_x;
            Nodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;
            Nodes[i].z = strtod(coord_z.c_str(), NULL) * L_z;
        }
        int numoftype;
        while (getline(inFile, strmesh)) {
            if (strmesh.find("# number of element types") > 0
            && strmesh.find("# number of element types") < strmesh.length()) {
                numoftype = stoi(strmesh.substr(0, strmesh.find(" ")));
                break;
            }
        }
        while (numoftype > 0) {
            numoftype--;
            string nameoftype;
            while (getline(inFile, strmesh)) {
                if (strmesh.find("# type name") > 0
                && strmesh.find("# type name") < strmesh.length()) {
                    stringstream sss;
                    sss.str(strmesh);
                    int three;
                    sss >> three >> nameoftype;
                    //nameoftype = strmesh.substr(2, strmesh.find(" "));
                    //cout << nameoftype << endl;
                    break;
                }
            }
            if (nameoftype == "quad") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    vector<int> bc;
                    bc.resize(4);
                    inFile >> bc[0] >> bc[1] >> bc[2] >> bc[3];
                    Boundaryindex.push_back(bc);
                    Boundary bound;
                    Boundaries.push_back(bound);
                    //fin_mesh>>boundaries[i].p1>>boundaries[i].p2;
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Boundaries.size() - numc; i < Boundaries.size(); i++) {
                    inFile >> Boundaries[i].index;
                }
            }
            else if (nameoftype == "tri") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    vector<int> bc;
                    bc.resize(3);
                    inFile >> bc[0] >> bc[1] >> bc[2] ;
                    Boundaryindex.push_back(bc);
                    Boundary bound;
                    Boundaries.push_back(bound);
                    //fin_mesh>>boundaries[i].p1>>boundaries[i].p2;
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Boundaries.size() - numc; i < Boundaries.size(); i++) {
                    inFile >> Boundaries[i].index;
                }
            }
            else if (nameoftype == "hex") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    Cell c1;
                    c1.vertexes.resize(8);
                    inFile >> c1.vertexes[0] >> c1.vertexes[1] >> c1.vertexes[2] >> c1.vertexes[3] >> c1.vertexes[4]
                    >> c1.vertexes[5] >> c1.vertexes[6] >> c1.vertexes[7];
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Elements.size() - numc; i < Elements.size(); i++) {
                    inFile>> Elements[i].index;
                }
            }
            else if (nameoftype == "tet") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    Cell c1;
                    c1.vertexes.resize(4);
                    inFile >> c1.vertexes[0] >> c1.vertexes[1] >> c1.vertexes[2] >> c1.vertexes[3] ;
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Elements.size() - numc; i < Elements.size(); i++) {
                    inFile>> Elements[i].index;
                }
            }
            else if (nameoftype == "pyr") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    Cell c1;
                    c1.vertexes.resize(5);
                    inFile >> c1.vertexes[0] >> c1.vertexes[1] >> c1.vertexes[2] >> c1.vertexes[3] >> c1.vertexes[4];
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Elements.size() - numc; i < Elements.size(); i++) {
                    inFile>> Elements[i].index;
                }
            }
            else if (nameoftype == "prism") {
                int numc;
                while (getline(inFile, strmesh)) {
                    if (strmesh.find("# number of elements") > 0
                    && strmesh.find("# number of elements") < strmesh.length()) {
                        numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                getline(inFile, strmesh);
                for (int i = 0; i < numc; i++) {
                    Cell c1;
                    c1.vertexes.resize(6);
                    inFile >> c1.vertexes[0] >> c1.vertexes[1] >> c1.vertexes[2] >> c1.vertexes[3] >> c1.vertexes[4]
                    >> c1.vertexes[5];
                    Elements.push_back(c1);
                }
                for (int i = 0; i < 4; i++)
                    getline(inFile, strmesh);
                for (int i = Elements.size() - numc; i < Elements.size(); i++) {
                    inFile>> Elements[i].index;
                }
            }



        }
        //set_index and matter
        for (int i = 0; i <Elements.size()  ; ++i) {
            Elements[i].index=Elements[i].index-1;
            if(Elements[i].index>=bands->geo_matter_index.size())
            {
                cout<<"Error: some region does not have material, check PHONON and GEOMETRY"<<endl;
                exit(0);
            }
            Elements[i].matter=bands->geo_matter_index[Elements[i].index];
        }
        //get_volume

        for (int i = 0; i <Elements.size(); ++i) {
            Elements[i].center.x=0;
            Elements[i].center.y=0;
            Elements[i].center.z=0;
            for (int j=0;j<Elements[i].vertexes.size();j++) {
                Elements[i].center=Elements[i].center+Nodes[Elements[i].vertexes[j]]/Elements[i].vertexes.size();
            }
            if (Elements[i].vertexes.size()==4)
            {
                vector<Point> a;
                a.resize(4);
                for (int j = 0; j <Elements[i].vertexes.size() ; ++j) {
                    a[j]=Nodes[Elements[i].vertexes[j]];
                }
                Elements[i].volume=get_tet_volume(a);
            } else if (Elements[i].vertexes.size()==5)
            {
                Elements[i].volume=0;
                vector<Point> a;
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[0]];
                a[1]=Nodes[Elements[i].vertexes[1]];
                a[2]=Nodes[Elements[i].vertexes[2]];
                a[3]=Nodes[Elements[i].vertexes[4]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[1]];
                a[1]=Nodes[Elements[i].vertexes[2]];
                a[2]=Nodes[Elements[i].vertexes[3]];
                a[3]=Nodes[Elements[i].vertexes[4]];
                Elements[i].volume+=get_tet_volume(a);
            }else if (Elements[i].vertexes.size()==6)
            {
                Elements[i].volume=0;
                vector<Point> a;
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[0]];
                a[1]=Nodes[Elements[i].vertexes[1]];
                a[2]=Nodes[Elements[i].vertexes[2]];
                a[3]=Nodes[Elements[i].vertexes[3]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[1]];
                a[1]=Nodes[Elements[i].vertexes[2]];
                a[2]=Nodes[Elements[i].vertexes[3]];
                a[3]=Nodes[Elements[i].vertexes[4]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[2]];
                a[1]=Nodes[Elements[i].vertexes[3]];
                a[2]=Nodes[Elements[i].vertexes[4]];
                a[3]=Nodes[Elements[i].vertexes[5]];
                Elements[i].volume+=get_tet_volume(a);
            }else if (Elements[i].vertexes.size()==8)
            {
                Elements[i].volume=0;
                vector<Point> a;
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[0]];
                a[1]=Nodes[Elements[i].vertexes[1]];
                a[2]=Nodes[Elements[i].vertexes[2]];
                a[3]=Nodes[Elements[i].vertexes[4]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[1]];
                a[1]=Nodes[Elements[i].vertexes[2]];
                a[2]=Nodes[Elements[i].vertexes[3]];
                a[3]=Nodes[Elements[i].vertexes[7]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[2]];
                a[1]=Nodes[Elements[i].vertexes[4]];
                a[2]=Nodes[Elements[i].vertexes[6]];
                a[3]=Nodes[Elements[i].vertexes[7]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[1]];
                a[1]=Nodes[Elements[i].vertexes[2]];
                a[2]=Nodes[Elements[i].vertexes[4]];
                a[3]=Nodes[Elements[i].vertexes[7]];
                Elements[i].volume+=get_tet_volume(a);
                a.resize(4);
                a[0]=Nodes[Elements[i].vertexes[1]];
                a[1]=Nodes[Elements[i].vertexes[4]];
                a[2]=Nodes[Elements[i].vertexes[5]];
                a[3]=Nodes[Elements[i].vertexes[7]];
                Elements[i].volume+=get_tet_volume(a);
            }
        }
        //set_face
        for (int i = 0; i <Elements.size(); ++i)
        {
            if (Elements[i].vertexes.size()==4)
            {
                Elements[i].faces.resize(4);
                for (int j = 0; j < 4 ; ++j) {
                    Elements[i].faces[j].vertexes.resize(3);
                }
                Elements[i].faces[0].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2]=Elements[i].vertexes[2];

                Elements[i].faces[1].vertexes[0]=Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1]=Elements[i].vertexes[2];
                Elements[i].faces[1].vertexes[2]=Elements[i].vertexes[3];

                Elements[i].faces[2].vertexes[0]=Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[1]=Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[2]=Elements[i].vertexes[0];

                Elements[i].faces[3].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[2]=Elements[i].vertexes[3];
            }
            else if (Elements[i].vertexes.size()==5)
            {
                Elements[i].faces.resize(5);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2]=Elements[i].vertexes[2];
                Elements[i].faces[0].vertexes[3]=Elements[i].vertexes[3];

                Elements[i].faces[1].vertexes.resize(3);
                Elements[i].faces[1].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[1].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[2]=Elements[i].vertexes[4];

                Elements[i].faces[2].vertexes.resize(3);
                Elements[i].faces[2].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1]=Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[2]=Elements[i].vertexes[4];

                Elements[i].faces[3].vertexes.resize(3);
                Elements[i].faces[3].vertexes[0]=Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[1]=Elements[i].vertexes[3];
                Elements[i].faces[3].vertexes[2]=Elements[i].vertexes[4];

                Elements[i].faces[4].vertexes.resize(3);
                Elements[i].faces[4].vertexes[0]=Elements[i].vertexes[2];
                Elements[i].faces[4].vertexes[1]=Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[2]=Elements[i].vertexes[4];
            }
            else if (Elements[i].vertexes.size()==6)
            {
                Elements[i].faces.resize(5);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2]=Elements[i].vertexes[3];
                Elements[i].faces[0].vertexes[3]=Elements[i].vertexes[4];

                Elements[i].faces[1].vertexes.resize(4);
                Elements[i].faces[1].vertexes[0]=Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1]=Elements[i].vertexes[2];
                Elements[i].faces[1].vertexes[2]=Elements[i].vertexes[4];
                Elements[i].faces[1].vertexes[3]=Elements[i].vertexes[5];

                Elements[i].faces[2].vertexes.resize(4);
                Elements[i].faces[2].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1]=Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[2]=Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[3]=Elements[i].vertexes[5];

                Elements[i].faces[3].vertexes.resize(3);
                Elements[i].faces[3].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[2]=Elements[i].vertexes[2];

                Elements[i].faces[4].vertexes.resize(3);
                Elements[i].faces[4].vertexes[0]=Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[1]=Elements[i].vertexes[4];
                Elements[i].faces[4].vertexes[2]=Elements[i].vertexes[5];
            }
            else if (Elements[i].vertexes.size()==8)
            {
                Elements[i].faces.resize(6);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2]=Elements[i].vertexes[2];
                Elements[i].faces[0].vertexes[3]=Elements[i].vertexes[3];

                Elements[i].faces[1].vertexes.resize(4);
                Elements[i].faces[1].vertexes[0]=Elements[i].vertexes[4];
                Elements[i].faces[1].vertexes[1]=Elements[i].vertexes[5];
                Elements[i].faces[1].vertexes[2]=Elements[i].vertexes[6];
                Elements[i].faces[1].vertexes[3]=Elements[i].vertexes[7];

                Elements[i].faces[2].vertexes.resize(4);
                Elements[i].faces[2].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1]=Elements[i].vertexes[1];
                Elements[i].faces[2].vertexes[2]=Elements[i].vertexes[4];
                Elements[i].faces[2].vertexes[3]=Elements[i].vertexes[5];

                Elements[i].faces[3].vertexes.resize(4);
                Elements[i].faces[3].vertexes[0]=Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1]=Elements[i].vertexes[2];
                Elements[i].faces[3].vertexes[2]=Elements[i].vertexes[4];
                Elements[i].faces[3].vertexes[3]=Elements[i].vertexes[6];

                Elements[i].faces[4].vertexes.resize(4);
                Elements[i].faces[4].vertexes[0]=Elements[i].vertexes[1];
                Elements[i].faces[4].vertexes[1]=Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[2]=Elements[i].vertexes[5];
                Elements[i].faces[4].vertexes[3]=Elements[i].vertexes[7];

                Elements[i].faces[5].vertexes.resize(4);
                Elements[i].faces[5].vertexes[0]=Elements[i].vertexes[2];
                Elements[i].faces[5].vertexes[1]=Elements[i].vertexes[3];
                Elements[i].faces[5].vertexes[2]=Elements[i].vertexes[6];
                Elements[i].faces[5].vertexes[3]=Elements[i].vertexes[7];
            }

        }
        //set_face_area_center
        for (int i = 0; i <Elements.size(); ++i) {
            for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                Elements[i].faces[j].center.x=0;
                Elements[i].faces[j].center.y=0;
                Elements[i].faces[j].center.z=0;
                for (int k = 0; k <Elements[i].faces[j].vertexes.size() ; ++k) {
                    Elements[i].faces[j].center=Elements[i].faces[j].center+Nodes[Elements[i].faces[j].vertexes[k]]/Elements[i].faces[j].vertexes.size();
                }
                if (Elements[i].faces[j].vertexes.size()==4)
                {
                    Elements[i].faces[j].area=0;
                    vector<Point> tri;
                    tri.resize(3);
                    tri[0]=Nodes[Elements[i].faces[j].vertexes[0]];
                    tri[1]=Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[2]=Nodes[Elements[i].faces[j].vertexes[2]];
                    Elements[i].faces[j].area+=get_tri_area(tri);
                    //cout<<get_tri_area(tri)<<endl;
                    tri.resize(3);
                    tri[0]=Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[1]=Nodes[Elements[i].faces[j].vertexes[2]];
                    tri[2]=Nodes[Elements[i].faces[j].vertexes[3]];
                    Elements[i].faces[j].area+=get_tri_area(tri);
                    //cout<<get_tri_area(tri)<<endl;
                } else
                {
                    vector<Point> tri;
                    tri.resize(3);
                    tri[0]=Nodes[Elements[i].faces[j].vertexes[0]];
                    tri[1]=Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[2]=Nodes[Elements[i].faces[j].vertexes[2]];
                    Elements[i].faces[j].area=get_tri_area(tri);
                }

            }
        }

        //get_face_norm
        for (int i = 0; i <Elements.size() ; ++i) {
            for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                vec a,b;
                a=Nodes[Elements[i].faces[j].vertexes[0]]-Nodes[Elements[i].faces[j].vertexes[1]];
                b=Nodes[Elements[i].faces[j].vertexes[0]]-Nodes[Elements[i].faces[j].vertexes[2]];
                Elements[i].faces[j].norm= cross_product(a,b);
                double length= get_length(Elements[i].faces[j].norm);
                Elements[i].faces[j].norm=Elements[i].faces[j].norm/length;
                a=Elements[i].faces[j].center-Elements[i].center;
                if (dot_product(a,Elements[i].faces[j].norm)<0)
                {
                    Elements[i].faces[j].norm.x=-Elements[i].faces[j].norm.x;
                    Elements[i].faces[j].norm.y=-Elements[i].faces[j].norm.y;
                    Elements[i].faces[j].norm.z=-Elements[i].faces[j].norm.z;
                }
            }
        }
        //
        for (int i = 0; i < Elements.size() ; ++i) {
            for (int j = 0; j <Elements[i].vertexes.size() ; ++j) {
                Nodes[Elements[i].vertexes[j]].cells.push_back(i);
            }
        }
        for (int i = 0; i <Elements.size() ; ++i) {
            for (int j = 0; j <Elements[i].faces.size() ; ++j) {
                Elements[i].faces[j].neighbor=-1;
                Elements[i].faces[j].bound=-1;
                for (int k = 0; k <Nodes[Elements[i].faces[j].vertexes[0]].cells.size() ; ++k) {
                    int index=Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index!=i)
                    {
                        vector<int> vertexes=Elements[index].vertexes;
                        if (ishave(Elements[i].faces[j].vertexes[1],vertexes)==1&&ishave(Elements[i].faces[j].vertexes[2],vertexes)==1&&ishave(Elements[i].faces[j].vertexes[0],vertexes)==1)
                        {
                            Elements[i].faces[j].neighbor=index;
                        }
                    }


                }
            }
        }
        for (int i = 0; i <Boundaryindex.size() ; ++i) {
            for (int j = 0; j <Nodes[Boundaryindex[i][0]].cells.size() ; ++j) {
                int index=Nodes[Boundaryindex[i][0]].cells[j];
                for (int k = 0; k <Elements[index].faces.size() ; ++k) {
                    if (ishave(Boundaryindex[i][0],Elements[index].faces[k].vertexes)==1&&ishave(Boundaryindex[i][1],Elements[index].faces[k].vertexes)==1&&ishave(Boundaryindex[i][2],Elements[index].faces[k].vertexes)==1)
                    {
                        Boundaries[i].cellindex.push_back(index);
                        Boundaries[i].faceindex.push_back(k);
                    }
                }
            }
            Boundaries[i].macro_Temp.resize(Boundaries[i].cellindex.size());
            Boundaries[i].emit_temp.resize(Boundaries[i].cellindex.size());
            for (int j = 0; j <Boundaries[i].cellindex.size() ; ++j) {
                Boundaries[i].macro_Temp[j].matter=Elements[Boundaries[i].cellindex[j]].matter;
                Boundaries[i].macro_Temp[j].Temp=0;
            }
        }


    }
    vector<Boundary> Boundaries1=Boundaries;
    Boundaries.resize(0);
    for (int i = 0; i <Boundaries1.size() ; ++i) {
        for (int j = 0; j < bcs->boundaryConditions.size() ; ++j) {
            if (Boundaries1[i].index==bcs->boundaryConditions[j].index)
            {
                Boundaries1[i].index=j;
                Boundaries.push_back(Boundaries1[i]);
            }
        }
    }
    for (int i = 0; i <Boundaries.size() ; ++i) {
        for (int j = 0; j <Boundaries[i].cellindex.size() ; ++j) {
            Elements[Boundaries[i].cellindex[j]].faces[Boundaries[i].faceindex[j]].bound=i;
        }
    }
    for (int i = 0; i <Elements.size() ; ++i) {
        for (int j = 0; j <Elements[i].faces.size() ; ++j) {
            if (Elements[i].faces[j].bound==-1&&Elements[i].faces[j].neighbor==-1)
            {
                cout<<i<<" "<<j<<endl;
                cout<<"some faces are neither boundary nor inner face, please check whether your boundary condition set is correct or not"<<endl;
                exit(1);
            }
        }
    }
    for (auto & Boundarie : Boundaries) {
        for (int j = 0; j <Boundarie.cellindex.size() ; ++j) {
            for (int nodeindex : Elements[Boundarie.cellindex[j]].faces[Boundarie.faceindex[j]].vertexes)
            {
                if (!ishave(nodeindex, Boundnodes)) {
                    Boundnodes.push_back(nodeindex);
                    Boundnodesindex.push_back(Boundarie.index);
                }
            }
        }
    }
    Boundnodes_cell.resize(Boundnodes.size());
    for (int i = 0; i < Boundnodes.size() ; ++i) {
        for (int j = 0; j <Boundaries.size() ; ++j) {
            if (ishave(Boundnodes[i],Elements[Boundaries[j].cellindex[0]].faces[Boundaries[j].faceindex[0]].vertexes))
            {
                Boundnodes_cell[i].push_back(j);
            }
        }
    }
    for (auto & Boundarie : Boundaries) {
        Boundarie.connection=-1;
        if (bcs->boundaryConditions[Boundarie.index].type<0)
        {
            double distance=1000;
            for (int j=0; j<Boundaries.size() ;j++) {
                if (bcs->boundaryConditions[Boundarie.index].index==-bcs->boundaryConditions[Boundaries[j].index].type-1)
                {
                    double distance1=get_distance( Elements[Boundarie.cellindex[0]].faces[Boundarie.faceindex[0]].center,Elements[Boundaries[j].cellindex[0]].faces[Boundaries[j].faceindex[0]].center);
                    if (distance1<distance)
                    {
                        distance=distance1;
                        Boundarie.connection=j;
                    }
                }
            }
        }
    }
    Boundnodesconnect.resize(Boundnodes.size());
    for (int i=0;i<Boundnodes.size();i++) {
        Boundnodesconnect[i]=-1;
        if (bcs->boundaryConditions[Boundnodesindex[i]].type<0)
        {
            double distance=100;
            for (int j=0; j<Boundnodes.size() ;j++) {
                if (bcs->boundaryConditions[Boundnodesindex[i]].index==-bcs->boundaryConditions[Boundnodesindex[j]].type-1)
                {
                    double distance1=Nodes[Boundnodes[i]].distance(Nodes[Boundnodes[j]]);
                    if (distance1<distance)
                    {
                        distance=distance1;
                        Boundnodesconnect[i]=j;
                    }
                }

            }
        }
    }
    for (auto &Boundarie:Boundaries) {
        Boundarie.type=bcs->boundaryConditions[Boundarie.index].type;
        Boundarie.Temperature=bcs->boundaryConditions[Boundarie.index].temperature;
        Boundarie.trans.resize(Boundarie.cellindex.size());
        for (int i = 0; i <Boundarie.cellindex.size() ; ++i) {
            for (int j = 0; j < bcs->boundaryConditions[Boundarie.index].matter_trans.size(); ++j) {
                if (Elements[Boundarie.cellindex[i]].matter==bcs->boundaryConditions[Boundarie.index].matter_trans[j].matter)
                {
                    Boundarie.trans[i]=bcs->boundaryConditions[Boundarie.index].matter_trans[j].trans;
                }
            }

        }
    }
    for (int i = 0; i <Boundaries.size() ; ++i) {
        double ie1=Boundaries[i].cellindex[0];
        double iface1=Boundaries[i].faceindex[0];
        for (int j = 0; j <Elements[ie1].faces[iface1].vertexes.size() ; ++j) {
            for (int k = 0; k <Nodes[Elements[ie1].faces[iface1].vertexes[j]].cells.size() ; ++k) {
                int cellindex=Nodes[Elements[ie1].faces[iface1].vertexes[j]].cells[k];
                for (int l = 0; l <Elements[cellindex].faces.size() ; ++l) {
                    if (Elements[cellindex].faces[l].bound>=0)
                    {
                        if (Boundaries[Elements[cellindex].faces[l].bound].index==Boundaries[i].index)
                        {
                            if(!ishave(Elements[cellindex].faces[l].bound,Boundaries[i].neighbors))
                            Boundaries[i].neighbors.push_back(Elements[cellindex].faces[l].bound);
                        }
                    }
                }

            }
        }
    }
    vector<vector<int>> isinfaces(this->Elements.size());
    for (int i = 0; i <Elements.size()  ; ++i) {
        isinfaces[i].resize(Elements[i].faces.size(),1);
    }
    for (int i = 0; i <Elements.size()  ; ++i) {
        for (int j = 0; j <Elements[i].faces.size() ; ++j) {
            if (isinfaces[i][j]==1)
            {
                faces F1;
                F1.cellindex.push_back(i);
                F1.faceindex.push_back(j);
                if (Elements[i].faces[j].bound==-1)
                {
                    F1.cellindex.push_back(Elements[i].faces[j].neighbor);
                    for (int k = 0; k <Elements[Elements[i].faces[j].neighbor].faces.size() ; ++k) {
                        if (Elements[Elements[i].faces[j].neighbor].faces[k].neighbor==i)
                        {
                            F1.faceindex.push_back(k);
                            isinfaces[Elements[i].faces[j].neighbor][k]=0;
                        }
                    }

                }
                F1.Boundindex=Elements[i].faces[j].bound;

                Element_Faces.push_back(F1);
            }
        }
    }
    //cout<<Element_Faces.size()<<endl;

}

BTEMesh::BTEMesh(int Dimension_Geometry,double L_x,double L_y,double L_z,std::vector<double> &nodeX, std::vector<double> &nodeY, std::vector<double> &nodeZ,
                 std::vector<std::vector<int>> &volumeElements, std::vector<int> &volumeElementIndex,
                 std::vector<std::vector<int>> &boundaryElements, std::vector<int> &boundaryElementIndex) {
    this->Dimension = Dimension_Geometry;
    this->L_x=L_x;
    this->L_y=L_y;
    this->L_z=L_z;
    Nodes.resize(nodeX.size());

    for (int i = 0; i < nodeX.size(); ++i) {
        Nodes[i].x=nodeX[i];
        Nodes[i].y=nodeY[i];
        Nodes[i].z=nodeZ[i];
    }
    Boundaries.resize(boundaryElements.size());

    Boundaryindex.resize(boundaryElements.size());

    for (int i = 0; i < boundaryElements.size(); ++i) {
        Boundaryindex[i]=boundaryElements[i];
        Boundary bound;
        bound.index=boundaryElementIndex[i];
        Boundaries[i]=bound;
    }

    Elements.resize(volumeElements.size());
    for (int i = 0; i < Elements.size(); ++i) {
        Elements[i].vertexes=volumeElements[i];
        Elements[i].index=volumeElementIndex[i];
    }

}

void BTEMesh::setMeshParams(BTEBand *bands)
{
    // if (Dimension_Geometry == 1)
    if (Dimension == 1)
    {
        for (int i = 0; i < Elements.size(); ++i)
        {
            if(Elements[i].index>=bands->geo_matter_index.size())
            {
                cout<<"Error: some region does not have material, check PHONON and GEOMETRY"<<endl;
                exit(0);
            }
            Elements[i].matter = bands->geo_matter_index[Elements[i].index];
        }
        // get_volume
        for (int i = 0; i < Elements.size(); ++i)
        {
            Elements[i].center.x = 1.0 / 2 * (Nodes[Elements[i].vertexes[0]].x + Nodes[Elements[i].vertexes[1]].x);
            Elements[i].faces.resize(2);
            Elements[i].volume = Nodes[Elements[i].vertexes[1]].x-Nodes[Elements[i].vertexes[0]].x;
            Elements[i].faces[0].norm.x = -1;
            Elements[i].faces[1].norm.x = 1;
            Elements[i].faces[0].area = 1;
            Elements[i].faces[1].area = 1;
            Elements[i].faces[0].vertexes.push_back(Elements[i].vertexes[0]);
            Elements[i].faces[1].vertexes.push_back(Elements[i].vertexes[1]);
            Elements[i].faces[0].center.x = Nodes[Elements[i].vertexes[0]].x;
            Elements[i].faces[1].center.x = Nodes[Elements[i].vertexes[1]].x;
            Elements[i].faces[0].direction.y = 1;
            Elements[i].faces[1].direction.y = -1;
        }


        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].vertexes.size(); ++j)
            {
                Nodes[Elements[i].vertexes[j]].cells.push_back(i);
            }
        }

        // get_neighbor
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].faces[j].neighbor = -1;
                for (int k = 0; k < Nodes[Elements[i].faces[j].vertexes[0]].cells.size(); ++k)
                {
                    int index = Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index != i)
                    {
                        if (ishave(Elements[i].faces[j].vertexes[0], Elements[index].vertexes) == 1)
                        {
                            Elements[i].faces[j].neighbor = index;
                        }
                    }
                }
                Elements[i].faces[j].bound = -1;
            }
        }
        for (int i = 0; i < Boundaryindex.size(); ++i)
        {
            for (int j = 0; j < Nodes[Boundaryindex[i][0]].cells.size(); ++j)
            {
                int index = Nodes[Boundaryindex[i][0]].cells[j];
                for (int k = 0; k < Elements[index].faces.size(); ++k)
                {
                    if (Boundaryindex[i][0] == Elements[index].faces[k].vertexes[0] )
                    {
                            Boundaries[i].cellindex.push_back(index);
                            Boundaries[i].faceindex.push_back(k);
                    }
                }
            }
            Boundaries[i].macro_Temp.resize(Boundaries[i].cellindex.size());
            Boundaries[i].emit_temp.resize(Boundaries[i].cellindex.size());
            for (int j = 0; j < Boundaries[i].cellindex.size(); ++j)
            {
                Boundaries[i].macro_Temp[j].matter = Elements[Boundaries[i].cellindex[j]].matter;
                Boundaries[i].macro_Temp[j].Temp = 0;
            }
        }
    }
    else if (Dimension == 2)
    {
        for (int i = 0; i < Elements.size(); ++i)
        {
            if(Elements[i].index>=bands->geo_matter_index.size())
            {
                cout<<"Error: some region does not have material, check PHONON and GEOMETRY"<<endl;
                exit(0);
            }
            Elements[i].matter = bands->geo_matter_index[Elements[i].index];
        }

        // get_volume
        for (int i = 0; i < Elements.size(); ++i)
        {
            if (Elements[i].vertexes.size() == 3)
            {
                vector<Point> tri;
                tri.push_back(Nodes[Elements[i].vertexes[0]]);
                tri.push_back(Nodes[Elements[i].vertexes[1]]);
                tri.push_back(Nodes[Elements[i].vertexes[2]]);
                Elements[i].volume = get_tri_area(tri);
                Elements[i].center = ((Nodes[Elements[i].vertexes[0]] + Nodes[Elements[i].vertexes[1]]) + Nodes[Elements[i].vertexes[2]]) / 3;
            }
            else if (Elements[i].vertexes.size() == 4)
            {
                vector<Point> tri(3);
                tri[0] = Nodes[Elements[i].vertexes[0]];
                tri[1] = Nodes[Elements[i].vertexes[1]];
                tri[2] = Nodes[Elements[i].vertexes[2]];
                Elements[i].volume = get_tri_area(tri);
                tri[0] = Nodes[Elements[i].vertexes[1]];
                tri[1] = Nodes[Elements[i].vertexes[2]];
                tri[2] = Nodes[Elements[i].vertexes[3]];
                Elements[i].volume += get_tri_area(tri);
                Elements[i].center = ((Nodes[Elements[i].vertexes[0]] + Nodes[Elements[i].vertexes[1]]) + Nodes[Elements[i].vertexes[2]] + Nodes[Elements[i].vertexes[3]]) / 4;
            }
        }
        // get_faces
        for (int i = 0; i < Elements.size(); ++i)
        {
            if (Elements[i].vertexes.size() == 4)
            {
                Elements[i].faces.resize(4);
                for (int j = 0; j < Elements[i].faces.size(); ++j)
                {
                    Elements[i].faces[j].vertexes.resize(2);
                }
                Elements[i].faces[0].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1] = Elements[i].vertexes[1];

                Elements[i].faces[1].vertexes[0] = Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1] = Elements[i].vertexes[3];

                Elements[i].faces[2].vertexes[0] = Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[1] = Elements[i].vertexes[2];

                Elements[i].faces[3].vertexes[0] = Elements[i].vertexes[2];
                Elements[i].faces[3].vertexes[1] = Elements[i].vertexes[0];

                for (int j = 0; j < 4; ++j)
                {
                    Elements[i].faces[j].area = get_distance(Nodes[Elements[i].faces[j].vertexes[0]], Nodes[Elements[i].faces[j].vertexes[1]]);
                    Elements[i].faces[j].center = (Nodes[Elements[i].faces[j].vertexes[0]] + Nodes[Elements[i].faces[j].vertexes[1]]) / 2;
                    vec a = (Nodes[Elements[i].faces[j].vertexes[0]] - Nodes[Elements[i].faces[j].vertexes[1]]) / Elements[i].faces[j].area;
                    Elements[i].faces[j].norm.x = a.y;
                    Elements[i].faces[j].norm.y = -a.x;
                    Elements[i].faces[j].norm.z = 0;
                    vec a1 = Elements[i].center - Elements[i].faces[j].center;
                    if (a1 * Elements[i].faces[j].norm > 0)
                    {
                        Elements[i].faces[j].norm.x = -Elements[i].faces[j].norm.x;
                        Elements[i].faces[j].norm.y = -Elements[i].faces[j].norm.y;
                    }
                }
            }
            else if (Elements[i].vertexes.size() == 3)
            {
                Elements[i].faces.resize(3);
                for (int j = 0; j < Elements[i].faces.size(); ++j)
                {
                    Elements[i].faces[j].vertexes.resize(2);
                    Elements[i].faces[j].vertexes[0] = Elements[i].vertexes[(j + 1) % 3];
                    Elements[i].faces[j].vertexes[1] = Elements[i].vertexes[(j + 2) % 3];
                    Elements[i].faces[j].area = get_distance(Nodes[Elements[i].faces[j].vertexes[0]], Nodes[Elements[i].faces[j].vertexes[1]]);
                    Elements[i].faces[j].center = (Nodes[Elements[i].faces[j].vertexes[0]] + Nodes[Elements[i].faces[j].vertexes[1]]) / 2;
                    vec a = (Nodes[Elements[i].faces[j].vertexes[0]] - Nodes[Elements[i].faces[j].vertexes[1]]) / Elements[i].faces[j].area;
                    Elements[i].faces[j].norm.x = a.y;
                    Elements[i].faces[j].norm.y = -a.x;
                    Elements[i].faces[j].norm.z = 0;
                    vec a1 = Elements[i].center - Elements[i].faces[j].center;
                    if (a1 * Elements[i].faces[j].norm > 0)
                    {
                        Elements[i].faces[j].norm.x = -Elements[i].faces[j].norm.x;
                        Elements[i].faces[j].norm.y = -Elements[i].faces[j].norm.y;
                    }
                }
            }
        }
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].vertexes.size(); ++j)
            {
                Nodes[Elements[i].vertexes[j]].cells.push_back(i);
            }
        }

        // get_neighbor
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].faces[j].neighbor = -1;
                for (int k = 0; k < Nodes[Elements[i].faces[j].vertexes[0]].cells.size(); ++k)
                {
                    int index = Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index != i)
                    {
                        if (ishave(Elements[i].faces[j].vertexes[0], Elements[index].vertexes) == 1 &&
                            ishave(Elements[i].faces[j].vertexes[1], Elements[index].vertexes) == 1)
                        {
                            Elements[i].faces[j].neighbor = index;
                        }
                    }
                }
                Elements[i].faces[j].bound = -1;
            }
        }
        for (int i = 0; i < Boundaryindex.size(); ++i)
        {
            for (int j = 0; j < Nodes[Boundaryindex[i][0]].cells.size(); ++j)
            {
                int index = Nodes[Boundaryindex[i][0]].cells[j];
                for (int k = 0; k < Elements[index].faces.size(); ++k)
                {
                    if (Boundaryindex[i][1] == Elements[index].faces[k].vertexes[0] ||
                        Boundaryindex[i][1] == Elements[index].faces[k].vertexes[1])
                    {
                        if (Boundaryindex[i][0] == Elements[index].faces[k].vertexes[0] ||
                            Boundaryindex[i][0] == Elements[index].faces[k].vertexes[1])
                        {
                            Boundaries[i].cellindex.push_back(index);
                            Boundaries[i].faceindex.push_back(k);
                        }
                    }
                }
            }
            Boundaries[i].macro_Temp.resize(Boundaries[i].cellindex.size());
            Boundaries[i].emit_temp.resize(Boundaries[i].cellindex.size());
            for (int j = 0; j < Boundaries[i].cellindex.size(); ++j)
            {
                Boundaries[i].macro_Temp[j].matter = Elements[Boundaries[i].cellindex[j]].matter;
                Boundaries[i].macro_Temp[j].Temp = 0;
            }
        }
    }
    else if (Dimension == 3)
    {
        for (int i = 0; i < Elements.size(); ++i)
        {

            Elements[i].matter = bands->geo_matter_index[Elements[i].index];
        }
        // get_volume

        for (int i = 0; i < Elements.size(); ++i)
        {
            Elements[i].center.x = 0;
            Elements[i].center.y = 0;
            Elements[i].center.z = 0;
            for (int j = 0; j < Elements[i].vertexes.size(); j++)
            {
                Elements[i].center = Elements[i].center + Nodes[Elements[i].vertexes[j]] / Elements[i].vertexes.size();
            }
            if (Elements[i].vertexes.size() == 4)
            {
                vector<Point> a;
                a.resize(4);
                for (int j = 0; j < Elements[i].vertexes.size(); ++j)
                {
                    a[j] = Nodes[Elements[i].vertexes[j]];
                }
                Elements[i].volume = get_tet_volume(a);
            }
            else if (Elements[i].vertexes.size() == 5)
            {
                Elements[i].volume = 0;
                vector<Point> a;
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[0]];
                a[1] = Nodes[Elements[i].vertexes[1]];
                a[2] = Nodes[Elements[i].vertexes[2]];
                a[3] = Nodes[Elements[i].vertexes[4]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[1]];
                a[1] = Nodes[Elements[i].vertexes[2]];
                a[2] = Nodes[Elements[i].vertexes[3]];
                a[3] = Nodes[Elements[i].vertexes[4]];
                Elements[i].volume += get_tet_volume(a);
            }
            else if (Elements[i].vertexes.size() == 6)
            {
                Elements[i].volume = 0;
                vector<Point> a;
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[0]];
                a[1] = Nodes[Elements[i].vertexes[1]];
                a[2] = Nodes[Elements[i].vertexes[2]];
                a[3] = Nodes[Elements[i].vertexes[3]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[1]];
                a[1] = Nodes[Elements[i].vertexes[2]];
                a[2] = Nodes[Elements[i].vertexes[3]];
                a[3] = Nodes[Elements[i].vertexes[4]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[2]];
                a[1] = Nodes[Elements[i].vertexes[3]];
                a[2] = Nodes[Elements[i].vertexes[4]];
                a[3] = Nodes[Elements[i].vertexes[5]];
                Elements[i].volume += get_tet_volume(a);
            }
            else if (Elements[i].vertexes.size() == 8)
            {
                Elements[i].volume = 0;
                vector<Point> a;
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[0]];
                a[1] = Nodes[Elements[i].vertexes[1]];
                a[2] = Nodes[Elements[i].vertexes[2]];
                a[3] = Nodes[Elements[i].vertexes[4]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[1]];
                a[1] = Nodes[Elements[i].vertexes[2]];
                a[2] = Nodes[Elements[i].vertexes[3]];
                a[3] = Nodes[Elements[i].vertexes[7]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[2]];
                a[1] = Nodes[Elements[i].vertexes[4]];
                a[2] = Nodes[Elements[i].vertexes[6]];
                a[3] = Nodes[Elements[i].vertexes[7]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[1]];
                a[1] = Nodes[Elements[i].vertexes[2]];
                a[2] = Nodes[Elements[i].vertexes[4]];
                a[3] = Nodes[Elements[i].vertexes[7]];
                Elements[i].volume += get_tet_volume(a);
                a.resize(4);
                a[0] = Nodes[Elements[i].vertexes[1]];
                a[1] = Nodes[Elements[i].vertexes[4]];
                a[2] = Nodes[Elements[i].vertexes[5]];
                a[3] = Nodes[Elements[i].vertexes[7]];
                Elements[i].volume += get_tet_volume(a);
            }
        }
        // set_face
        for (int i = 0; i < Elements.size(); ++i)
        {
            if (Elements[i].vertexes.size() == 4)
            {
                Elements[i].faces.resize(4);
                for (int j = 0; j < 4; ++j)
                {
                    Elements[i].faces[j].vertexes.resize(3);
                }
                Elements[i].faces[0].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2] = Elements[i].vertexes[2];

                Elements[i].faces[1].vertexes[0] = Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1] = Elements[i].vertexes[2];
                Elements[i].faces[1].vertexes[2] = Elements[i].vertexes[3];

                Elements[i].faces[2].vertexes[0] = Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[1] = Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[2] = Elements[i].vertexes[0];

                Elements[i].faces[3].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[2] = Elements[i].vertexes[3];
            }
            else if (Elements[i].vertexes.size() == 5)
            {
                Elements[i].faces.resize(5);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2] = Elements[i].vertexes[2];
                Elements[i].faces[0].vertexes[3] = Elements[i].vertexes[3];

                Elements[i].faces[1].vertexes.resize(3);
                Elements[i].faces[1].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[1].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[2] = Elements[i].vertexes[4];

                Elements[i].faces[2].vertexes.resize(3);
                Elements[i].faces[2].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1] = Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[2] = Elements[i].vertexes[4];

                Elements[i].faces[3].vertexes.resize(3);
                Elements[i].faces[3].vertexes[0] = Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[1] = Elements[i].vertexes[3];
                Elements[i].faces[3].vertexes[2] = Elements[i].vertexes[4];

                Elements[i].faces[4].vertexes.resize(3);
                Elements[i].faces[4].vertexes[0] = Elements[i].vertexes[2];
                Elements[i].faces[4].vertexes[1] = Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[2] = Elements[i].vertexes[4];
            }
            else if (Elements[i].vertexes.size() == 6)
            {
                Elements[i].faces.resize(5);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2] = Elements[i].vertexes[3];
                Elements[i].faces[0].vertexes[3] = Elements[i].vertexes[4];

                Elements[i].faces[1].vertexes.resize(4);
                Elements[i].faces[1].vertexes[0] = Elements[i].vertexes[1];
                Elements[i].faces[1].vertexes[1] = Elements[i].vertexes[2];
                Elements[i].faces[1].vertexes[2] = Elements[i].vertexes[4];
                Elements[i].faces[1].vertexes[3] = Elements[i].vertexes[5];

                Elements[i].faces[2].vertexes.resize(4);
                Elements[i].faces[2].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1] = Elements[i].vertexes[2];
                Elements[i].faces[2].vertexes[2] = Elements[i].vertexes[3];
                Elements[i].faces[2].vertexes[3] = Elements[i].vertexes[5];

                Elements[i].faces[3].vertexes.resize(3);
                Elements[i].faces[3].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[3].vertexes[2] = Elements[i].vertexes[2];

                Elements[i].faces[4].vertexes.resize(3);
                Elements[i].faces[4].vertexes[0] = Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[1] = Elements[i].vertexes[4];
                Elements[i].faces[4].vertexes[2] = Elements[i].vertexes[5];
            }
            else if (Elements[i].vertexes.size() == 8)
            {
                Elements[i].faces.resize(6);

                Elements[i].faces[0].vertexes.resize(4);
                Elements[i].faces[0].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[0].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[0].vertexes[2] = Elements[i].vertexes[2];
                Elements[i].faces[0].vertexes[3] = Elements[i].vertexes[3];

                Elements[i].faces[1].vertexes.resize(4);
                Elements[i].faces[1].vertexes[0] = Elements[i].vertexes[4];
                Elements[i].faces[1].vertexes[1] = Elements[i].vertexes[5];
                Elements[i].faces[1].vertexes[2] = Elements[i].vertexes[6];
                Elements[i].faces[1].vertexes[3] = Elements[i].vertexes[7];

                Elements[i].faces[2].vertexes.resize(4);
                Elements[i].faces[2].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[2].vertexes[1] = Elements[i].vertexes[1];
                Elements[i].faces[2].vertexes[2] = Elements[i].vertexes[4];
                Elements[i].faces[2].vertexes[3] = Elements[i].vertexes[5];

                Elements[i].faces[3].vertexes.resize(4);
                Elements[i].faces[3].vertexes[0] = Elements[i].vertexes[0];
                Elements[i].faces[3].vertexes[1] = Elements[i].vertexes[2];
                Elements[i].faces[3].vertexes[2] = Elements[i].vertexes[4];
                Elements[i].faces[3].vertexes[3] = Elements[i].vertexes[6];

                Elements[i].faces[4].vertexes.resize(4);
                Elements[i].faces[4].vertexes[0] = Elements[i].vertexes[1];
                Elements[i].faces[4].vertexes[1] = Elements[i].vertexes[3];
                Elements[i].faces[4].vertexes[2] = Elements[i].vertexes[5];
                Elements[i].faces[4].vertexes[3] = Elements[i].vertexes[7];

                Elements[i].faces[5].vertexes.resize(4);
                Elements[i].faces[5].vertexes[0] = Elements[i].vertexes[2];
                Elements[i].faces[5].vertexes[1] = Elements[i].vertexes[3];
                Elements[i].faces[5].vertexes[2] = Elements[i].vertexes[6];
                Elements[i].faces[5].vertexes[3] = Elements[i].vertexes[7];
            }
        }
        // set_face_area_center
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].faces[j].center.x = 0;
                Elements[i].faces[j].center.y = 0;
                Elements[i].faces[j].center.z = 0;
                for (int k = 0; k < Elements[i].faces[j].vertexes.size(); ++k)
                {
                    Elements[i].faces[j].center = Elements[i].faces[j].center + Nodes[Elements[i].faces[j].vertexes[k]] / Elements[i].faces[j].vertexes.size();
                }
                if (Elements[i].faces[j].vertexes.size() == 4)
                {
                    Elements[i].faces[j].area = 0;
                    vector<Point> tri;
                    tri.resize(3);
                    tri[0] = Nodes[Elements[i].faces[j].vertexes[0]];
                    tri[1] = Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[2] = Nodes[Elements[i].faces[j].vertexes[2]];
                    Elements[i].faces[j].area += get_tri_area(tri);
                    // cout<<get_tri_area(tri)<<endl;
                    tri.resize(3);
                    tri[0] = Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[1] = Nodes[Elements[i].faces[j].vertexes[2]];
                    tri[2] = Nodes[Elements[i].faces[j].vertexes[3]];
                    Elements[i].faces[j].area += get_tri_area(tri);
                    // cout<<get_tri_area(tri)<<endl;
                }
                else
                {
                    vector<Point> tri;
                    tri.resize(3);
                    tri[0] = Nodes[Elements[i].faces[j].vertexes[0]];
                    tri[1] = Nodes[Elements[i].faces[j].vertexes[1]];
                    tri[2] = Nodes[Elements[i].faces[j].vertexes[2]];
                    Elements[i].faces[j].area = get_tri_area(tri);
                }
            }
        }

        // get_face_norm
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                vec a, b;
                a = Nodes[Elements[i].faces[j].vertexes[0]] - Nodes[Elements[i].faces[j].vertexes[1]];
                b = Nodes[Elements[i].faces[j].vertexes[0]] - Nodes[Elements[i].faces[j].vertexes[2]];
                Elements[i].faces[j].norm = cross_product(a, b);
                double length = get_length(Elements[i].faces[j].norm);
                Elements[i].faces[j].norm = Elements[i].faces[j].norm / length;
                a = Elements[i].faces[j].center - Elements[i].center;
                if (dot_product(a, Elements[i].faces[j].norm) < 0)
                {
                    Elements[i].faces[j].norm.x = -Elements[i].faces[j].norm.x;
                    Elements[i].faces[j].norm.y = -Elements[i].faces[j].norm.y;
                    Elements[i].faces[j].norm.z = -Elements[i].faces[j].norm.z;
                }
            }
        }
        //
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].vertexes.size(); ++j)
            {
                Nodes[Elements[i].vertexes[j]].cells.push_back(i);
            }
        }
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].faces[j].neighbor = -1;
                Elements[i].faces[j].bound = -1;
                for (int k = 0; k < Nodes[Elements[i].faces[j].vertexes[0]].cells.size(); ++k)
                {
                    int index = Nodes[Elements[i].faces[j].vertexes[0]].cells[k];
                    if (index != i)
                    {
                        vector<int> vertexes = Elements[index].vertexes;
                        if (ishave(Elements[i].faces[j].vertexes[1], vertexes) == 1 && ishave(Elements[i].faces[j].vertexes[2], vertexes) == 1 && ishave(Elements[i].faces[j].vertexes[0], vertexes) == 1)
                        {
                            Elements[i].faces[j].neighbor = index;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < Boundaryindex.size(); ++i)
        {
            for (int j = 0; j < Nodes[Boundaryindex[i][0]].cells.size(); ++j)
            {
                int index = Nodes[Boundaryindex[i][0]].cells[j];
                for (int k = 0; k < Elements[index].faces.size(); ++k)
                {
                    if (ishave(Boundaryindex[i][0], Elements[index].faces[k].vertexes) == 1 && ishave(Boundaryindex[i][1], Elements[index].faces[k].vertexes) == 1 && ishave(Boundaryindex[i][2], Elements[index].faces[k].vertexes) == 1)
                    {
                        Boundaries[i].cellindex.push_back(index);
                        Boundaries[i].faceindex.push_back(k);
                    }
                }
            }
            Boundaries[i].macro_Temp.resize(Boundaries[i].cellindex.size());
            Boundaries[i].emit_temp.resize(Boundaries[i].cellindex.size());
            for (int j = 0; j < Boundaries[i].cellindex.size(); ++j)
            {
                Boundaries[i].macro_Temp[j].matter = Elements[Boundaries[i].cellindex[j]].matter;
                Boundaries[i].macro_Temp[j].Temp = 0;
            }
        }
    }
}
void BTEMesh::setMeshParams1(BTEBoundaryCondition *bcs)
{
    vector<Boundary> Boundaries1 = Boundaries;
    vector<Boundary> Boundaries2 =Boundaries;
    Boundaries.resize(0);
    for (int i = 0; i < Boundaries1.size(); ++i)
    {
        for (int j = 0; j < bcs->boundaryConditions.size(); ++j)
        {
            if (Boundaries1[i].index == bcs->boundaryConditions[j].index)
            {
                Boundaries2[i].index = j;
                Boundaries.push_back(Boundaries2[i]);
            }
        }
    }
    for (int i = 0; i < Boundaries.size(); ++i)
    {
        for (int j = 0; j < Boundaries[i].cellindex.size(); ++j)
        {
            Elements[Boundaries[i].cellindex[j]].faces[Boundaries[i].faceindex[j]].bound = i;
        }
    }
    for (int i = 0; i < Elements.size(); ++i)
    {
        for (int j = 0; j < Elements[i].faces.size(); ++j)
        {
            if (Elements[i].faces[j].bound == -1 && Elements[i].faces[j].neighbor == -1)
            {
                cout << i << " " << j << endl;
                cout<<"some faces are neither boundary nor inner face, please check whether your boundary condition set is correct or not"<<endl;
                exit(1);
            }
        }
    }
    for (auto &Boundarie : Boundaries)
    {
        for (int j = 0; j < Boundarie.cellindex.size(); ++j)
        {
            for (int nodeindex : Elements[Boundarie.cellindex[j]].faces[Boundarie.faceindex[j]].vertexes)
            {
                if (!ishave(nodeindex, Boundnodes))
                {
                    Boundnodes.push_back(nodeindex);
                    Boundnodesindex.push_back(Boundarie.index);
                }
            }
        }
    }
    Boundnodes_cell.resize(Boundnodes.size());
    for (int i = 0; i < Boundnodes.size(); ++i)
    {
        for (int j = 0; j < Boundaries.size(); ++j)
        {
            if (ishave(Boundnodes[i], Elements[Boundaries[j].cellindex[0]].faces[Boundaries[j].faceindex[0]].vertexes))
            {
                Boundnodes_cell[i].push_back(j);
            }
        }
    }
    for (auto &Boundarie : Boundaries)
    {
        Boundarie.connection = -1;
        if (bcs->boundaryConditions[Boundarie.index].type < 0)
        {
            double distance = 1000;
            for (int j = 0; j < Boundaries.size(); j++)
            {
                if (bcs->boundaryConditions[Boundarie.index].index == -bcs->boundaryConditions[Boundaries[j].index].type - 1)
                {
                    double distance1 = get_distance(Elements[Boundarie.cellindex[0]].faces[Boundarie.faceindex[0]].center, Elements[Boundaries[j].cellindex[0]].faces[Boundaries[j].faceindex[0]].center);
                    if (distance1 < distance)
                    {
                        distance = distance1;
                        Boundarie.connection = j;
                    }
                }
            }
        }
    }
    Boundnodesconnect.resize(Boundnodes.size());
    for (int i = 0; i < Boundnodes.size(); i++)
    {
        Boundnodesconnect[i] = -1;
        if (bcs->boundaryConditions[Boundnodesindex[i]].type < 0)
        {
            double distance = 100;
            for (int j = 0; j < Boundnodes.size(); j++)
            {
                if (bcs->boundaryConditions[Boundnodesindex[i]].index == -bcs->boundaryConditions[Boundnodesindex[j]].type - 1)
                {
                    double distance1 = Nodes[Boundnodes[i]].distance(Nodes[Boundnodes[j]]);
                    if (distance1 < distance)
                    {
                        distance = distance1;
                        Boundnodesconnect[i] = j;
                    }
                }
            }
        }
    }
    for (auto &Boundarie : Boundaries)
    {
        Boundarie.type = bcs->boundaryConditions[Boundarie.index].type;
        Boundarie.Temperature = bcs->boundaryConditions[Boundarie.index].temperature;
        Boundarie.trans.resize(Boundarie.cellindex.size());
        for (int i = 0; i < Boundarie.cellindex.size(); ++i)
        {
            for (int j = 0; j < bcs->boundaryConditions[Boundarie.index].matter_trans.size(); ++j)
            {
                if (Elements[Boundarie.cellindex[i]].matter == bcs->boundaryConditions[Boundarie.index].matter_trans[j].matter)
                {
                    Boundarie.trans[i] = bcs->boundaryConditions[Boundarie.index].matter_trans[j].trans;
                }
            }
        }
    }
    for (int i = 0; i < Boundaries.size(); ++i)
    {
        double ie1 = Boundaries[i].cellindex[0];
        double iface1 = Boundaries[i].faceindex[0];
        for (int j = 0; j < Elements[ie1].faces[iface1].vertexes.size(); ++j)
        {
            for (int k = 0; k < Nodes[Elements[ie1].faces[iface1].vertexes[j]].cells.size(); ++k)
            {
                int cellindex = Nodes[Elements[ie1].faces[iface1].vertexes[j]].cells[k];
                for (int l = 0; l < Elements[cellindex].faces.size(); ++l)
                {
                    if (Elements[cellindex].faces[l].bound >= 0)
                    {
                        if (Boundaries[Elements[cellindex].faces[l].bound].index == Boundaries[i].index)
                        {
                            if(!ishave(Elements[cellindex].faces[l].bound,Boundaries[i].neighbors)&&Elements[cellindex].faces[l].bound!=i)
                            Boundaries[i].neighbors.push_back(Elements[cellindex].faces[l].bound);
                        }
                    }
                }
            }
        }
    }
    vector<vector<int>> isinfaces(this->Elements.size());
    for (int i = 0; i < Elements.size(); ++i)
    {
        isinfaces[i].resize(Elements[i].faces.size(), 1);
    }
    for (int i = 0; i < Elements.size(); ++i)
    {
        for (int j = 0; j < Elements[i].faces.size(); ++j)
        {
            if (isinfaces[i][j] == 1)
            {
                faces F1;
                F1.cellindex.push_back(i);
                F1.faceindex.push_back(j);
                Elements[i].faces[j].index=Element_Faces.size();
                Elements[i].faces[j].index=Element_Faces.size();
                if (Elements[i].faces[j].bound == -1)
                {
                    F1.cellindex.push_back(Elements[i].faces[j].neighbor);
                    for (int k = 0; k < Elements[Elements[i].faces[j].neighbor].faces.size(); ++k)
                    {
                        if (Elements[Elements[i].faces[j].neighbor].faces[k].neighbor == i)
                        {
                            F1.faceindex.push_back(k);
                            isinfaces[Elements[i].faces[j].neighbor][k] = 0;
                            Elements[Elements[i].faces[j].neighbor].faces[k].index=Element_Faces.size();
                            Elements[Elements[i].faces[j].neighbor].faces[k].index=Element_Faces.size();
                        }
                    }
                }

                F1.Boundindex = Elements[i].faces[j].bound;

                Element_Faces.push_back(F1);
            }
        }
    }
    // cout<<Element_Faces.size()<<endl;
}
void BTEMesh::BTEMesh_initialTemp(ifstream &initialTemp,double Tref)//jiaxuan
{
    if (!initialTemp.is_open())
    {
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].initial_temperature = Tref;
            }
        }
    }
    else
    {
        cout<<"**************************"<<endl;
        cout<<"Begin to read initial temperature file !"<<endl;
        if (Dimension==1) {
            string line;
            string strmesh;
            int numtempnode=0;
            while(getline(initialTemp, strmesh)){  //
                if(strmesh.find("Nodes:")>0  //
                   && strmesh.find("Nodes:")<strmesh.length()){ //
                    int flag = 0;
                    for (int i = 0; i<strmesh.length(); i++){
                        if (strmesh[i]>='0' && strmesh[i]<='9'){
                            numtempnode = numtempnode * 10 + strmesh[i] - '0';
                            flag = 1;
                        }
                        if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                            break;
                        }
                    }
                    break;
                }
            }

            cout<<"Number of initial temperature nodes: "<<numtempnode<<endl;
            this->Tempnodes.resize(numtempnode);
            for (int i = 0; i < numtempnode; i++) {
                getline(initialTemp, strmesh);
                stringstream sss;
                sss.str(strmesh);
                string coord_x, coord_temp;
                sss >> coord_x >> coord_temp;
                Tempnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x;
                Tempnodes[i].y = 0;
                Tempnodes[i].z = 0;
                Tempnodes[i].initialtemp = strtod(coord_temp.c_str(), NULL);
            }

            double distances[Tempnodes.size()]={0.0};
            double selectMindistance[2][2]={0.0};
            for (int i=0; i<Nodes.size(); i++){
                int k=0;
                for (int j=0 ;j<Tempnodes.size(); j++){
                    if (Nodes[i].x==Tempnodes[j].x){
                        Nodes[i].initialtemp=Tempnodes[j].initialtemp;
                        k=1;
                        break;
                    }
                }

                if (k==0){
                    for (int j=0 ;j<Tempnodes.size(); j++){
                        distances[j] = get_distance(Tempnodes[j],Nodes[i]);
                    }
                    for (int numNearest = 0; numNearest<2; numNearest++){
                        int minPosition = min_element(distances,distances+Tempnodes.size()) - distances;
                        selectMindistance[numNearest][0] = 1/distances[minPosition];
                        selectMindistance[numNearest][1] = Tempnodes[minPosition].initialtemp;
                        distances[minPosition]=1000;
                    }
                    double totalRevDistance = 0;
                    for (int numNearest = 0; numNearest<2; numNearest++){
                        Nodes[i].initialtemp += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                        totalRevDistance += selectMindistance[numNearest][0];
                    }
                    Nodes[i].initialtemp = Nodes[i].initialtemp/totalRevDistance;
                }
            }

            for (int i=0; i<Elements.size(); i++){
                Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp)/2;
            }
            cout<<"Successfully read from COORDINATE initial temperature file !"<<endl;
        }

        if (Dimension==2) {
            string line;
            string strmesh;
            int numtempnode=0;
            while(getline(initialTemp, strmesh)){  //
                if(strmesh.find("Nodes:")>0  //
                   && strmesh.find("Nodes:")<strmesh.length()){ //
                    int flag = 0;
                    for (int i = 0; i<strmesh.length(); i++){
                        if (strmesh[i]>='0' && strmesh[i]<='9'){
                            numtempnode = numtempnode * 10 + strmesh[i] - '0';
                            flag = 1;
                        }
                        if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                            break;
                        }
                    }
                    break;
                }
            }
            cout<<"Number of initial temperature nodes: "<<numtempnode<<endl;

            this->Tempnodes.resize(numtempnode);
            for (int i = 0; i < numtempnode; i++) {
                getline(initialTemp, strmesh);
                stringstream sss;
                sss.str(strmesh);
                string coord_x, coord_y, coord_temp;
                sss >> coord_x >> coord_y >> coord_temp;
                Tempnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x;
                Tempnodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;
                Tempnodes[i].z = 0;
                Tempnodes[i].initialtemp = strtod(coord_temp.c_str(), NULL);
            }

            double distances[Tempnodes.size()]={0.0};
            double selectMindistance[4][2]={0.0};
            for (int i=0; i<Nodes.size(); i++){
                int k=0;
                for (int j=0 ;j<Tempnodes.size(); j++){
                    if (Nodes[i].x==Tempnodes[j].x && Nodes[i].y==Tempnodes[j].y){
                        Nodes[i].initialtemp=Tempnodes[j].initialtemp;
                        k=1;
                        break;
                    }
                }
                if (k==0){
                    for (int j=0 ;j<Tempnodes.size(); j++){
                        distances[j] = get_distance(Tempnodes[j],Nodes[i]);
                    }
                    for (int numNearest = 0; numNearest<4; numNearest++){
                        int minPosition = min_element(distances,distances+Tempnodes.size()) - distances;
                        selectMindistance[numNearest][0] = 1/distances[minPosition];
                        selectMindistance[numNearest][1] = Tempnodes[minPosition].initialtemp;
                        distances[minPosition]=1000;
                    }
                    // calculate initial temperature of the node
                    double totalRevDistance = 0;
                    for (int numNearest = 0; numNearest<4; numNearest++){
                        Nodes[i].initialtemp += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                        totalRevDistance += selectMindistance[numNearest][0];
                    }
                    Nodes[i].initialtemp = Nodes[i].initialtemp/totalRevDistance;
                }
            }

            //jiaxuan: set Elements[i].initial_temperature
            for (int i=0; i<Elements.size(); i++){
                if (Elements[i].vertexes.size()==3){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp)/3;
                }
                if (Elements[i].vertexes.size()==4){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp + Nodes[Elements[i].vertexes[3]].initialtemp)/4;
                }
                //cout<<Elements[i].center.x<<" "<<Elements[i].center.y<<" "<<Elements[i].initial_temperature<<endl;
            }
            cout<<"Successfully read from COORDINATE initial temperature file !"<<endl;
        }

        if (Dimension==3) {
            string line;
            string strmesh;
            int numtempnode=0;
            while(getline(initialTemp, strmesh)){  //
                if(strmesh.find("Nodes:")>0  //
                   && strmesh.find("Nodes:")<strmesh.length()){ //
                    int flag = 0;
                    for (int i = 0; i<strmesh.length(); i++){
                        if (strmesh[i]>='0' && strmesh[i]<='9'){
                            numtempnode = numtempnode * 10 + strmesh[i] - '0';
                            flag = 1;
                        }
                        if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                            break;
                        }
                    }
                    break;
                }
            }

            cout<<"Number of initial temperature nodes: "<<numtempnode<<endl;
            this->Tempnodes.resize(numtempnode);
            for (int i = 0; i < numtempnode; i++) {
                getline(initialTemp, strmesh);
                stringstream sss;
                sss.str(strmesh);
                string coord_x, coord_y, coord_z, coord_temp;
                sss >> coord_x >> coord_y >> coord_z >> coord_temp;
                Tempnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x;
                Tempnodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;
                Tempnodes[i].z = strtod(coord_z.c_str(), NULL) * L_z;;
                Tempnodes[i].initialtemp = strtod(coord_temp.c_str(), NULL);
            }

            double distances[Tempnodes.size()]={0.0};
            double selectMindistance[8][2]={0.0};
            for (int i=0; i<Nodes.size(); i++){
                int k=0;
                for (int j=0 ;j<Tempnodes.size(); j++){
                    if (Nodes[i].x==Tempnodes[j].x && Nodes[i].y==Tempnodes[j].y && Nodes[i].z==Tempnodes[j].z){
                        Nodes[i].initialtemp=Tempnodes[j].initialtemp;
                        k=1;
                        break;
                    }
                }

                if (k==0){
                    for (int j=0 ;j<Tempnodes.size(); j++){
                        distances[j] = get_distance(Tempnodes[j],Nodes[i]);
                    }
                    for (int numNearest = 0; numNearest<8; numNearest++){
                        int minPosition = min_element(distances,distances+Tempnodes.size()) - distances;
                        selectMindistance[numNearest][0] = 1/distances[minPosition];
                        selectMindistance[numNearest][1] = Tempnodes[minPosition].initialtemp;
                        distances[minPosition]=1000;
                    }
                    double totalRevDistance = 0;
                    for (int numNearest = 0; numNearest<8; numNearest++){
                        Nodes[i].initialtemp += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                        totalRevDistance += selectMindistance[numNearest][0];
                    }
                    Nodes[i].initialtemp = Nodes[i].initialtemp/totalRevDistance;
                }
            }

            for (int i=0; i<Elements.size(); i++){
                if (Elements[i].vertexes.size()==4){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp + Nodes[Elements[i].vertexes[3]].initialtemp)/4;
                }
                if (Elements[i].vertexes.size()==5){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp + Nodes[Elements[i].vertexes[3]].initialtemp + Nodes[Elements[i].vertexes[4]].initialtemp)/5;
                }
                if (Elements[i].vertexes.size()==6){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp + Nodes[Elements[i].vertexes[3]].initialtemp + Nodes[Elements[i].vertexes[4]].initialtemp + Nodes[Elements[i].vertexes[5]].initialtemp)/6;
                }
                if (Elements[i].vertexes.size()==8){
                    Elements[i].initial_temperature=(Nodes[Elements[i].vertexes[0]].initialtemp + Nodes[Elements[i].vertexes[1]].initialtemp + Nodes[Elements[i].vertexes[2]].initialtemp + Nodes[Elements[i].vertexes[3]].initialtemp + Nodes[Elements[i].vertexes[4]].initialtemp + Nodes[Elements[i].vertexes[5]].initialtemp + Nodes[Elements[i].vertexes[6]].initialtemp + Nodes[Elements[i].vertexes[7]].initialtemp)/8;
                }
            }
            cout<<"Successfully read from COORDINATE initial temperature file !"<<endl;
        }

        cout<<"**************************"<<endl;
    }

}
void BTEMesh::BTEMesh_heatin(ifstream &inHeat, double Uniform_Heat, std::string heat_type) //yufei
{
    if (!inHeat.is_open())
    {
        for (int i = 0; i < Elements.size(); ++i)
        {
            for (int j = 0; j < Elements[i].faces.size(); ++j)
            {
                Elements[i].heat_source = Uniform_Heat;
            }
        }
    }
    else
    {
        cout<<"**************************"<<endl;
        cout<<"Begin to read heat source file !"<<endl;

        if (heat_type == "REGION"){
            int num_of_matter;
            inHeat >> num_of_matter;
            vector<int> index(num_of_matter);
            vector<double> heatsource(num_of_matter);
            string str;
            getline(inHeat, str);
            for (int i = 0; i < num_of_matter; ++i)
            {
                inHeat >> index[i] >> heatsource[i];
            }
            for (int i = 0; i < Elements.size(); ++i)
            {
                if (Elements[i].index > num_of_matter)
                {
                    cout << "Wrong heat file !";
                    exit(0);
                }

                Elements[i].heat_source = heatsource[Elements[i].index];
            }
            cout<<"Successfully read from REGION heat file !"<<endl;
        }


        if (heat_type == "COORDINATE"){
                    //yufei: for 1D mesh
            if (Dimension==1) {
                    string line;
                    string strmesh;
                    int numheatnode=0;
                    while(getline(inHeat, strmesh)){  // 
                        if(strmesh.find("Nodes:")>0  // 
                        && strmesh.find("Nodes:")<strmesh.length()){ // 
                            //yufei: extract the number of heatnodes:
                            int flag = 0;
                            for (int i = 0; i<strmesh.length(); i++){
                                if (strmesh[i]>='0' && strmesh[i]<='9'){
                                    numheatnode = numheatnode * 10 + strmesh[i] - '0';
                                    flag = 1;
                                }
                                if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                                    break;
                                }
                            }
                            break;
                        }
                    }

                    cout<<"Number of heat nodes: "<<numheatnode<<endl;
                    this->Heatnodes.resize(numheatnode); 
                    for (int i = 0; i < numheatnode; i++) {
                        getline(inHeat, strmesh);
                        stringstream sss; 
                        sss.str(strmesh);  
                        string coord_x, coord_heat;
                        sss >> coord_x >> coord_heat; 
                        Heatnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x; 
                        Heatnodes[i].y = 0;  
                        Heatnodes[i].z = 0;
                        Heatnodes[i].heat = strtod(coord_heat.c_str(), NULL);  
                        //cout << "Heatnodes" << Heatnodes[i].x << " " <<  Heatnodes[i].y << " " << Heatnodes[i].z << " " << Heatnodes[i].heat << endl;
                    }

                    //yufei: test whether heatnode is consistent with node of mesh
                    double distances[Heatnodes.size()]={0.0};
                    double selectMindistance[2][2]={0.0}; //yufei:this is to store the min distance and heat, for 1D, consider 2 nearest neibour
                    for (int i=0; i<Nodes.size(); i++){    
                        int k=0;
                        for (int j=0 ;j<Heatnodes.size(); j++){
                            if (Nodes[i].x==Heatnodes[j].x){
                                Nodes[i].heat=Heatnodes[j].heat;
                                k=1; 
                                break; 
                            }
                        }
                        //yufei: if cannot find the same mesh, calculate the distance
                        if (k==0){
                            for (int j=0 ;j<Heatnodes.size(); j++){ 
                                distances[j] = get_distance(Heatnodes[j],Nodes[i]);
                                //cout << "distance" << distances[j] << endl;
                            } 
                            // yufei: find the nearest 2 heatnode
                            for (int numNearest = 0; numNearest<2; numNearest++){
                                int minPosition = min_element(distances,distances+Heatnodes.size()) - distances;
                                selectMindistance[numNearest][0] = 1/distances[minPosition]; 
                                selectMindistance[numNearest][1] = Heatnodes[minPosition].heat;
                                distances[minPosition]=1000;
                                //cout << "selectMindinstace" << selectMindistance[numNearest][0] << " " << selectMindistance[numNearest][1] << endl;
                            }
                            // calculate heat of the node
                            double totalRevDistance = 0;
                            for (int numNearest = 0; numNearest<2; numNearest++){
                                Nodes[i].heat += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                                totalRevDistance += selectMindistance[numNearest][0];
                                //cout<< "TotalRevDistance" << totalRevDistance << endl;
                            }
                            Nodes[i].heat = Nodes[i].heat/totalRevDistance;
                        }
                    }

                //yufei: for 1D mesh, vertex = 2
                    for (int i=0; i<Elements.size(); i++){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat)/2;
                    }
                    cout<<"Successfully read from COORDINATE heat file !"<<endl;
                }


                //yufei: for 2D mesh
            if (Dimension==2) {
                    string line;
                    string strmesh;
                    int numheatnode=0;
                    while(getline(inHeat, strmesh)){  // 
                        if(strmesh.find("Nodes:")>0  // 
                        && strmesh.find("Nodes:")<strmesh.length()){ // 
                            //yufei: extract the number of heatnodes:
                            int flag = 0;
                            for (int i = 0; i<strmesh.length(); i++){
                                if (strmesh[i]>='0' && strmesh[i]<='9'){
                                    numheatnode = numheatnode * 10 + strmesh[i] - '0';
                                    flag = 1;
                                }
                                if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                                    break;
                                }
                            }
                            break;
                        }
                    }

                    cout<<"Number of heat nodes: "<<numheatnode<<endl;
                    this->Heatnodes.resize(numheatnode); 
                    for (int i = 0; i < numheatnode; i++) {
                        getline(inHeat, strmesh);
                        stringstream sss; 
                        sss.str(strmesh);  
                        string coord_x, coord_y, coord_heat;
                        sss >> coord_x >> coord_y >> coord_heat; 
                        Heatnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x; 
                        Heatnodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;  
                        Heatnodes[i].z = 0;
                        Heatnodes[i].heat = strtod(coord_heat.c_str(), NULL); 
                    }

                    //yufei: test whether heatnode is consistent with node of mesh
                    double distances[Heatnodes.size()]={0.0};
                    double selectMindistance[4][2]={0.0}; //yufei:this is to store the min distance and heat, for 2D, consider 4 nearest neibour
                    for (int i=0; i<Nodes.size(); i++){    
                        int k=0;
                        for (int j=0 ;j<Heatnodes.size(); j++){
                            if (Nodes[i].x==Heatnodes[j].x && Nodes[i].y==Heatnodes[j].y){
                                Nodes[i].heat=Heatnodes[j].heat;
                                k=1; 
                                break; 
                            }
                        }
                        //yufei: if cannot find the same mesh, calculate the distance
                        if (k==0){
                            for (int j=0 ;j<Heatnodes.size(); j++){ 
                                distances[j] = get_distance(Heatnodes[j],Nodes[i]);
                            } 
                            // yufei: find the nearest 4 heatnode
                            for (int numNearest = 0; numNearest<4; numNearest++){
                                int minPosition = min_element(distances,distances+Heatnodes.size()) - distances;
                                selectMindistance[numNearest][0] = 1/distances[minPosition]; 
                                selectMindistance[numNearest][1] = Heatnodes[minPosition].heat;
                                distances[minPosition]=1000;
                            }
                            // calculate heat of the node
                            double totalRevDistance = 0;
                            for (int numNearest = 0; numNearest<4; numNearest++){
                                Nodes[i].heat += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                                totalRevDistance += selectMindistance[numNearest][0];
                            }
                            Nodes[i].heat = Nodes[i].heat/totalRevDistance;
                        }
                    }

                //yufei: for 2D mesh, vertex = 3 & 4
                    for (int i=0; i<Elements.size(); i++){
                        if (Elements[i].vertexes.size()==3){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat)/3;
                        }
                        if (Elements[i].vertexes.size()==4){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat + Nodes[Elements[i].vertexes[3]].heat)/4;
                        }
                    }
                    cout<<"Successfully read from COORDINATE heat file !"<<endl;
                }


                //yufei: for 3D mesh
                if (Dimension==3) {
                    string line;
                    string strmesh;
                    int numheatnode=0;
                    while(getline(inHeat, strmesh)){  // 
                        if(strmesh.find("Nodes:")>0  // 
                        && strmesh.find("Nodes:")<strmesh.length()){ // 
                            //yufei: extract the number of heatnodes:
                            int flag = 0;
                            for (int i = 0; i<strmesh.length(); i++){
                                if (strmesh[i]>='0' && strmesh[i]<='9'){
                                    numheatnode = numheatnode * 10 + strmesh[i] - '0';
                                    flag = 1;
                                }
                                if (flag == 1 && (strmesh[i]<'0'||strmesh[i]>'9')){
                                    break;
                                }
                            }
                            break;
                        }
                    }

                    cout<<"Number of heat nodes: "<<numheatnode<<endl;
                    this->Heatnodes.resize(numheatnode); 
                    for (int i = 0; i < numheatnode; i++) {
                        getline(inHeat, strmesh);
                        stringstream sss; 
                        sss.str(strmesh);  
                        string coord_x, coord_y, coord_z, coord_heat;
                        sss >> coord_x >> coord_y >> coord_z >> coord_heat; 
                        Heatnodes[i].x = strtod(coord_x.c_str(), NULL) * L_x; 
                        Heatnodes[i].y = strtod(coord_y.c_str(), NULL) * L_y;  
                        Heatnodes[i].z = strtod(coord_z.c_str(), NULL) * L_z;;
                        Heatnodes[i].heat = strtod(coord_heat.c_str(), NULL);  
                        //cout << "Heatnodes" << Heatnodes[i].x << " " <<  Heatnodes[i].y << " " << Heatnodes[i].z << " " << Heatnodes[i].heat << endl;
                    }

                    //yufei: test whether heatnode is consistent with node of mesh
                    double distances[Heatnodes.size()]={0.0};
                    double selectMindistance[8][2]={0.0}; //yufei:this is to store the min distance and heat, for 3D, consider 8 nearest neibour
                    for (int i=0; i<Nodes.size(); i++){    
                        int k=0;
                        for (int j=0 ;j<Heatnodes.size(); j++){
                            if (Nodes[i].x==Heatnodes[j].x && Nodes[i].y==Heatnodes[j].y && Nodes[i].z==Heatnodes[j].z){
                                Nodes[i].heat=Heatnodes[j].heat;
                                k=1; 
                                break; 
                            }
                        }
                        //yufei: if cannot find the same mesh, calculate the distance
                        if (k==0){
                            for (int j=0 ;j<Heatnodes.size(); j++){ 
                                distances[j] = get_distance(Heatnodes[j],Nodes[i]);
                                //cout << "distance" << distances[j] << endl;
                            } 
                            // yufei: find the nearest 8 heatnode
                            for (int numNearest = 0; numNearest<8; numNearest++){
                                int minPosition = min_element(distances,distances+Heatnodes.size()) - distances;
                                selectMindistance[numNearest][0] = 1/distances[minPosition]; 
                                selectMindistance[numNearest][1] = Heatnodes[minPosition].heat;
                                distances[minPosition]=1000;
                                //cout << "selectMindinstace" << selectMindistance[numNearest][0] << " " << selectMindistance[numNearest][1] << endl;
                            }
                            // calculate heat of the node
                            double totalRevDistance = 0;
                            for (int numNearest = 0; numNearest<8; numNearest++){
                                Nodes[i].heat += selectMindistance[numNearest][0]*selectMindistance[numNearest][1];
                                totalRevDistance += selectMindistance[numNearest][0];
                               // cout<< "TotalRevDistance" << totalRevDistance << endl;
                            }
                            Nodes[i].heat = Nodes[i].heat/totalRevDistance;
                        }
                    }

                //yufei: for 3D mesh, vertex = 4 & 5 & 6 & 8
                    for (int i=0; i<Elements.size(); i++){
                        if (Elements[i].vertexes.size()==4){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat + Nodes[Elements[i].vertexes[3]].heat)/4;
                        }
                        if (Elements[i].vertexes.size()==5){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat + Nodes[Elements[i].vertexes[3]].heat + Nodes[Elements[i].vertexes[4]].heat)/5;
                        }
                        if (Elements[i].vertexes.size()==6){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat + Nodes[Elements[i].vertexes[3]].heat + Nodes[Elements[i].vertexes[4]].heat + Nodes[Elements[i].vertexes[5]].heat)/6;
                        }
                        if (Elements[i].vertexes.size()==8){
                            Elements[i].heat_source=(Nodes[Elements[i].vertexes[0]].heat + Nodes[Elements[i].vertexes[1]].heat + Nodes[Elements[i].vertexes[2]].heat + Nodes[Elements[i].vertexes[3]].heat + Nodes[Elements[i].vertexes[4]].heat + Nodes[Elements[i].vertexes[5]].heat + Nodes[Elements[i].vertexes[6]].heat + Nodes[Elements[i].vertexes[7]].heat)/8;
                        }
                    }
                    cout<<"Successfully read from COORDINATE heat file !"<<endl;
                }
        }
        cout<<"**************************"<<endl;
    }
}
BTEMesh::~BTEMesh()
{
    // cout << "~BTEMesh is activated !!" << endl;
}
