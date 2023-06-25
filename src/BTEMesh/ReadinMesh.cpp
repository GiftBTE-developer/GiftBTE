//
// Created by yuehu on 2023/4/1.
//

#include "BTEMesh/ReadinMesh.h"
using namespace std;

ReadIn::ReadIn(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z, std::string mesh_type) {
    if(mesh_type=="1D")
    {
        OneD(Dimension_Geometry,inFile,L_x,L_y,L_z);
    } else if(mesh_type=="COMSOL")
    {
        COMSOL(Dimension_Geometry,inFile,L_x,L_y,L_z);
    } else if(mesh_type=="MSH")
    {
        MSH(Dimension_Geometry,inFile,L_x,L_y,L_z);
    } else
    {
        cout<<"Error: does not support this mesh type"<<endl;
        exit(0);
    }
}

void ReadIn::OneD(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z) {
    if (Dimension_Geometry==1)
    {
        if (!inFile.is_open()) {
            std::cout << "Error: Mesh file not open" << std::endl;
            exit(1);

        }
        string line;

        getline(inFile,line);
        int num=0;

        int num_geo;
        inFile>>num_geo>>line;
        getline(inFile,line);
        vector<double> distribution(num_geo,0);
        vector<double> meshdistribution(num_geo,0);
        for (int i = 0; i <num_geo ; ++i) {
            getline(inFile,line);
            std::stringstream ss(line);
            int temp;
            ss>>temp>>distribution[i]>>meshdistribution[i];
            num+=meshdistribution[i];
        }
        nodeX.resize(num+1,0);
        nodeY.resize(num+1,0);
        nodeZ.resize(num+1,0);
        volumeElements.resize(num);
        volumeElementIndex.resize(num);

        int index=0;
        for (int i = 0; i <meshdistribution[0]+1 ; ++i) {

            nodeX[index]=(L_x*distribution[0])/meshdistribution[0]*i;
            index++;
        }
        for (int ig = 1; ig <num_geo ; ++ig) {
            for (int i = 1; i <meshdistribution[ig]+1 ; ++i) {

                nodeX[index]=(L_x*(distribution[ig]-distribution[ig-1]))/meshdistribution[ig]*i+L_x*distribution[ig-1];
                index++;
            }
        }

        for (int i = 0; i < num; ++i) {
            volumeElements[i].push_back(i);
            volumeElements[i].push_back(i+1);

        }

        for (int i = 0; i < num ; ++i) {
            if (1.0/2*(nodeX[i]+nodeX[i+1])<distribution[0]*L_x)
            {
                volumeElementIndex[i]=0;
            }
            for (int j = 1; j < num_geo ; ++j) {
                if (1.0/2*(nodeX[i]+nodeX[i+1])<distribution[j]*L_x&&1.0/2*(nodeX[i]+nodeX[i+1])>distribution[j-1]*L_x)
                {
                    volumeElementIndex[i]=j;
                }
            }
        }

        boundaryElements.resize(num_geo+1);
        boundaryElementIndex.resize(num_geo+1);

        for (int i = 0; i < num_geo; ++i) {
            int startnode=0;
            int endnode=distribution[0];
            for (int j = 1; j <=i ; ++j) {
                startnode=startnode+meshdistribution[j-1];
                endnode=endnode+meshdistribution[j];
            }
            boundaryElements[i].push_back(startnode);
        }
        boundaryElements[num_geo].push_back(num);


        for (int i = 0; i < boundaryElementIndex.size() ; ++i) {
            boundaryElementIndex[i]=i+1;//jiaxuan
        }
    }
    else
    {
        cout<<"Error: oneD is only for one dimensional mesh"<<endl;
        exit(0);
    }

}

void ReadIn::COMSOL(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z) {
    if (Dimension_Geometry == 1)//jiaxuan
    {
        if (!inFile.is_open())
        {
            cout << "Error: Mesh file not open" << endl;
            exit(1);
        }
        string line;
        string strmesh;

        getline(inFile,strmesh);
        if (strmesh.find("COMSOL") > 0 && strmesh.find("COMSOL") < strmesh.length()) {
            int dim;
            while (getline(inFile, strmesh)) {
                if (strmesh.find("# sdim") > 0 && strmesh.find("# sdim") < strmesh.length()) {
                    dim = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            if (dim != Dimension_Geometry) {
                cout << "Error: please check the dimension is correct or not" << endl;
                exit(1);
            }
            int numnode;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of mesh") > 0 && strmesh.find("# number of mesh") < strmesh.length())
                {
                    numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            nodeX.resize(numnode);
            nodeY.resize(numnode);
            nodeZ.resize(numnode);
            for (int i = 0; i < 3; i++)
                getline(inFile, strmesh);
            for (int i = 0; i < numnode; i++)
            {
                inFile>>nodeX[i];
                nodeZ[i]=0;
                nodeY[i]=0;
                nodeX[i]=nodeX[i]*L_x;
            }
            int numoftype;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of element types") > 0 && strmesh.find("# number of element types") < strmesh.length())
                {
                    numoftype = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            while (numoftype > 0)
            {
                numoftype--;
                int numvert;
                while (getline(inFile, strmesh))
                {
                    if (strmesh.find("# number of vertices per element") > 0 && strmesh.find("# number of vertices per element") < strmesh.length() || strmesh.find("# number of nodes per element") > 0 && strmesh.find("# number of nodes per element") < strmesh.length())
                    {
                        numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                if (numvert == 1)
                {
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    boundaryElements.resize(numvert);
                    getline(inFile, strmesh);
                    boundaryElementIndex.resize(numvert);
                    for (int i = 0; i < numvert; i++)
                    {
                        boundaryElements[i].resize(1);
                        //inFile >> boundaryElements[i][0] >> boundaryElements[i][1];
                        inFile >> boundaryElements[i][0];
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = 0; i < numvert; i++)
                    {
                        inFile >>  boundaryElementIndex[i];
                        boundaryElementIndex[i] = boundaryElementIndex[i] + 1;//jiaxuan
                    }

                }
                else if (numvert == 2)
                {
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numvert; ++i)
                    {
                        vector<int> c1;
                        c1.resize(2);
                        inFile >> c1[0] >> c1[1];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numvert; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
            }
        }else
        {
            cout<<"Error: mesh file format error"<<endl;
            exit(1);
        }

    }
    else if (Dimension_Geometry == 2)
    {
        // vector<vector<int>> Boundaryindex;
        if (!inFile.is_open())
        {
            cout << "Error: Mesh file not open" << endl;
            exit(1);
        }


        string strmesh;

        getline(inFile,strmesh);
        if (strmesh.find("COMSOL") > 0 && strmesh.find("COMSOL") < strmesh.length())
        {
            int dim;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# sdim") > 0 && strmesh.find("# sdim") < strmesh.length())
                {
                    dim = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            if(dim!=Dimension_Geometry)
            {
                cout<<"Error: please check the dimension is correct or not"<<endl;
                exit(1);
            }
            int numnode;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of mesh") > 0 && strmesh.find("# number of mesh") < strmesh.length())
                {
                    numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            nodeX.resize(numnode);
            nodeY.resize(numnode);
            nodeZ.resize(numnode);
            for (int i = 0; i < 3; i++)
                getline(inFile, strmesh);
            for (int i = 0; i < numnode; i++)
            {
                inFile>>nodeX[i]>>nodeY[i];
                nodeZ[i]=0;
                nodeX[i]=nodeX[i]*L_x;
                nodeY[i]=nodeY[i]*L_y;
            }
            int numoftype;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of element types") > 0 && strmesh.find("# number of element types") < strmesh.length())
                {
                    numoftype = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            while (numoftype > 0)
            {
                numoftype--;
                int numvert;
                while (getline(inFile, strmesh))
                {
                    if (strmesh.find("# number of vertices per element") > 0 && strmesh.find("# number of vertices per element") < strmesh.length() || strmesh.find("# number of nodes per element") > 0 && strmesh.find("# number of nodes per element") < strmesh.length())
                    {
                        numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                        break;
                    }
                }
                if (numvert == 2)
                {
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    boundaryElements.resize(numvert);
                    getline(inFile, strmesh);
                    boundaryElementIndex.resize(numvert);
                    for (int i = 0; i < numvert; i++)
                    {
                        boundaryElements[i].resize(2);
                        inFile >> boundaryElements[i][0] >> boundaryElements[i][1];
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = 0; i < numvert; i++)
                    {
                        inFile >>  boundaryElementIndex[i];
                        boundaryElementIndex[i] = boundaryElementIndex[i] + 1;//jiaxuan
                    }
                }
                else if (numvert == 3)
                {
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numvert; ++i)
                    {
                        vector<int> c1;
                        c1.resize(3);
                        inFile >> c1[0] >> c1[1] >> c1[2];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numvert; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
                else if (numvert == 4)
                {
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numvert = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numvert; ++i)
                    {
                        vector<int> c1;
                        c1.resize(4);
                        inFile >> c1[0] >> c1[1] >> c1[2] >> c1[3];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numvert; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
            }
        } else
        {
            cout<<"Error: mesh file format error"<<endl;
            exit(1);
        }

    }
    else if (Dimension_Geometry == 3)
    {
        if (!inFile.is_open())
        {
            cout << "Error: Mesh file not open" << endl;
            exit(1);
        }
        std::string strmesh;
        getline(inFile,strmesh);
        if (strmesh.find("COMSOL") > 0 && strmesh.find("COMSOL") < strmesh.length())
        {
            int dim;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# sdim") > 0 && strmesh.find("# sdim") < strmesh.length())
                {
                    dim = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            if(dim!=Dimension_Geometry)
            {
                cout<<"Error: please check the dimension is correct or not"<<endl;
                exit(1);
            }
            int numnode;

            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of mesh") > 0 && strmesh.find("# number of mesh") < strmesh.length())
                {
                    numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            nodeX.resize(numnode);
            nodeY.resize(numnode);
            nodeZ.resize(numnode);
            for (int i = 0; i < 3; i++)
                getline(inFile, strmesh);
            for (int i = 0; i < numnode; i++)
            {
                inFile>>nodeX[i]>>nodeY[i]>>nodeZ[i];

                nodeX[i]=nodeX[i]*L_x;
                nodeY[i]=nodeY[i]*L_y;
                nodeZ[i]=nodeZ[i]*L_z;
            }
            int numoftype;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("# number of element types") > 0 && strmesh.find("# number of element types") < strmesh.length())
                {
                    numoftype = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }
            while (numoftype > 0)
            {
                numoftype--;
                string nameoftype;
                while (getline(inFile, strmesh))
                {
                    if (strmesh.find("# type name") > 0 && strmesh.find("# type name") < strmesh.length())
                    {
                        stringstream sss;
                        sss.str(strmesh);
                        int three;
                        sss >> three >> nameoftype;
                        break;
                    }
                }
                if (nameoftype == "quad")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> bc;
                        bc.resize(4);
                        inFile >> bc[0] >> bc[1] >> bc[2] >> bc[3];
                        //   Boundaryindex.resize(numvert);
                        boundaryElements.push_back(bc);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = boundaryElements.size() - numc; i < boundaryElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        boundaryElementIndex.push_back(m+1);//jiaxuan
                    }
                }
                else if (nameoftype == "tri")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> bc;
                        bc.resize(3);
                        inFile >> bc[0] >> bc[1] >> bc[2];
                        boundaryElements.push_back(bc);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = boundaryElements.size() - numc; i < boundaryElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        boundaryElementIndex.push_back(m+1);//jiaxuan
                    }
                }
                else if (nameoftype == "hex")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> c1;
                        c1.resize(8);
                        inFile >> c1[0] >> c1[1] >> c1[2] >> c1[3] >> c1[4] >> c1[5] >> c1[6] >> c1[7];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numc; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
                else if (nameoftype == "tet")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> c1;
                        c1.resize(4);
                        inFile >> c1[0] >> c1[1] >> c1[2] >> c1[3];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numc; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
                else if (nameoftype == "pyr")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> c1;
                        c1.resize(8);
                        inFile >> c1[0] >> c1[1] >> c1[2] >> c1[3] >> c1[4];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numc; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
                else if (nameoftype == "prism")
                {
                    int numc;
                    while (getline(inFile, strmesh))
                    {
                        if (strmesh.find("# number of elements") > 0 && strmesh.find("# number of elements") < strmesh.length())
                        {
                            numc = stoi(strmesh.substr(0, strmesh.find(" ")));
                            break;
                        }
                    }
                    getline(inFile, strmesh);
                    for (int i = 0; i < numc; i++)
                    {
                        vector<int> c1;
                        c1.resize(6);
                        inFile >> c1[0] >> c1[1] >> c1[2] >> c1[3] >> c1[4] >> c1[5];
                        volumeElements.push_back(c1);
                    }
                    for (int i = 0; i < 4; i++)
                        getline(inFile, strmesh);
                    for (int i = volumeElements.size() - numc; i < volumeElements.size(); i++)
                    {
                        int m;
                        inFile>>m;
                        volumeElementIndex.push_back(m);
                    }
                }
            }
        } else
        {
            cout<<"Error: mesh file format error"<<endl;
            exit(1);
        }
        if(volumeElements.size()==0)
        {
            cout<<"Error: no element, please check the dimension is correct or not"<<endl;
            exit(1);
        }
    }
    for (int i = 0; i < volumeElementIndex.size() ; ++i) {
        volumeElementIndex[i]=volumeElementIndex[i]-1;
    }
}


void ReadIn::MSH(int Dimension_Geometry, ifstream &inFile, double L_x, double L_y, double L_z) {
    if (Dimension_Geometry == 2)
    {
        std::string strmesh;
        getline(inFile, strmesh);
        if (strmesh.find("MeshFormat") > 0 && strmesh.find("MeshFormat") < strmesh.length())
        {
            int numnode;
            int numcell;
            int numblock;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("Nodes") > 0 && strmesh.find("Nodes") < strmesh.length())
                {
                    inFile>>numblock>>numnode;
                    getline(inFile, strmesh);


                    nodeX.resize(numnode);
                    nodeY.resize(numnode);
                    nodeZ.resize(numnode);
                    for (int i = 0; i < numblock; ++i) {
                        int dim;
                        int tag_all;
                        int para;
                        int numnodeinblock;

                        inFile>>dim>>tag_all>>para>>numnodeinblock;
                        vector<int> tag;
                        tag.resize(numnodeinblock);
                        for (int j = 0; j < numnodeinblock; ++j) {
                            inFile>>tag[j];
                        }
                        for (int j = 0; j < numnodeinblock; ++j) {
                            Point a;
                            inFile>>a.x>>a.y>>a.z;
                            nodeX[tag[j]-1]=a.x*L_x;
                            nodeY[tag[j]-1]=a.y*L_y;
                            nodeZ[tag[j]-1]=a.z*L_z;
                        }
                    }
                    //cout<<"node ok"<<endl;
                    //numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }

            while (getline(inFile, strmesh))
            {
                if (strmesh.find("Elements") > 0 && strmesh.find("Elements") < strmesh.length())
                {
                    inFile>>numblock>>numcell;
                    getline(inFile, strmesh);
                    for (int i = 0; i < numblock; ++i) {
                        int dim;
                        int tag_all;
                        int type;
                        int numnodeinblock;
                        inFile>>dim>>tag_all>>type>>numnodeinblock;
                        getline(inFile, strmesh);
                        for (int j = 0; j < numnodeinblock; ++j) {
                            if(type==1)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(2);
                                inFile>>tag_local>>a[0]>>a[1];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                boundaryElementIndex.push_back(tag_all);//jiaxuan
                                boundaryElements.push_back(a);
                            } else if (type==2)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(3);
                                inFile>>tag_local>>a[0]>>a[1]>>a[2];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                volumeElements.push_back(a);
                                volumeElementIndex.push_back(tag_all-1);
                            } else if (type==3)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(4);
                                inFile>>tag_local>>a[0]>>a[1]>>a[3]>>a[2];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                a[3]=a[3]-1;
                                volumeElements.push_back(a);
                                volumeElementIndex.push_back(tag_all-1);
                            } else
                            {
                                getline(inFile, strmesh);
                            }
                            if(dim==3)
                            {
                                cout<<"Error: find three dimensional element,please check the dimension is correct or not"<<endl;
                                exit(1);
                            }
                        }
                    }
                    //cout<<"node ok"<<endl;
                    //numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }

            }
        } else
        {
            cout<<"Error: mesh file format error"<<endl;
            exit(1);
        }

    }
    else if(Dimension_Geometry == 3)
    {

        std::string strmesh;
        getline(inFile, strmesh);
        if (strmesh.find("MeshFormat") > 0 && strmesh.find("MeshFormat") < strmesh.length())
        {
            int numnode;
            int numcell;
            int numblock;
            while (getline(inFile, strmesh))
            {
                if (strmesh.find("Nodes") > 0 && strmesh.find("Nodes") < strmesh.length())
                {
                    inFile>>numblock>>numnode;
                    getline(inFile, strmesh);


                    nodeX.resize(numnode);
                    nodeY.resize(numnode);
                    nodeZ.resize(numnode);
                    for (int i = 0; i < numblock; ++i) {
                        int dim;
                        int tag_all;
                        int para;
                        int numnodeinblock;

                        inFile>>dim>>tag_all>>para>>numnodeinblock;
                        vector<int> tag;
                        tag.resize(numnodeinblock);
                        for (int j = 0; j < numnodeinblock; ++j) {
                            inFile>>tag[j];
                        }
                        for (int j = 0; j < numnodeinblock; ++j) {
                            Point a;
                            inFile>>a.x>>a.y>>a.z;
                            nodeX[tag[j]-1]=a.x*L_x;
                            nodeY[tag[j]-1]=a.y*L_y;
                            nodeZ[tag[j]-1]=a.z*L_z;
                        }
                    }
                    //cout<<"node ok"<<endl;
                    //numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }
            }

            while (getline(inFile, strmesh))
            {
                if (strmesh.find("Elements") > 0 && strmesh.find("Elements") < strmesh.length())
                {
                    inFile>>numblock>>numcell;
                    getline(inFile, strmesh);
                    for (int i = 0; i < numblock; ++i) {
                        int dim;
                        int tag_all;
                        int type;
                        int numnodeinblock;
                        inFile>>dim>>tag_all>>type>>numnodeinblock;
                        getline(inFile, strmesh);
                        for (int j = 0; j < numnodeinblock; ++j) {
                            if (type==2)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(3);
                                inFile>>tag_local>>a[0]>>a[1]>>a[2];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                boundaryElements.push_back(a);
                                boundaryElementIndex.push_back(tag_all);//jiaxuan
                            } else if (type==3)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(4);
                                inFile>>tag_local>>a[0]>>a[1]>>a[3]>>a[2];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                a[3]=a[3]-1;
                                boundaryElements.push_back(a);
                                boundaryElementIndex.push_back(tag_all);//jiaxuan
                            } else if(type==4)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(4);
                                inFile>>tag_local>>a[0]>>a[1]>>a[2]>>a[3];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                a[3]=a[3]-1;
                                volumeElements.push_back(a);
                                volumeElementIndex.push_back(tag_all-1);
                            }
                            else if(type==5)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(8);
                                inFile>>tag_local>>a[0]>>a[1]>>a[3]>>a[2]>>a[4]>>a[5]>>a[7]>>a[6];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                a[3]=a[3]-1;
                                a[4]=a[4]-1;
                                a[5]=a[5]-1;
                                a[6]=a[6]-1;
                                a[7]=a[7]-1;
                                volumeElements.push_back(a);
                                volumeElementIndex.push_back(tag_all-1);
                            }
                            else if(type==6)
                            {
                                int tag_local;
                                vector<int> a;
                                a.resize(6);
                                inFile>>tag_local>>a[3]>>a[4]>>a[5]>>a[0]>>a[1]>>a[2];
                                a[0]=a[0]-1;
                                a[1]=a[1]-1;
                                a[2]=a[2]-1;
                                a[3]=a[3]-1;
                                a[4]=a[4]-1;
                                a[5]=a[5]-1;
                                volumeElements.push_back(a);
                                volumeElementIndex.push_back(tag_all-1);
                            }
                            else
                            {
                                getline(inFile, strmesh);
                            }
                        }
                    }
                    //cout<<"node ok"<<endl;
                    //numnode = stoi(strmesh.substr(0, strmesh.find(" ")));
                    break;
                }

            }
        } else
        {
            cout<<"Error: mesh file format error"<<endl;
            exit(1);
        }
        if(volumeElements.size()==0)
        {
            cout<<"Error: no element,please check the dimension is correct or not"<<endl;
            exit(1);
        }

    }
}

ReadIn::~ReadIn() {

}

