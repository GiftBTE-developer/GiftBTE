//
// Created by yuehu on 2023/3/6.
//

#include "BTEMesh/COMSOL.h"
using namespace std;
COMSOL::COMSOL(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z) {
    if (Dimension_Geometry == 1)
    {
        if (!inFile.is_open())
        {
            cout << "DEBUG: Mesh file not open" << endl;
            exit(1);
        }
        string line;

        getline(inFile, line);
        int num;
        inFile >> num >> line;
        volumeElements.resize(num);
        nodeX.resize(num+1);
        nodeY.resize(num+1);
        nodeZ.resize(num+1);
        for (int i = 0; i < num + 1; ++i)
        {

            nodeX[i] = L_x / num * i;
            nodeY[i]=0;
            nodeZ[i]=0;
        }

        for (int i = 0; i < num; ++i)
        {
            volumeElements[i].push_back(i);
            volumeElements[i].push_back(i+1);
        }
    }
    else if (Dimension_Geometry == 2)
    {
        // vector<vector<int>> Boundaryindex;
        if (!inFile.is_open())
        {
            cout << "DEBUG: Mesh file not open" << endl;
            exit(1);
        }
        string line;
        string strmesh;
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
            getline(inFile, strmesh);
            stringstream sss;
            sss.str(strmesh);
            string coord_x, coord_y;
            sss >> coord_x >> coord_y;
            nodeX[i] = strtod(coord_x.c_str(), NULL) * L_x;
            nodeY[i] = strtod(coord_y.c_str(), NULL) * L_y;
            nodeZ[i] = 0;
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
    }
    else if (Dimension_Geometry == 3)
    {

        std::string strmesh;
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
            getline(inFile, strmesh);
            stringstream sss;
            sss.str(strmesh);
            string coord_x, coord_y, coord_z;
            sss >> coord_x >> coord_y >> coord_z;
            nodeX[i] = strtod(coord_x.c_str(), NULL) * L_x;
            nodeY[i] = strtod(coord_y.c_str(), NULL) * L_y;
            nodeZ[i] = strtod(coord_z.c_str(), NULL) * L_z;
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
                    // nameoftype = strmesh.substr(2, strmesh.find(" "));
                    // cout << nameoftype << endl;
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
                    boundaryElementIndex.push_back(m);
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
                    boundaryElementIndex.push_back(m);
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


    }
    for (int i = 0; i < volumeElementIndex.size() ; ++i) {
            volumeElementIndex[i]=volumeElementIndex[i]-1;
    }
}

COMSOL::~COMSOL() {

}
