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
    else if (Dimension_Geometry == 3){
        //delete_Jia
    }
    for (int i = 0; i < volumeElementIndex.size() ; ++i) {
            volumeElementIndex[i]=volumeElementIndex[i]-1;
    }
}

COMSOL::~COMSOL() {

}
