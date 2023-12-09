//
// Created by yuehu on 2023/3/27.
//

#include "BTEMesh/oneD.h"
using namespace std;
OneD::OneD(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z) {
    if (Dimension_Geometry==1)
    {
        if (!inFile.is_open()) {
            std::cout << "DEBUG: Mesh file not open" << std::endl;
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
            boundaryElementIndex[i]=i;
        }
    }
    else
    {
        cout<<"oneD is only for one dimensional mesh"<<endl;
        exit(0);
    }

}

OneD::~OneD() {

}
