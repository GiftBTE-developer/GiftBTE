//
// Created by yuehu on 2023/3/6.
//

#include "BTEMesh/Distribute.h"

#include "utility/utility.h"
#include <fstream>
DistributeMesh::DistributeMesh(int Dimension_Geometry, std::ifstream &inFile, double L_x, double L_y, double L_z,
                               std::string mesh_type,BTEBand * bands, BTEBoundaryCondition *bcs,
                               std::ifstream &inHeat,double Uniform_Heat, std::string heat_type,std::ifstream &initialTemp) {  //yufei jiaxaun

    ReadIn readIn(Dimension_Geometry, inFile, L_x, L_y, L_z, mesh_type);
    nodeXFourier=readIn.nodeX;
    nodeYFourier=readIn.nodeY;
    nodeZFourier=readIn.nodeZ;
    volumeElementsFourier=readIn.volumeElements;
    volumeElementIndexFourier=readIn.volumeElementIndex;
    boundaryElementsFourier=readIn.boundaryElements;
    boundaryElementIndexFourier=readIn.boundaryElementIndex;

    BTEMesh Fouriermesh1(Dimension_Geometry,L_x,L_y,L_z,nodeXFourier,nodeYFourier,nodeZFourier,volumeElementsFourier
                         ,volumeElementIndexFourier,boundaryElementsFourier,boundaryElementIndexFourier); 
    Fouriermesh1.setMeshParams(bands);
    Fouriermesh1.setMeshParams1(bcs);
    Fouriermesh1.BTEMesh_heatin(inHeat, Uniform_Heat, heat_type); //yufei
    Fouriermesh1.BTEMesh_initialTemp(initialTemp);//jiaxuan
    FourierMeshes=Fouriermesh1;


}
void DistributeMesh::_build_BTEMesh(int Dimension_Geometry,double L_x, double L_y, double L_z,BTEBand * bands, BTEBoundaryCondition *bcs, std::ifstream &inHeat,double Uniform_Heat,std::string mesh_solver_file)
{
    int num_geo=0;
    vector<int> Allregion;

    for (int i = 0; i < volumeElementsFourier.size(); ++i) {
        if(!ishave(volumeElementIndexFourier[i],Allregion))
        {
            Allregion.push_back(volumeElementIndexFourier[i]);
        }
    }

    ifstream inmesh_solver(mesh_solver_file);

    num_geo=Allregion.size();
    vector<string> Mesh_solver(num_geo);
    for (int i = 0; i < num_geo; ++i) {
        Mesh_solver[i]="BTE";
    }
    if(inmesh_solver.is_open())
    {
        for (int i = 0; i < num_geo; ++i) {
            int m;
            inmesh_solver>>m>>Mesh_solver[m];
        }
    }
    vector<vector<int>> nodeElementIndex(nodeXFourier.size());
    vector<vector<int>> nodeGeoIndex(nodeXFourier.size());
    vector<vector<int>> boundaryVolumeIndex(boundaryElementsFourier.size());
    vector<vector<int>> boundaryGeoIndex(boundaryElementsFourier.size());
    vector<vector<pair<int,int>>> boundaryEachIndex(boundaryElementsFourier.size());



    for (int i = 0; i <volumeElementsFourier.size(); ++i) {
        for (int j = 0; j < volumeElementsFourier[i].size(); ++j) {
            nodeElementIndex[volumeElementsFourier[i][j]].push_back(i);
            if(!ishave(volumeElementIndexFourier[i], nodeGeoIndex[volumeElementsFourier[i][j]]))
            {
                nodeGeoIndex[volumeElementsFourier[i][j]].push_back(volumeElementIndexFourier[i]);
            }
        }
    }

    for (int i = 0; i < boundaryElementsFourier.size(); ++i) {
        vector<int> result=nodeElementIndex[boundaryElementsFourier[i][0]];
        for (int j = 1; j < boundaryElementsFourier[i].size(); ++j) {
            vector<int> result1=result;
            result.resize(0);
            for (int k = 0; k < nodeElementIndex[boundaryElementsFourier[i][j]].size(); ++k) {
                if(ishave(nodeElementIndex[boundaryElementsFourier[i][j]][k], result1))
                {
                    result.push_back(nodeElementIndex[boundaryElementsFourier[i][j]][k]);
                }
            }
        }
        boundaryVolumeIndex[i]=result;
        for (int j = 0; j < result.size(); ++j) {
            boundaryGeoIndex[i].push_back(volumeElementIndexFourier[result[j]]);
        }
    }

    vector<int> projection(nodeGeoIndex.size());

    for (int j = 0; j < nodeGeoIndex.size(); ++j) {
        int is=0;
        for (int k = 0; k < nodeGeoIndex[j].size(); ++k) {
            if(Mesh_solver[nodeGeoIndex[j][k]]=="BTE")
            {
                is=1;
            }
        }
        projection[j]=-1;
        if(is==1)
        {
            nodeXBTE.push_back(nodeXFourier[j]);
            nodeYBTE.push_back(nodeYFourier[j]);
            nodeZBTE.push_back(nodeZFourier[j]);
            projection[j]=nodeXBTE.size()-1;
            nodeProjection.push_back(j);
        }
    }
    nodeProjection=projection;


    for (int j = 0; j < volumeElementIndexFourier.size(); ++j) {
        vector<int> volumeElement(volumeElementsFourier[j].size());
        for (int k = 0; k < volumeElementsFourier[j].size(); ++k) {
            volumeElement[k]=projection[volumeElementsFourier[j][k]];
        }
        if(Mesh_solver[volumeElementIndexFourier[j]]=="BTE")
        {
            volumeElementsBTE.push_back(volumeElement);
            volumeElementIndexBTE.push_back(volumeElementIndexFourier[j]);
            elementProjection.push_back(j);
        }
    }

    for (int j = 0; j < boundaryElementsFourier.size(); ++j) {
        vector<int> boundaryElement(boundaryElementsFourier[j].size());
        for (int k = 0; k < boundaryElementsFourier[j].size(); ++k) {
            boundaryElement[k]=projection[boundaryElementsFourier[j][k]];
        }
        int is=0;
        for (int i = 0; i < boundaryGeoIndex[j].size(); ++i) {
            if(Mesh_solver[boundaryGeoIndex[j][i]]=="BTE")
            {
                is=1;
            }
        }
        if(is==1)
        {
            boundaryElementsBTE.push_back(boundaryElement);
            boundaryElementIndexBTE.push_back(boundaryElementIndexFourier[j]);
            boundaryProjection.push_back(j);
        }

    }

    BTEMesh BTEmesh1(Dimension_Geometry,L_x,L_y,L_z,nodeXBTE,nodeYBTE,nodeZBTE,volumeElementsBTE
                         ,volumeElementIndexBTE,boundaryElementsBTE,boundaryElementIndexBTE);
    BTEmesh1.setMeshParams(bands);
    vector<int> isset;
    isset.resize(BTEmesh1.Boundaries.size());
    vector<Boundary> Boundaries1 = BTEmesh1.Boundaries;
    vector<Boundary> Boundaries;
    for (int i = 0; i < Boundaries1.size(); ++i)
    {
        for (int j = 0; j < bcs->boundaryConditions.size(); ++j)
        {
            if (Boundaries1[i].index == bcs->boundaryConditions[j].index)
            {
               isset[i]=1;
            }
        }
    }
    for (int i = 0; i < Boundaries1.size(); ++i) {
        if(Boundaries1[i].cellindex.size()==2)
        {
            isset[i]=1;
        }
    }
    vector<int> nodefine;
    for (int i = 0; i < Boundaries1.size(); ++i) {
        if(isset[i]==0)
        {
            if(!ishave(Boundaries1[i].index,nodefine))
            {
                nodefine.push_back(Boundaries1[i].index);
            }

        }
    }

    BTEBoundaryCondition bcs1;
    for (int i = 0; i < bcs->boundaryConditions.size(); ++i) {
        bcs1.boundaryConditions.push_back(bcs->boundaryConditions[i]);
    }
    for (int i = 0; i < nodefine.size(); ++i) {
        BoundaryCondition bc;
        bc.index=nodefine[i];
        bc.type=1;
        bc.temperature=0;
        bcs1.boundaryConditions.push_back(bc);
    }
    BTEmesh1.setMeshParams1(&bcs1);
    for (int i = 0; i < BTEmesh1.Elements.size(); ++i) {
        BTEmesh1.Elements[i].heat_source=FourierMeshes.Elements[elementProjection[i]].heat_source;
        BTEmesh1.Elements[i].initial_temperature=FourierMeshes.Elements[elementProjection[i]].initial_temperature;//jiaxuan
    }
    //BTEmesh1.BTEMesh_heatin(inHeat, Uniform_Heat);
    BTEMeshes=BTEmesh1;
    //cout<<endl;
}

DistributeMesh::~DistributeMesh() {

}
