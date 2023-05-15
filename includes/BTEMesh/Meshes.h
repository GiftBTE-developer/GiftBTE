//
// Created by yuehu on 2023/3/6.
//

#ifndef STATICBTESOLVER_MESHES_H
#define STATICBTESOLVER_MESHES_H
#include "utility/utility.h"
#include "BTEMesh.h"
#include "Distribute.h"
class Meshes{
    std::vector<BTEMesh> BTEMeshes;
    Meshes(DistributeMesh & distributeMesh);
    ~Meshes();
};

#endif //STATICBTESOLVER_MESHES_H
