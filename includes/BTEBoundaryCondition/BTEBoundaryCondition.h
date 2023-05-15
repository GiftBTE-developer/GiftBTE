//
// Created by yuehu on 2021/9/8.
//

#ifndef STATICBTESOLVER_BTEBOUNDARYCONDITION_H
#define STATICBTESOLVER_BTEBOUNDARYCONDITION_H
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "BTEBand/BTEBand.h"
#include "utility/utility.h"

class transcoeff {
public:
    int matter;
    std::vector<double> trans;
};

class BoundaryCondition {
public:
    int index;
    int type;
    double temperature;
    std::string Name_Interface_file;
    std::vector<transcoeff> matter_trans;

    explicit BoundaryCondition(int index = 0, int type = 0, int temp = 0, std::string name = "") : index(index),
                                                                                                   type(type),
                                                                                                   temperature(temp),
                                                                                                   Name_Interface_file(name) {}

    friend std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc)
    {
        os << "Index: " << bc.index
           << ", Condition: " << bc.type
           << ", Temperature: " << bc.temperature;
        return os;
    }
};

class BTEBoundaryCondition {    
public:
    std::vector<BoundaryCondition> boundaryConditions;

    BTEBoundaryCondition() = default;
    explicit BTEBoundaryCondition(std::ifstream &inFile, BTEBand *bands);
    BTEBoundaryCondition(const BTEBoundaryCondition &other) = default;
    BTEBoundaryCondition &operator=(const BTEBoundaryCondition &other) = default;
    // void copyVector();
    ~BTEBoundaryCondition();
};
#endif // STATICBTESOLVER_BTEBOUNDARYCONDITION_H
