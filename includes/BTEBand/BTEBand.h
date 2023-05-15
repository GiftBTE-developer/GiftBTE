//
// Created by yuehu on 2021/9/8.
//

#ifndef STATICBTESOLVER_BTEBAND_H
#define STATICBTESOLVER_BTEBAND_H
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "utility/utility.h"

class Band {
public:
    std::vector<double> group_velocity;
    std::vector<double> relaxation_time;
    std::vector<double> heat_capacity;
    std::vector<double> lattice_ratio;
    std::vector<double> heat_ratio;

    explicit Band(double vg = 0, double rt = 0, double ct = 0, double Lr = 0, double Hr = 0) : group_velocity(vg),
                                                                                               relaxation_time(rt),
                                                                                               heat_capacity(ct),
                                                                                               lattice_ratio(Lr),
                                                                                               heat_ratio(Hr) {}

    friend std::ostream &operator<<(std::ostream &os, const Band &band)
    {
        for (int i = 0; i < band.group_velocity.size(); ++i)
        {
            os << i << ": "
               << "Group velocity: " << band.group_velocity[i]
               << ", relaxation time: " << band.relaxation_time[i]
               << ", heat_capacity: " << band.heat_capacity[i]
               << ", lattice_ratio: " << band.lattice_ratio[i] 
               << "heat_ratio, : " << band.heat_ratio[i]
               << std::endl;
        }

        return os;
    }
};

class BTEBand {   
public:
    int nband;
    std::vector<Band> bands;
    std::vector<int> geo_matter_index; // in each geometry, which matter is in this part,should be correpsonding to the meshes
    std::vector<double> kappabulk;
    std::vector<double> capacitybulk;

    explicit BTEBand(std::ifstream &inFile, int Dimension_Material);
    BTEBand(const BTEBand &band) = delete;
    // void copyVector();
    ~BTEBand();
};
#endif // STATICBTESOLVER_BTEBAND_H
