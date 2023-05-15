//
// Created by yuehu on 2021/9/10.
//
#include "BTEAngle/BTEAngle.h"
#include <iostream>
//     double SolidAngle;
//     vector<vec> direction;
//     vector<double> weight;
//     vector<vec> control_direction;
using namespace std;
BTEAngle::BTEAngle(int ntheta, int nphi, int Dimension_material, int Angle_method)
{
    if (Angle_method == 1) //CADOM
    {
        if (Dimension_material == 2)
        {
            SolidAngle = 2 * PI;
            weight.resize(4 * nphi);
            control_direction.resize(4 * nphi);
            direction.resize(4 * nphi);
            double delta_phi = 0.5 * PI / nphi;
            vector<double> phi(4 * nphi, 0);
            phi[0] = 0.5 * delta_phi;
            for (int np = 1; np < phi.size(); np++)
            {
                phi[np] = phi[np - 1] + delta_phi;
            }

            for (int np = 0; np < phi.size(); np++)
            {
                weight[np] = SolidAngle / phi.size();
                control_direction[np].x = weight[np] * cos(phi[np]);
                control_direction[np].y = weight[np] * sin(phi[np]);
                control_direction[np].z = 0;
                direction[np].x = cos(phi[np]);
                direction[np].y = sin(phi[np]);
                direction[np].z = 0;
            }
        }
        else if (Dimension_material == 3)
        {
            double sum = 0;
            SolidAngle = 4 * PI;
            double dtheta = 0.5 * PI / ntheta;
            double dphi = 0.5 * PI / nphi;
            weight.resize(8 * ntheta * nphi);
            control_direction.resize(8 * ntheta * nphi);
            direction.resize(8 * ntheta * nphi);
            vector<double> theta(2 * ntheta);
            //cout<<"2"<<endl;
            vector<double> phi(4 * nphi);
            theta[0] = 0.5 * dtheta;
            phi[0] = 0.5 * dphi;
            for (int nt = 1; nt < 2 * ntheta; nt++)
                theta[nt] = theta[nt - 1] + dtheta;
            for (int np = 1; np < 4 * nphi; np++)
                phi[np] = phi[np - 1] + dphi;
            for (int nt = 0; nt < 2 * ntheta; nt++)
            {
                for (int np = 0; np < 4 * nphi; np++)
                {
                    int nf = np + nt * 4 * nphi;
                    weight[nf] = 2. * abs(sin(theta[nt])) * sin(0.5 * dtheta) * dphi;
                    control_direction[nf].x = cos(phi[np]) * sin(0.5 * dphi) * (dtheta - cos(2 * theta[nt]) * sin(dtheta));
                    control_direction[nf].y = sin(phi[np]) * sin(0.5 * dphi) * (dtheta - cos(2 * theta[nt]) * sin(dtheta));
                    control_direction[nf].z = 0.5 * dphi * sin(2 * theta[nt]) * sin(dtheta);
                    direction[nf].x = sin(theta[nt]) * cos(phi[np]);
                    direction[nf].y = sin(theta[nt]) * sin(phi[np]);
                    direction[nf].z = cos(theta[nt]);
                    sum += weight[nf];
                }
            }
            //cout<<sum<<endl;
        }
    }
    else //GaussDOM
    {
        if (Dimension_material == 2)
        {
            SolidAngle = 2 * PI;
            weight.resize(4 * nphi);
            control_direction.resize(4 * nphi);
            direction.resize(4 * nphi);
            vector<double> cost, sint;
            vector<double> W;
            vector<double> gauss = GaussIntegrationPoints(0, PI / 2, nphi);
            ;
            for (int dir_index = 0; dir_index < nphi; dir_index++)
            {
                cost.push_back(cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = nphi - 1; dir_index >= 0; dir_index--)
            {
                cost.push_back(-cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = 0; dir_index < nphi; dir_index++)
            {
                cost.push_back(-cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(-pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = nphi - 1; dir_index >= 0; dir_index--)
            {
                cost.push_back(cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(-pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = 0; dir_index < 4 * nphi; dir_index++)
            {
                weight[dir_index] = W[dir_index];
                direction[dir_index].x = cost[dir_index];
                direction[dir_index].y = sint[dir_index];
                direction[dir_index].z = 0;
                control_direction[dir_index].x = cost[dir_index] * weight[dir_index];
                control_direction[dir_index].y = sint[dir_index] * weight[dir_index];
                control_direction[dir_index].z = 0;
            }
        }
        else if (Dimension_material == 3)
        {
            double sum = 0;
            SolidAngle = 4 * PI;

            weight.resize(8 * ntheta * nphi);
            control_direction.resize(8 * ntheta * nphi);
            direction.resize(8 * ntheta * nphi);
            vector<double> cost, cosp, sint, sinp;
            vector<double> W, Wphi;
            vector<double> gauss = GaussIntegrationPoints(-1, 0, ntheta);
            vector<double> gaussp = GaussIntegrationPoints(0, PI / 2, nphi);
            for (int dir_index = 0; dir_index < ntheta; dir_index++)
            {
                cost.push_back(-gauss[2 * dir_index]);
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = ntheta - 1; dir_index >= 0; dir_index--)
            {
                cost.push_back(gauss[2 * dir_index]);
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }

            for (int i = 0; i < nphi; i++)
            {
                cosp.push_back(cos(gaussp[2 * i]));
                Wphi.push_back(gaussp[2 * i + 1]);
                sinp.push_back(pow(1 - cosp[i] * cosp[i], 0.5));
            }
            for (int i = nphi - 1; i >= 0; i--)
            {
                cosp.push_back(-cos(gaussp[2 * i]));
                Wphi.push_back(gaussp[2 * i + 1]);
                sinp.push_back(pow(1 - cosp[i] * cosp[i], 0.5));
            }
            for (int i = 0; i < nphi; i++)
            {
                cosp.push_back(-cos(gaussp[2 * i]));
                Wphi.push_back(gaussp[2 * i + 1]);
                sinp.push_back(-pow(1 - cosp[i] * cosp[i], 0.5));
            }
            for (int i = nphi - 1; i >= 0; i--)
            {
                cosp.push_back(cos(gaussp[2 * i]));
                Wphi.push_back(gaussp[2 * i + 1]);
                sinp.push_back(-pow(1 - cosp[i] * cosp[i], 0.5));
            }
            for (int nt = 0; nt < 2 * ntheta; nt++)
            {
                for (int np = 0; np < 4 * nphi; np++)
                {
                    int nf = np + nt * 4 * nphi;
                    weight[nf] = W[nt] * Wphi[np];
                    control_direction[nf].x = sint[nt] * cosp[np] * weight[nf];
                    control_direction[nf].y = sint[nt] * sinp[np] * weight[nf];
                    control_direction[nf].z = cost[nt] * weight[nf];
                    direction[nf].x = sint[nt] * cosp[np];
                    direction[nf].y = sint[nt] * sinp[np];
                    direction[nf].z = cost[nt];
                    sum += weight[nf];
                }
            }
            //cout<<sum<<endl;
        }
    }
}
BTEAngle::~BTEAngle()
{
    cout << "BTEAngle析构函数被调用" << endl;
}
