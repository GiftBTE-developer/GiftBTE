#include "utility/utility.h"

__global__ void copy_data_cuda(int numcell, double *Temp_old, double *Temp, double *Re_macro, double *Temp_sythetic, double *Heat_Flux_x, double *Heat_Flux_y, double *Heat_Flux_z)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Temp_old[i] = Temp[i];
        Re_macro[i] = 0;
        // Temp[i]=0;
        Temp_sythetic[i] = 0;
        Heat_Flux_x[i] = 0;
        Heat_Flux_y[i] = 0;
        Heat_Flux_z[i] = 0;
    }
}

__global__ void get_gradientn(int numcell, int iband, double *gradientn, const int *face_num, const double *face_norm, const double *face_center, const int *face_neighbor, const double *Cell_center,const double *ee_n, const double *ee_bound, const int *cell_face_bound)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (face_neighbor[j + i * 6] >= 0)
            {
                double ax, ay, az;
                ax = Cell_center[face_neighbor[j + i * 6] * 3 + 0] - Cell_center[i * 3 + 0];
                ay = Cell_center[face_neighbor[j + i * 6] * 3 + 1] - Cell_center[i * 3 + 1];
                az = Cell_center[face_neighbor[j + i * 6] * 3 + 2] - Cell_center[i * 3 + 2];
                // ax = Cell_centerX[face_neighbor[j + i * 6]] - Cell_centerX[i];
                // ay = Cell_centerY[face_neighbor[j + i * 6]] - Cell_centerY[i];
                // az = Cell_centerZ[face_neighbor[j + i * 6]] - Cell_centerZ[i];
                gradientn[i * 6 + j] = (ee_n[face_neighbor[j + i * 6]] - ee_n[i]) *
                                       (face_norm[0 + j * 3 + i * 3 * 6] * ax + face_norm[1 + j * 3 + i * 3 * 6] * ay + face_norm[2 + j * 3 + i * 3 * 6] * az) / (ax * ax + ay * ay + az * az);
            }
            else
            {
                double ax, ay, az;
                ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 0];
                ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 1];
                az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 2];
                gradientn[i * 6 + j] = (ee_bound[cell_face_bound[i * 6 + j]] - ee_n[i]) *
                                       (face_norm[0 + j * 3 + i * 3 * 6] * ax + face_norm[1 + j * 3 + i * 3 * 6] * ay + face_norm[2 + j * 3 + i * 3 * 6] * az) / (ax * ax + ay * ay + az * az);
            }
        }
    }
}
__global__ void get_gradientn_split(int numcell, int iband, double *gradientn, const int *face_num, const double *face_normX,const double *face_normY,const double *face_normZ, const double *face_centerX,const double *face_centerY,const double *face_centerZ, const int *face_neighbor, const double *Cell_centerX,const double *Cell_centerY,const double *Cell_centerZ,const double *ee_n, const double *ee_bound, const int *cell_face_bound)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (face_neighbor[j + i * 6] >= 0)
            {
                double ax, ay, az;
                // ax = Cell_center[face_neighbor[j + i * 6] * 3 + 0] - Cell_center[i * 3 + 0];
                // ay = Cell_center[face_neighbor[j + i * 6] * 3 + 1] - Cell_center[i * 3 + 1];
                // az = Cell_center[face_neighbor[j + i * 6] * 3 + 2] - Cell_center[i * 3 + 2];
                ax = Cell_centerX[face_neighbor[j + i * 6]] - Cell_centerX[i];
                ay = Cell_centerY[face_neighbor[j + i * 6]] - Cell_centerY[i];
                az = Cell_centerZ[face_neighbor[j + i * 6]] - Cell_centerZ[i];
                // gradientn[i * 6 + j] = (ee_n[face_neighbor[j + i * 6]] - ee_n[i]) *
                //                        (face_norm[0 + j * 3 + i * 3 * 6] * ax + face_norm[1 + j * 3 + i * 3 * 6] * ay + face_norm[2 + j * 3 + i * 3 * 6] * az) / (ax * ax + ay * ay + az * az);
                gradientn[i * 6 + j] = (ee_n[face_neighbor[j + i * 6]] - ee_n[i]) *
                                       (face_normX[j + i * 6] * ax + face_normY[j + i * 6] * ay + face_normZ[j + i * 6] * az) / (ax * ax + ay * ay + az * az);
            }
            else
            {
                double ax, ay, az;
                // ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 0];
                // ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 1];
                // az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 2];
                ax = face_centerX[j  + i * 6] - Cell_centerX[i];
                ay = face_centerY[j  + i * 6] - Cell_centerY[i];
                az = face_centerZ[j  + i * 6] - Cell_centerZ[i];
                // gradientn[i * 6 + j] = (ee_bound[cell_face_bound[i * 6 + j]] - ee_n[i]) *
                //                        (face_norm[0 + j * 3 + i * 3 * 6] * ax + face_norm[1 + j * 3 + i * 3 * 6] * ay + face_norm[2 + j * 3 + i * 3 * 6] * az) / (ax * ax + ay * ay + az * az);
                gradientn[i * 6 + j] = (ee_bound[cell_face_bound[i * 6 + j]] - ee_n[i]) *
                                       (face_normX[j + i * 6] * ax + face_normY[ j  + i * 6] * ay + face_normZ[j + i * 6] * az) / (ax * ax + ay * ay + az * az);
            }
        }
    }
}

__global__ void get_macro_R_bound(int numcell, int inf, int iband, double kappabulk, double *local_R, const int *face_num, const double *face_norm, const double *face_center,
                                  const int *face_neighbor, const double *ss1, const int *boundtype, const int *face_vertex_size,
                                  const int *face_vertexes, const double *nodes, const double *gradientn,
                                  const double *face_area, const double *weight, const double *tau_n, const double *vg_n, const double *C_n, const double *R_n,
                                  const double *ee_n, const double *Cell_center, double *J, double *JT, double *Jmul, double *JT1, double *d, double *gradientx1, double *gradienty1, double *gradientz1)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (boundtype[-face_neighbor[j + i * 6] - 1] != 3)
            {
                // int index1 = i * 6 * 3;
                // int index2 = i * 6;
                double gradientx, gradienty, gradientz;
                gradientx = gradientx1[i];
                gradienty = gradienty1[i];
                gradientz = gradientz1[i];
                local_R[i] -=
                    face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] *
                    vg_n[iband] * (ss1[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1[inf * 3 + 0] * gradientx + ss1[inf * 3 + 1] * gradienty + ss1[inf * 3 + 2] * gradientz);
                local_R[i] +=
                    face_area[j + i * 6] * weight[inf] * R_n[iband] / C_n[iband] *
                    (gradientx * face_norm[0 + j * 3 + i * 3 * 6] + gradienty * face_norm[1 + j * 3 + i * 3 * 6] + gradienty * face_norm[2 + j * 3 + i * 3 * 6]) * kappabulk;
            }
        }
    }
}
__global__ void get_macro_R_bound_split(int numcell, int inf, int iband, double kappabulk, double *local_R, const int *face_num, const double *face_norm, const double *face_center,
                                  const int *face_neighbor, const double *ss1X,const double *ss1Y,const double *ss1Z, const int *boundtype, const int *face_vertex_size,
                                  const int *face_vertexes, const double *nodes, const double *gradientn,
                                  const double *face_area, const double *weight, const double *tau_n, const double *vg_n, const double *C_n, const double *R_n,
                                  const double *ee_n, const double *Cell_center, double *J, double *JT, double *Jmul, double *JT1, double *d, double *gradientx1, double *gradienty1, double *gradientz1)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (boundtype[-face_neighbor[j + i * 6] - 1] != 3)
            {
                // int index1 = i * 6 * 3;
                // int index2 = i * 6;
                double gradientx, gradienty, gradientz;
                gradientx = gradientx1[i];
                gradienty = gradienty1[i];
                gradientz = gradientz1[i];
                local_R[i] -=
                    face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] *
                    // vg_n[iband] * (ss1[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1[inf * 3 + 0] * gradientx + ss1[inf * 3 + 1] * gradienty + ss1[inf * 3 + 2] * gradientz);
                    vg_n[iband] * (ss1X[inf] * face_norm[0 + j * 3 + i * 3 * 6] + ss1Y[inf] * face_norm[1 + j * 3 + i * 3 * 6] + ss1Z[inf] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1X[inf] * gradientx + ss1Y[inf] * gradienty + ss1Z[inf] * gradientz);
                local_R[i] +=
                    face_area[j + i * 6] * weight[inf] * R_n[iband] / C_n[iband] *
                    (gradientx * face_norm[0 + j * 3 + i * 3 * 6] + gradienty * face_norm[1 + j * 3 + i * 3 * 6] + gradienty * face_norm[2 + j * 3 + i * 3 * 6]) * kappabulk;
            }
        }
    }
}
__global__ void get_local_R_1(int numcell, int inf, int iband, double kappabulk, double *local_R, const double *face_norm,
                              const int *face_neighbor, const double *ss1, const int *face_vertex_size,
                              const int *face_vertexes, const double *nodes, const double *ee_vertex, const double *gradientn,
                              const double *face_area, const double *weight, const double *tau_n, const double *vg_n, const double *C_n, const double *R_n, double *J, double *JT, double *Jmul, double *JT1, double *d, const int *boundtype, const int *cell_face_bound)
{
    // int itotal = threadIdx.x+blockIdx.x*blockDim.x;

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    // for (; itotal < numcell*6 ; itotal+=blockDim.x*gridDim.x)  {
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        // int j=itotal%6;
        // int i=(itotal-j)/6;
        for (int j = 0; j < 6; ++j)
        {
            if (boundtype[cell_face_bound[j + i * 6]] != 3)
            {
                double gradientx = 0, gradienty = 0, gradientz = 0;
                int index1 = (i)*4 * 3;
                int index2 = (i)*4 * 3;
                int index3 = (i)*4;
                J[index1 + 0 * 3 + 0] = face_norm[0 + j * 3 + i * 3 * 6];
                J[index1 + 0 * 3 + 1] = face_norm[1 + j * 3 + i * 3 * 6];
                J[index1 + 0 * 3 + 2] = face_norm[2 + j * 3 + i * 3 * 6];

                d[index3 + 0] = 0;
                for (int m = 0; m < face_vertex_size[j + i * 6] - 1; m++)
                {
                    J[index1 + (m + 1) * 3 + 0] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 0] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 0];
                    J[index1 + (m + 1) * 3 + 1] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 1] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 1];
                    J[index1 + (m + 1) * 3 + 2] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 2] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 2];
                    d[index3 + m + 1] = ee_vertex[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - ee_vertex[face_vertexes[m + j * 4 + i * 6 * 4]];
                }
                for (int m = 0; m < face_vertex_size[j + i * 6]; ++m)
                {
                    for (int n = 0; n < 3; ++n)
                    {
                        JT[index1 + n * face_vertex_size[j + i * 6] + m] = J[index1 + m * 3 + n];
                    }
                }
                for (int m = 0; m < 3; ++m)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        Jmul[index2 + m * 3 + k] = 0;
                        for (int n = 0; n < face_vertex_size[j + i * 6]; ++n)
                        {
                            Jmul[index2 + m * 3 + k] += JT[index1 + m * face_vertex_size[j + i * 6] + n] * J[index1 + n * 3 + k];
                        }
                        // cout<<Jmul[m*3+k]<<" ";
                    }
                    // cout<<endl;
                }
                // cout<<endl;

                double detJ = Jmul[index2 + 0 * 3 + 0] * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]) + Jmul[index2 + 0 * 3 + 1] * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]) + Jmul[index2 + 0 * 3 + 2] * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);

                J[index1 + 0 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]);
                J[index1 + 1 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 2]);
                J[index1 + 2 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 1 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 1 * 3 + 1]);
                J[index1 + 0 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]);
                J[index1 + 1 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 0]);
                J[index1 + 2 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 2] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 2]);
                J[index1 + 0 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);
                J[index1 + 1 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 1]);
                J[index1 + 2 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 1] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 1]);

                for (int m = 0; m < 3; ++m)
                {
                    for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
                    {
                        JT1[index1 + m * face_vertex_size[j + i * 6] + k] = 0;
                        for (int n = 0; n < 3; ++n)
                        {
                            JT1[index1 + m * face_vertex_size[j + i * 6] + k] += J[index1 + m * 3 + n] * JT[index1 + n * face_vertex_size[j + i * 6] + k];
                        }
                        // cout<<JT1[m*face_vertex_size[j+i*6]+k]<<" ";
                    }
                    // cout<<endl;
                }
                for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
                {
                    gradientx += JT1[index1 + 0 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                    gradienty += JT1[index1 + 1 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                    gradientz += JT1[index1 + 2 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                }

                // VectorXd gradient1(3);
                gradientx = gradientn[i * 6 + j] * face_norm[0 + j * 3 + i * 3 * 6] + gradientx;
                gradienty = gradientn[i * 6 + j] * face_norm[1 + j * 3 + i * 3 * 6] + gradienty;
                gradientz = gradientn[i * 6 + j] * face_norm[2 + j * 3 + i * 3 * 6] + gradientz;
                // gradientx=gradientn[i*6+j]*face_norm[0+j*3+i*3*6];
                // gradienty=gradientn[i*6+j]*face_norm[1+j*3+i*3*6];
                // gradientz=gradientn[i*6+j]*face_norm[2+j*3+i*3*6];
                local_R[i] -= face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] * vg_n[iband] *
                              (ss1[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1[inf * 3 + 0] * gradientx + ss1[inf * 3 + 1] * gradienty + ss1[inf * 3 + 2] * gradientz);
                local_R[i] += face_area[j + i * 6] * weight[inf] * R_n[iband] / C_n[iband] * gradientn[i * 6 + j] * kappabulk;
            }
        }
    }
}
__global__ void get_local_R_1_split(int numcell, int inf, int iband, double kappabulk, double *local_R, const double *face_normX, const double *face_normY, const double *face_normZ,
                              const int *face_neighbor, const double *ss1X,const double *ss1Y,const double *ss1Z, const int *face_vertex_size,
                              const int *face_vertexes, const double *nodesX,const double *nodesY,const double *nodesZ, const double *ee_vertex, const double *gradientn,
                              const double *face_area, const double *weight, const double *tau_n, const double *vg_n, const double *C_n, const double *R_n, double *J, double *JT, double *Jmul, double *JT1, double *d, const int *boundtype, const int *cell_face_bound)
{
    // int itotal = threadIdx.x+blockIdx.x*blockDim.x;

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    // for (; itotal < numcell*6 ; itotal+=blockDim.x*gridDim.x)  {
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        // int j=itotal%6;
        // int i=(itotal-j)/6;
        for (int j = 0; j < 6; ++j)
        {
            if (boundtype[cell_face_bound[j + i * 6]] != 3)
            {
                double gradientx = 0, gradienty = 0, gradientz = 0;
                int index1 = (i)*4 * 3;
                int index2 = (i)*4 * 3;
                int index3 = (i)*4;
                J[index1 + 0 * 3 + 0] = face_normX[j  + i * 6];
                J[index1 + 0 * 3 + 1] = face_normY[j  + i * 6];
                J[index1 + 0 * 3 + 2] = face_normZ[j  + i * 6];

                d[index3 + 0] = 0;
                for (int m = 0; m < face_vertex_size[j + i * 6] - 1; m++)
                {
                    // J[index1 + (m + 1) * 3 + 0] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 0] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 0];
                    // J[index1 + (m + 1) * 3 + 1] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 1] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 1];
                    // J[index1 + (m + 1) * 3 + 2] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 2] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 2];
                    J[index1 + (m + 1) * 3 + 0] = nodesX[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - nodesX[face_vertexes[m + j * 4 + i * 6 * 4]];
                    J[index1 + (m + 1) * 3 + 1] = nodesY[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - nodesY[face_vertexes[m + j * 4 + i * 6 * 4]];
                    J[index1 + (m + 1) * 3 + 2] = nodesZ[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - nodesZ[face_vertexes[m + j * 4 + i * 6 * 4]];
                    d[index3 + m + 1] = ee_vertex[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - ee_vertex[face_vertexes[m + j * 4 + i * 6 * 4]];
                }
                for (int m = 0; m < face_vertex_size[j + i * 6]; ++m)
                {
                    for (int n = 0; n < 3; ++n)
                    {
                        JT[index1 + n * face_vertex_size[j + i * 6] + m] = J[index1 + m * 3 + n];
                    }
                }
                for (int m = 0; m < 3; ++m)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        Jmul[index2 + m * 3 + k] = 0;
                        for (int n = 0; n < face_vertex_size[j + i * 6]; ++n)
                        {
                            Jmul[index2 + m * 3 + k] += JT[index1 + m * face_vertex_size[j + i * 6] + n] * J[index1 + n * 3 + k];
                        }
                        // cout<<Jmul[m*3+k]<<" ";
                    }
                    // cout<<endl;
                }
                // cout<<endl;

                double detJ = Jmul[index2 + 0 * 3 + 0] * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]) + Jmul[index2 + 0 * 3 + 1] * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]) + Jmul[index2 + 0 * 3 + 2] * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);

                J[index1 + 0 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]);
                J[index1 + 1 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 2]);
                J[index1 + 2 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 1 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 1 * 3 + 1]);
                J[index1 + 0 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]);
                J[index1 + 1 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 0]);
                J[index1 + 2 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 2] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 2]);
                J[index1 + 0 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);
                J[index1 + 1 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 1]);
                J[index1 + 2 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 1] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 1]);

                for (int m = 0; m < 3; ++m)
                {
                    for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
                    {
                        JT1[index1 + m * face_vertex_size[j + i * 6] + k] = 0;
                        for (int n = 0; n < 3; ++n)
                        {
                            JT1[index1 + m * face_vertex_size[j + i * 6] + k] += J[index1 + m * 3 + n] * JT[index1 + n * face_vertex_size[j + i * 6] + k];
                        }
                        // cout<<JT1[m*face_vertex_size[j+i*6]+k]<<" ";
                    }
                    // cout<<endl;
                }
                for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
                {
                    gradientx += JT1[index1 + 0 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                    gradienty += JT1[index1 + 1 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                    gradientz += JT1[index1 + 2 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
                }

                // VectorXd gradient1(3);
                gradientx = gradientn[i * 6 + j] * face_normX[j  + i * 6] + gradientx;
                gradienty = gradientn[i * 6 + j] * face_normY[j  + i * 6] + gradienty;
                gradientz = gradientn[i * 6 + j] * face_normZ[j  + i * 6] + gradientz;
                // gradientx=gradientn[i*6+j]*face_norm[0+j*3+i*3*6];
                // gradienty=gradientn[i*6+j]*face_norm[1+j*3+i*3*6];
                // gradientz=gradientn[i*6+j]*face_norm[2+j*3+i*3*6];
                // local_R[i] -= face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] * vg_n[iband] *
                //               (ss1[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1[inf * 3 + 0] * gradientx + ss1[inf * 3 + 1] * gradienty + ss1[inf * 3 + 2] * gradientz);
                local_R[i] -= face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] * vg_n[iband] *
                              (ss1X[inf] * face_normX[j  + i * 6] + ss1Y[inf] * face_normY[j  + i * 6] + ss1Z[inf] * face_normZ[j  + i * 6]) * (ss1X[inf] * gradientx + ss1Y[inf] * gradienty + ss1Z[inf] * gradientz);
                local_R[i] += face_area[j + i * 6] * weight[inf] * R_n[iband] / C_n[iband] * gradientn[i * 6 + j] * kappabulk;
            }
        }
    }
}
__global__ void ini_macro_temp(int Boundnum, double *macro_temp)
{
    int ib = threadIdx.x + blockIdx.x * blockDim.x;
    for (; ib < Boundnum; ib += blockDim.x * gridDim.x)
    {
        macro_temp[ib] = 0;
    }
}
__global__ void get_macro_temp(int Boundnum, double *macro_temp, int iband, int inf, const double *ee_bound, const double *lattice_ratio, const double *heat_capacity, const double *weight)
{
    int ib = threadIdx.x + blockIdx.x * blockDim.x;
    for (; ib < Boundnum; ib += blockDim.x * gridDim.x)
    {
        macro_temp[ib] += lattice_ratio[iband] * ee_bound[ib] / heat_capacity[iband] * weight[inf];
    }
}
__global__ void _get_gradient(int numcell, int nftot, int nband, int iband, int inf, double *gradientx,
                              double *gradienty, double *gradientz, const int *face_num, const int *face_neighbor, const double *Cell_center, const double *ee_n, double *J, double *JT, double *Jmul, double *JT1, double *d)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        int index1 = i * 6 * 3;
        int index2 = i * 6;
        gradientx[i] = 0;
        gradienty[i] = 0;
        gradientz[i] = 0;
        int nump = 0;
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (face_neighbor[j + i * 6] >= 0)
            {
                nump++;
            }
        }
        nump = 0;
        for (int m = 0; m < face_num[i]; ++m)
        {

            if (face_neighbor[m + i * 6] >= 0)
            {
                J[index1 + nump * 3 + 0] = Cell_center[face_neighbor[m + i * 6] * 3 + 0] -
                                           Cell_center[i * 3 + 0];
                J[index1 + nump * 3 + 1] = Cell_center[face_neighbor[m + i * 6] * 3 + 1] -
                                           Cell_center[i * 3 + 1];
                J[index1 + nump * 3 + 2] = Cell_center[face_neighbor[m + i * 6] * 3 + 2] -
                                           Cell_center[i * 3 + 2];
                d[index2 + nump] = ee_n[face_neighbor[m + i * 6]] -
                                   ee_n[i];
                // cout<<J(nump,0)<<" "<<J(nump,1)<<" "<<J(nump,2)<<" "<<d(nump)<<endl;

                nump++;
            }
        }
        for (int m = 0; m < nump; ++m)
        {
            for (int n = 0; n < 3; ++n)
            {
                JT[index1 + n * nump + m] = J[index1 + m * 3 + n];
            }
        }
        for (int m = 0; m < 3; ++m)
        {
            for (int k = 0; k < 3; ++k)
            {
                Jmul[index1 + m * 3 + k] = 0;
                for (int n = 0; n < nump; ++n)
                {
                    Jmul[index1 + m * 3 + k] += JT[index1 + m * nump + n] * J[index1 + n * 3 + k];
                }
                // cout<<Jmul[m*3+k]<<" ";
            }
            // cout<<endl;
        }
        // cout<<endl;

        double detJ = Jmul[index1 + 0 * 3 + 0] * (Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 1]) + Jmul[index1 + 0 * 3 + 1] * (Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 2]) + Jmul[index1 + 0 * 3 + 2] * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 0]);
        if (detJ != 0)
        {
            J[index1 + 0 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 1]);
            J[index1 + 1 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 2 * 3 + 2]);
            J[index1 + 2 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 1 * 3 + 2] - Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 1 * 3 + 1]);
            J[index1 + 0 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 2]);
            J[index1 + 1 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 2 * 3 + 0]);
            J[index1 + 2 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 0 * 3 + 2] - Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 1 * 3 + 2]);
            J[index1 + 0 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 0]);
            J[index1 + 1 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 2 * 3 + 1]);
            J[index1 + 2 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 1 * 3 + 1] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 0 * 3 + 1]);

            for (int m = 0; m < 3; ++m)
            {
                for (int k = 0; k < nump; ++k)
                {
                    JT1[index1 + m * nump + k] = 0;
                    for (int n = 0; n < 3; ++n)
                    {
                        JT1[index1 + m * nump + k] += J[index1 + m * 3 + n] * JT[index1 + n * nump + k];
                    }
                    // cout<<JT1[m*face_vertex_size[j+i*6]+k]<<" ";
                }
                // cout<<endl;
            }
            for (int k = 0; k < nump; ++k)
            {
                gradientx[i] += JT1[index1 + 0 * nump + k] * d[index2 + k];
                gradienty[i] += JT1[index1 + 1 * nump + k] * d[index2 + k];
                gradientz[i] += JT1[index1 + 2 * nump + k] * d[index2 + k];
            }
            // VectorXd gradient1(3);
        }
    }
}
__global__ void _get_gradient_split(int numcell, int nftot, int nband, int iband, int inf, double *gradientx,
                              double *gradienty, double *gradientz, const int *face_num, const int *face_neighbor, const double *Cell_centerX,const double *Cell_centerY,const double *Cell_centerZ, const double *ee_n, double *J, double *JT, double *Jmul, double *JT1, double *d)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        int index1 = i * 6 * 3;
        int index2 = i * 6;
        gradientx[i] = 0;
        gradienty[i] = 0;
        gradientz[i] = 0;
        int nump = 0;
        for (int j = 0; j < face_num[i]; ++j)
        {
            if (face_neighbor[j + i * 6] >= 0)
            {
                nump++;
            }
        }
        nump = 0;
        for (int m = 0; m < face_num[i]; ++m)
        {

            if (face_neighbor[m + i * 6] >= 0)
            {
                // J[index1 + nump * 3 + 0] = Cell_center[face_neighbor[m + i * 6] * 3 + 0] -
                //                            Cell_center[i * 3 + 0];
                // J[index1 + nump * 3 + 1] = Cell_center[face_neighbor[m + i * 6] * 3 + 1] -
                //                            Cell_center[i * 3 + 1];
                // J[index1 + nump * 3 + 2] = Cell_center[face_neighbor[m + i * 6] * 3 + 2] -
                //                            Cell_center[i * 3 + 2];
                J[index1 + nump * 3 + 0] = Cell_centerX[face_neighbor[m + i * 6]] -
                                           Cell_centerX[i];
                J[index1 + nump * 3 + 1] = Cell_centerY[face_neighbor[m + i * 6]] -
                                           Cell_centerY[i];
                J[index1 + nump * 3 + 2] = Cell_centerZ[face_neighbor[m + i * 6]] -
                                           Cell_centerZ[i];
                d[index2 + nump] = ee_n[face_neighbor[m + i * 6]] -
                                   ee_n[i];
                // cout<<J(nump,0)<<" "<<J(nump,1)<<" "<<J(nump,2)<<" "<<d(nump)<<endl;

                nump++;
            }
        }
        for (int m = 0; m < nump; ++m)
        {
            for (int n = 0; n < 3; ++n)
            {
                JT[index1 + n * nump + m] = J[index1 + m * 3 + n];
            }
        }
        for (int m = 0; m < 3; ++m)
        {
            for (int k = 0; k < 3; ++k)
            {
                Jmul[index1 + m * 3 + k] = 0;
                for (int n = 0; n < nump; ++n)
                {
                    Jmul[index1 + m * 3 + k] += JT[index1 + m * nump + n] * J[index1 + n * 3 + k];
                }
                // cout<<Jmul[m*3+k]<<" ";
            }
            // cout<<endl;
        }
        // cout<<endl;

        double detJ = Jmul[index1 + 0 * 3 + 0] * (Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 1]) + Jmul[index1 + 0 * 3 + 1] * (Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 2]) + Jmul[index1 + 0 * 3 + 2] * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 0]);
        if (detJ != 0)
        {
            J[index1 + 0 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 1]);
            J[index1 + 1 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 2 * 3 + 2]);
            J[index1 + 2 * 3 + 0] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 1 * 3 + 2] - Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 1 * 3 + 1]);
            J[index1 + 0 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 2] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 2]);
            J[index1 + 1 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 2 * 3 + 2] - Jmul[index1 + 0 * 3 + 2] * Jmul[index1 + 2 * 3 + 0]);
            J[index1 + 2 * 3 + 1] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 0 * 3 + 2] - Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 1 * 3 + 2]);
            J[index1 + 0 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 2 * 3 + 1] - Jmul[index1 + 1 * 3 + 1] * Jmul[index1 + 2 * 3 + 0]);
            J[index1 + 1 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 1] * Jmul[index1 + 2 * 3 + 0] - Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 2 * 3 + 1]);
            J[index1 + 2 * 3 + 2] = 1.0 / detJ * (Jmul[index1 + 0 * 3 + 0] * Jmul[index1 + 1 * 3 + 1] - Jmul[index1 + 1 * 3 + 0] * Jmul[index1 + 0 * 3 + 1]);

            for (int m = 0; m < 3; ++m)
            {
                for (int k = 0; k < nump; ++k)
                {
                    JT1[index1 + m * nump + k] = 0;
                    for (int n = 0; n < 3; ++n)
                    {
                        JT1[index1 + m * nump + k] += J[index1 + m * 3 + n] * JT[index1 + n * nump + k];
                    }
                    // cout<<JT1[m*face_vertex_size[j+i*6]+k]<<" ";
                }
                // cout<<endl;
            }
            for (int k = 0; k < nump; ++k)
            {
                gradientx[i] += JT1[index1 + 0 * nump + k] * d[index2 + k];
                gradienty[i] += JT1[index1 + 1 * nump + k] * d[index2 + k];
                gradientz[i] += JT1[index1 + 2 * nump + k] * d[index2 + k];
            }
            // VectorXd gradient1(3);
        }
    }
}

__global__ void _BTE_iteration_cell(int numcell, int nftot, int nband, int iband, int inf, const double *group_velocity, double *relaxation_time, const double *heat_capacity, const double *heat_ratio, const int *face_num, const double *face_area, const double *face_norm, const double *face_center, const double *Cell_center, const double *Cell_volume, const double *ss, const double *weight, const double *ee_n, const double *gradientx, const double *gradienty, const double *gradientz,
                                    const int *face_neighbor, const double *Temp, const double *heat_source, int Boundnum, int *Bound_cell, int *Bound_face, int *cell_face_bound, double *Re)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Re[i] = 0;
        for (int j = 0; j < face_num[i]; ++j)
        {
            double temp = weight[inf] * group_velocity[iband] * relaxation_time[iband] * face_area[j + i * 6] / Cell_volume[i] *
                          (ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]);
            double e = 0;
            if ((ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) >= 0)
            {
                double ax, ay, az;
                ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 0];
                ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 1];
                az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 2];

                e = ee_n[i] + (ax * gradientx[i] + ay * gradienty[i] + az * gradientz[i]);
            }
            else
            {
                if (face_neighbor[j + i * 6] >= 0)
                {
                    int neiindex = face_neighbor[j + i * 6];
                    double ax, ay, az;
                    ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 0];
                    ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 1];
                    az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 2];
                    e = ee_n[face_neighbor[j + i * 6]] + (ax * gradientx[neiindex] + ay * gradienty[neiindex] + az * gradientz[neiindex]);
                }
            }
            Re[i] -= temp * e;
        }
        Re[i] += Temp[i] * heat_capacity[iband] / 4 / PI * weight[inf];
        Re[i] += heat_source[i] * heat_ratio[iband] / 4 / PI * weight[inf] * relaxation_time[iband];
        Re[i] -= ee_n[i] * weight[inf];
    }
}
__global__ void _BTE_iteration_cell_split(int numcell, int nftot, int nband, int iband, int inf, const double *group_velocity, double *relaxation_time, const double *heat_capacity, const double *heat_ratio, const int *face_num, const double *face_area, const double *face_normX, const double *face_normY,const double *face_normZ,const double *face_centerX,const double *face_centerY,const double *face_centerZ, const double *Cell_centerX,const double *Cell_centerY,const double *Cell_centerZ, const double *Cell_volume, const double *ssX,const double *ssY,const double *ssZ, const double *weight, const double *ee_n, const double *gradientx, const double *gradienty, const double *gradientz,
                                    const int *face_neighbor, const double *Temp, const double *heat_source, int Boundnum, int *Bound_cell, int *Bound_face, int *cell_face_bound, double *Re)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Re[i] = 0;
        for (int j = 0; j < face_num[i]; ++j)
        {
            // double temp = weight[inf] * group_velocity[iband] * relaxation_time[iband] * face_area[j + i * 6] / Cell_volume[i] *
            //               (ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]);
            double temp = weight[inf] * group_velocity[iband] * relaxation_time[iband] * face_area[j + i * 6] / Cell_volume[i] *
                          (ssX[inf] * face_normX[j  + i  * 6] + ssY[inf] * face_normY[ j + i * 6] + ssZ[inf] * face_normZ[j  + i * 6]);
            double e = 0;
            // if ((ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) >= 0)
            if ((ssX[inf] * face_normX[j  + i  * 6] + ssY[inf] * face_normY[j  + i  * 6] + ssZ[inf] * face_normZ[j  + i  * 6]) >= 0)
            {
                double ax, ay, az;
                // ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 0];
                // ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 1];
                // az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[i * 3 + 2];
                ax = face_centerX[j + i * 6] - Cell_centerX[i];
                ay = face_centerY[j + i * 6] - Cell_centerY[i];
                az = face_centerZ[j + i * 6] - Cell_centerZ[i];
                e = ee_n[i] + (ax * gradientx[i] + ay * gradienty[i] + az * gradientz[i]);
            }
            else
            {
                if (face_neighbor[j + i * 6] >= 0)
                {
                    int neiindex = face_neighbor[j + i * 6];
                    double ax, ay, az;
                    // ax = face_center[0 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 0];
                    // ay = face_center[1 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 1];
                    // az = face_center[2 + j * 3 + i * 6 * 3] - Cell_center[neiindex * 3 + 2];
                    ax = face_centerX[j + i * 6] - Cell_centerX[neiindex];
                    ay = face_centerY[j + i * 6] - Cell_centerY[neiindex];
                    az = face_centerZ[j + i * 6] - Cell_centerZ[neiindex];
                    e = ee_n[face_neighbor[j + i * 6]] + (ax * gradientx[neiindex] + ay * gradienty[neiindex] + az * gradientz[neiindex]);
                }
            }
            Re[i] -= temp * e;
        }
        Re[i] += Temp[i] * heat_capacity[iband] / 4 / PI * weight[inf];
        Re[i] += heat_source[i] * heat_ratio[iband] / 4 / PI * weight[inf] * relaxation_time[iband];
        Re[i] -= ee_n[i] * weight[inf];
    }
}

__global__ void _get_bound_ee_cuda(int numcell, int nftot, int nband, int iband, int inf, int Boundnum, const int *Bound_cell, int *Bound_face, double *ee_bound, double *face_center, double *Cell_center, const double *face_norm, const double *ss, const double *ee_n, double *gradientx, double *gradienty, double *gradientz)
{
    int ib = threadIdx.x + blockIdx.x * blockDim.x;
    for (; ib < Boundnum; ++ib)
    {
        int i = Bound_cell[ib];
        int j = Bound_face[ib];
        if ((ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) >= 0)
        {

            ee_bound[ib] = ee_n[i];
        }
    }
}
__global__ void _get_bound_ee_cuda_split(int numcell, int nftot, int nband, int iband, int inf, int Boundnum, const int *Bound_cell, int *Bound_face, double *ee_bound, double *face_center, double *Cell_center, const double *face_normX,const double *face_normY,const double *face_normZ, const double *ssX, const double *ssY, const double *ssZ, const double *ee_n, double *gradientx, double *gradienty, double *gradientz)
{
    int ib = threadIdx.x + blockIdx.x * blockDim.x;
    for (; ib < Boundnum; ++ib)
    {
        int i = Bound_cell[ib];
        int j = Bound_face[ib];
        // if ((ss[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) >= 0)
        if ((ssX[inf] * face_normX[j  + i  * 6] + ssY[inf] * face_normY[ j  + i  * 6] + ssZ[inf] * face_normZ[ j  + i  * 6]) >= 0)
        {

            ee_bound[ib] = ee_n[i];
        }
    }
}
__global__ void _BTE_result(int numcell, int nftot, int nband, int iband, int inf, double *ee_n, const double *sol)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        // cout<<sol1[i]<<endl;
        ee_n[i] = sol[i] + ee_n[i];
        // cout<<sol[i]<<endl;
    }
}

__global__ void _get_heat_flux_cuda(int numcell, int nftot, int nband, int iband, int inf, double *Heat_Flux_x, double *Heat_Flux_y,
                                    double *Heat_Flux_z, const double *ee_n, const double *weight, const double *group_velocity, const double *ss)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Heat_Flux_x[i] += group_velocity[iband] * ee_n[i] * ss[inf * 3 + 0] * weight[inf];
        Heat_Flux_y[i] += group_velocity[iband] * ee_n[i] * ss[inf * 3 + 1] * weight[inf];
        Heat_Flux_z[i] += group_velocity[iband] * ee_n[i] * ss[inf * 3 + 2] * weight[inf];
    }
}
__global__ void _get_heat_flux_cuda_split(int numcell, int nftot, int nband, int iband, int inf, double *Heat_Flux_x, double *Heat_Flux_y,
                                    double *Heat_Flux_z, const double *ee_n, const double *weight, const double *group_velocity, const double *ssX,const double *ssY,const double *ssZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Heat_Flux_x[i] += group_velocity[iband] * ee_n[i] * ssX[inf] * weight[inf];
        Heat_Flux_y[i] += group_velocity[iband] * ee_n[i] * ssY[inf] * weight[inf];
        Heat_Flux_z[i] += group_velocity[iband] * ee_n[i] * ssZ[inf] * weight[inf];
    }
}
__global__ void _recover_temperature_cuda(int numcell, int nftot, int nband, int iband, int inf, double *Temp,
                                          const double *lattice_ratio, const double *ee_n, const double *weight, const double *heat_capacity)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    for (; i < numcell; i += blockDim.x * gridDim.x)
    {
        Temp[i] += lattice_ratio[iband] * ee_n[i] * weight[inf] / heat_capacity[iband];
    }
}
// void get_local_R(int numcell, int inf, int iband, double kappabulk, double *local_R, const double *face_norm,
//                  const int *face_neighbor, const double *ss1, const int *face_vertex_size,
//                  const int *face_vertexes, const double *nodes, const double *ee_vertex, const double *gradientn,
//                  const double *face_area, const double *weight, const double *tau_n, const double *vg_n, const double *C_n, const double *R_n, double *J, double *JT, double *Jmul, double *JT1, double *d)
// {
//     // int itotal = threadIdx.x+blockIdx.x*blockDim.x;

//     // for (; itotal < numcell*6 ; itotal+=blockDim.x*gridDim.x)  {
//     for (int i = 0; i < numcell; i++)
//     {
//         // int j=itotal%6;
//         // int i=(itotal-j)/6;
//         for (int j = 0; j < 6; ++j)
//         {
//             if (face_neighbor[j + i * 6] >= 0)
//             {
//                 double gradientx = 0, gradienty = 0, gradientz = 0;
//                 /*int index1=(i)*4*3;
//                 int index2=(i)*4*3;
//                 int index3=(i)*4;*/
//                 int index1 = 0;
//                 int index2 = 0;
//                 int index3 = 0;
//                 // index1=0;
//                 // index2=0;
//                 // index3=0;
//                 J[index1 + 0 * 3 + 0] = face_norm[0 + j * 3 + i * 3 * 6];
//                 J[index1 + 0 * 3 + 1] = face_norm[1 + j * 3 + i * 3 * 6];
//                 J[index1 + 0 * 3 + 2] = face_norm[2 + j * 3 + i * 3 * 6];

//                 d[index3 + 0] = 0;
//                 for (int m = 0; m < face_vertex_size[j + i * 6] - 1; m++)
//                 {
//                     J[index1 + (m + 1) * 3 + 0] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 0] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 0];
//                     J[index1 + (m + 1) * 3 + 1] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 1] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 1];
//                     J[index1 + (m + 1) * 3 + 2] = nodes[face_vertexes[m + 1 + j * 4 + i * 6 * 4] * 3 + 2] - nodes[face_vertexes[m + j * 4 + i * 6 * 4] * 3 + 2];
//                     d[index3 + m + 1] = ee_vertex[face_vertexes[m + 1 + j * 4 + i * 6 * 4]] - ee_vertex[face_vertexes[m + j * 4 + i * 6 * 4]];
//                 }

//                 for (int m = 0; m < face_vertex_size[j + i * 6]; ++m)
//                 {
//                     for (int n = 0; n < 3; ++n)
//                     {
//                         JT[index1 + n * face_vertex_size[j + i * 6] + m] = J[index1 + m * 3 + n];
//                     }
//                 }
//                 for (int m = 0; m < 3; ++m)
//                 {
//                     for (int k = 0; k < 3; ++k)
//                     {
//                         Jmul[index2 + m * 3 + k] = 0;
//                         for (int n = 0; n < face_vertex_size[j + i * 6]; ++n)
//                         {
//                             Jmul[index2 + m * 3 + k] += JT[index1 + m * face_vertex_size[j + i * 6] + n] * J[index1 + n * 3 + k];
//                         }
//                         // cout<<Jmul[m*3+k]<<" ";
//                     }
//                     // cout<<endl;
//                 }
//                 // cout<<endl;

//                 double detJ = Jmul[index2 + 0 * 3 + 0] * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]) + Jmul[index2 + 0 * 3 + 1] * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]) + Jmul[index2 + 0 * 3 + 2] * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);

//                 J[index1 + 0 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 1]);
//                 J[index1 + 1 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 2]);
//                 J[index1 + 2 * 3 + 0] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 1 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 1 * 3 + 1]);
//                 J[index1 + 0 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 2] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 2]);
//                 J[index1 + 1 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 2] - Jmul[index2 + 0 * 3 + 2] * Jmul[index2 + 2 * 3 + 0]);
//                 J[index1 + 2 * 3 + 1] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 2] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 2]);
//                 J[index1 + 0 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 2 * 3 + 1] - Jmul[index2 + 1 * 3 + 1] * Jmul[index2 + 2 * 3 + 0]);
//                 J[index1 + 1 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 1] * Jmul[index2 + 2 * 3 + 0] - Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 2 * 3 + 1]);
//                 J[index1 + 2 * 3 + 2] = 1.0 / detJ * (Jmul[index2 + 0 * 3 + 0] * Jmul[index2 + 1 * 3 + 1] - Jmul[index2 + 1 * 3 + 0] * Jmul[index2 + 0 * 3 + 1]);

//                 for (int m = 0; m < 3; ++m)
//                 {
//                     for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
//                     {
//                         JT1[index1 + m * face_vertex_size[j + i * 6] + k] = 0;
//                         for (int n = 0; n < 3; ++n)
//                         {
//                             JT1[index1 + m * face_vertex_size[j + i * 6] + k] += J[index1 + m * 3 + n] * JT[index1 + n * face_vertex_size[j + i * 6] + k];
//                         }
//                         // cout<<JT1[m*face_vertex_size[j+i*6]+k]<<" ";
//                     }
//                     // cout<<endl;
//                 }
//                 for (int k = 0; k < face_vertex_size[j + i * 6]; ++k)
//                 {
//                     gradientx += JT1[index1 + 0 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
//                     gradienty += JT1[index1 + 1 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
//                     gradientz += JT1[index1 + 2 * face_vertex_size[j + i * 6] + k] * d[index3 + k];
//                 }
//                 // VectorXd gradient1(3);
//                 gradientx = gradientn[i * 6 + j] * face_norm[0 + j * 3 + i * 3 * 6] + gradientx;
//                 gradienty = gradientn[i * 6 + j] * face_norm[1 + j * 3 + i * 3 * 6] + gradienty;
//                 gradientz = gradientn[i * 6 + j] * face_norm[2 + j * 3 + i * 3 * 6] + gradientz;
//                 // gradientx=gradientn[i*6+j]*face_norm[0+j*3+i*3*6];
//                 // gradienty=gradientn[i*6+j]*face_norm[1+j*3+i*3*6];
//                 // gradientz=gradientn[i*6+j]*face_norm[2+j*3+i*3*6];
//                 local_R[i] -= face_area[j + i * 6] * weight[inf] * tau_n[iband] * vg_n[iband] * vg_n[iband] *
//                               (ss1[inf * 3 + 0] * face_norm[0 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 1] * face_norm[1 + j * 3 + i * 3 * 6] + ss1[inf * 3 + 2] * face_norm[2 + j * 3 + i * 3 * 6]) * (ss1[inf * 3 + 0] * gradientx + ss1[inf * 3 + 1] * gradienty + ss1[inf * 3 + 2] * gradientz);
//                 local_R[i] += face_area[j + i * 6] * weight[inf] * R_n[iband] / C_n[iband] * gradientn[i * 6 + j] * kappabulk;
//             }
//         }
//     }
// }

