//
// Created by yuehu on 2021/9/8.
//

#ifndef STATICBTESOLVER_UTILITY_H
#define STATICBTESOLVER_UTILITY_H

#include "cmath"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include "mpi.h"
#include "sstream"

#ifndef USE_GPU
//#include <petscksp.h>
#endif

#define PI M_PI
#define Eps pow(2, -52)
#define my_min(a, b) (((a) < (b)) ? (a) : (b))



class vec {
public:
    double x = 0;
    double y = 0;
    double z = 0;
    explicit vec(double xx = 0, double yy = 0, double zz = 0) : x{xx}, y{yy}, z{zz} {}
    double operator*(const vec &rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }
    vec operator*(const double rhs) const
    {
        return vec(x * rhs, y * rhs, z * rhs);
    }
    vec operator/(const double rhs) const
    {
        return vec(x / rhs, y / rhs, z / rhs);
    }

    vec operator-(const vec rhs) const
    {
        return vec(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }
};

class Point {
public:
    double x = 0;
    double y = 0;
    double z = 0;
    double heat = 0; //yufei adding
    std::vector<int> cells;
    Point &operator=(const Point &rhs) = default;
    explicit Point(double xx = 0, double yy = 0, double zz = 0) : x{xx}, y{yy}, z{zz} {}
    Point operator+(const Point &rhs) const
    {
        return Point(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z);
    }
    Point operator*(const double rhs) const
    {
        return Point(this->x * rhs, this->y * rhs, this->z * rhs);
    }
    Point operator/(const double rhs) const
    {
        return Point(this->x / rhs, this->y / rhs, this->z / rhs);
    }
    vec operator-(const Point &rhs) const
    {
        return vec(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    double distance(const Point &p) const
    {
        return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z));
    }
};

class Face {
public:
    double area;
    int neighbor;
    int bound;
    std::vector<int> vertexes;
    vec norm;
    vec direction;
    Point center;
    int index;
};

class Cell {
public:
    double volume;
    int index;
    int matter;
    double heat_source;
    std::vector<Face> faces;
    std::vector<int> vertexes;
    Point center;
};

class Boundarytemp {
public:
    int matter;
    double Temp;
};

class faces {
public:
    int Boundindex;
    std::vector<int> cellindex;
    std::vector<int> faceindex;
};

class Boundary {
public:
    int index;
    int connection;
    int type;
    double Temperature;
    std::vector<Boundarytemp> macro_Temp;
    std::vector<double> emit_temp;
    std::vector<int> cellindex;
    std::vector<int> faceindex;
    std::vector<int> neighbors;
    std::vector<std::vector<double>> trans;
};

class ContinuousArray {
public:
    // double *data{};
    // int dim1{};
    // int dim2{};
    double *data;
    int dim1;
    int dim2;

    ContinuousArray() = default;

    // ContinuousArray(int dim1, int dim2)
    // {
    //     init(dim1, dim2);
    // }
    ContinuousArray(int d1, int d2)
    {
        // this->dim1 = d1;
        // this->dim2 = d2;
        dim1 = d1;
        dim2 = d2;
        data = new double[d1 * d2];
    }

    ContinuousArray(const ContinuousArray &rhs)
    {
        // this->dim1 = rhs.dim1;
        // this->dim2 = rhs.dim2;
        // data = new double[this->dim1 * this->dim2];
        // std::copy(rhs.data, rhs.data + rhs.dim1 * rhs.dim2, data);
        dim1 = rhs.dim1;
        dim2 = rhs.dim2;
        data = new double[dim1 * dim2];
        std::copy(rhs.data, rhs.data + rhs.dim1 * rhs.dim2, data);
    }
 
    ContinuousArray &operator=(const ContinuousArray &rhs)
    {
        if (this != &rhs)
        {
            // delete[] this->data;
            dim1 = rhs.dim1;
            dim2 = rhs.dim2;
            data = new double[dim1 * dim2];
            std::copy(rhs.data, rhs.data + rhs.dim1 * rhs.dim2, data);
        }
        return *this;
    }

    void init(int d1, int d2)
    {
        // delete[] data;
        //  data = new double[d1 * d2]();
        //  this->dim1 = d1;
        //  this->dim2 = d2;
        dim1 = d1;
        dim2 = d2;
        data = new double[d1 * d2]();
    }

    double *get_ptr(int i, int j) const
    {
        // return data + (i * this->dim2 + j);
        return data + (i * dim2 + j);
    }
    void set(int i, int j, double m) const
    {
        // std::cout<<dim1<<" "<<dim2<<std::endl;
        // std::fill_n(this->data + (i * this->dim2 + j), 1, m);
        // std::cout<<"i and j result"<<std::endl;
        // std::fill_n(data + (i * dim2 + j), 1, m);
        // std::cout<<data[i * dim2 + j]<<std::endl;
        data[i * dim2 + j] = m;
        // std::cout<<data[i * dim2 + j]<<std::endl;
        // std::cout<<"i and j result"<<std::endl;
    }

    double get(int i, int j) const
    {
      
        // return *get_ptr(i, j);
        return data[i * dim2 + j];
    }

    void clear() const
    {
        // std::fill_n(this->data, this->dim1 * this->dim2, 0.0);
        std::fill_n(data, dim1 * dim2, 0.0);
    }

    void print() const
    {
        for (int i = 0; i < dim1; i++)
        {
            for (int j = 0; j < dim2; j++)
            {
                std::cout << data[i * dim2 + j] << "  ";
            }
            std::cout << std::endl;
        }
    }

    ~ContinuousArray()
    {
        delete[] data;
    }
};

double get_distance(Point &a, Point &b);
double get_length(vec &a);
double dot_product(vec &a, vec &b);
double dot_product(std::vector<double> &a, std::vector<double> &b);
double detJ(double a);
double get_tet_volume(std::vector<Point> &v1);
double get_tri_area(std::vector<Point> &tri);
double sgn(double a);
vec minus(Point &a, Point &b);
vec cross_product(vec &p1, vec &p2);
std::vector<double> GaussIntegrationPoints(double a, double b, int N);
size_t get_host_memory();
vec Get_gradient_2D(std::vector<std::vector<double>> &J, std::vector<double> &d);
vec Get_gradient_3D(std::vector<std::vector<double>> &J, std::vector<double> &d);
void Get_inverse_2D(double **Jfin, std::vector<std::vector<double>> &J);
void Get_inverse_3D(double **Jfin,std::vector<std::vector<double>> &J);
#ifndef USE_GPU
double *_solve_matrix(int *csrRowPtr, int *csrColInd, double *csrVal, std::vector<double> &Re, int N_cell);
double *_solve_matrix_1(int *csrRowPtr, int *csrColInd, double *csrVal, std::vector<double> &Re, int N_cell);
double *_solve_matrix(int *RowPtr, int *ColInd, double *Val, double *Re, int N_cell);
double *_solve_matrix_hypre(int *csrRowPtr, int *csrColInd, double *csrVal, std::vector<double> &Re, int N_cell, int world_rank, int num_proc);
double *_solve_matrix_hypre(int *RowPtr, int *ColInd, double *Val, double *Re, int N_cell, int world_rank, int num_proc);
bool ishave(int a,const std::vector<int>& array);

#endif

#endif 
// STATICBTESOLVER_UTILITY_H

// Iterations = 26
// Final Relative Residual Norm = 6.749481e-07

// internal error0.690644

// Iterations = 26
// Final Relative Residual Norm = 6.749481e-07

// internal error2.96443e-16
// ----------------------------------------------------------------------------------
// Iteration #28	 Margin per band per cell: Temperature:9.82116e-06 Heat Flux:9.59354e-06
//   Time taken by inner loop: 99.661 milliseconds
//   Time taken by gradient 1: 13.673 milliseconds
//   Time taken by BTE solver: 18.13 milliseconds
//   Time taken by Solving matrix: 16.31 milliseconds
//   Time taken by get_temp_flux: 0.08 milliseconds
//   Time taken by Boundary: 0.252 milliseconds 0.968 milliseconds

//   Time taken by non_fourier: 42.291 milliseconds
//   Time taken by macro: 18.009 milliseconds

//   Time taken by Solving matrix: 9.597 milliseconds
//   Time taken by transfer: 6.15 milliseconds 0.024 milliseconds
// ----------------------------------------------------------------------------------

// internal error3.40534e-16
// ----------------------------------------------------------------------------------
// Iteration #28    Margin per band per cell: Temperature:9.82116e-06 Heat Flux:9.5935e-06
//   Time taken by inner loop: 286.174 milliseconds
//   Time taken by gradient 1: 47.62 milliseconds
//   Time taken by BTE solver: 74.998 milliseconds
//   Time taken by Solving matrix: 67.303 milliseconds
//   Time taken by get_temp_flux: 0.391 milliseconds
//   Time taken by Boundary: 0.916 milliseconds 2.607 milliseconds

//   Time taken by non_fourier: 142.59 milliseconds
//   Time taken by macro: 16.731 milliseconds

//   Time taken by Solving matrix: 9.855 milliseconds
//   Time taken by transfer: 0.069 milliseconds 0.003 milliseconds
// ----------------------------------------------------------------------------------
// end_iter
// end_iter1
// Time taken by iteration: 9715.87 milliseconds
// BTEBand析构函数被调用
// ~BTEBoundaryCondition is activated!!
// ~BTEMesh is activated !!
// BTEAngle析构函数被调用
// ~StaticBTESynthetic is activated !!

/*
Iterations = 36
Final Relative Residual Norm = 9.476718e-07

internal error3.40534e-16
----------------------------------------------------------------------------------
Iteration #28    Margin per band per cell: Temperature:9.82116e-06 Heat Flux:9.5935e-06
  Time taken by inner loop: 282.131 milliseconds
  Time taken by gradient 1: 45.95 milliseconds
  Time taken by BTE solver: 71.698 milliseconds
  Time taken by Solving matrix: 65.165 milliseconds
  Time taken by get_temp_flux: 0.35 milliseconds
  Time taken by Boundary: 0.905 milliseconds 2.595 milliseconds

  Time taken by non_fourier: 143.307 milliseconds
  Time taken by macro: 16.955 milliseconds

  Time taken by Solving matrix: 9.951 milliseconds
  Time taken by transfer: 0.096 milliseconds 0.004 milliseconds
----------------------------------------------------------------------------------
end_iter
end_iter1
Time taken by iteration: 10091.8 milliseconds
BTEBand析构函数被调用
~BTEBoundaryCondition is activated!!
~BTEMesh is activated !!
BTEAngle析构函数被调用
~StaticBTESynthetic is activated !!
*/
// Iterations = 36
// Final Relative Residual Norm = 9.476718e-07

// internal error0.690642

// Iterations = 36
// Final Relative Residual Norm = 9.476718e-07

// internal error3.40534e-16
// ----------------------------------------------------------------------------------
// Iteration #28    Margin per band per cell: Temperature:9.82116e-06 Heat Flux:9.5935e-06
//   Time taken by inner loop: 281.056 milliseconds
//   Time taken by gradient 1: 44.194 milliseconds
//   Time taken by BTE solver: 73.682 milliseconds
//   Time taken by Solving matrix: 66.656 milliseconds
//   Time taken by get_temp_flux: 0.39 milliseconds
//   Time taken by Boundary: 1.021 milliseconds 2.626 milliseconds
//   Time taken by non_fourier: 141.883 milliseconds
//   Time taken by macro: 16.88 milliseconds
//   Time taken by Solving matrix: 9.966 milliseconds
//   Time taken by transfer: 0.098 milliseconds 0.005 milliseconds
// ----------------------------------------------------------------------------------
// end_iter
// end_iter1

// siyuanyihao
// 4.77944e-16


// GPU result
// Iteration #42	 Margin per band per cell: Temperature:5.89557e-05
// internal error0.690617
// ----------------------------------------------------------------------------------
// Iteration #43	 Margin per band per cell: Temperature:1.46325e-05
// internal error0.690634
// ----------------------------------------------------------------------------------
// Iteration #44	 Margin per band per cell: Temperature:4.07054e-05
// internal error0.690629
// ----------------------------------------------------------------------------------
// Iteration #45	 Margin per band per cell: Temperature:1.13693e-05
// internal error0.690641
// ----------------------------------------------------------------------------------
// Iteration #46	 Margin per band per cell: Temperature:2.81463e-05
// internal error0.690637
// ----------------------------------------------------------------------------------
// Iteration #47	 Margin per band per cell: Temperature:8.70658e-06
