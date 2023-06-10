//
// Created by yuehu on 2021/9/8.
//
#include "utility/utility.h"

using namespace std;
bool ishave(int a, const vector<int> &array)
{
    for (int i : array)
        if (a == i)
        {
            return true;
        }
    return false;
}

//获取两点之间的距离
double get_distance(Point &a, Point &b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}
//获取向量的长度
double get_length(vec &a)
{
    return sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
}
//获取两个向量的内积
double dot_product(vec &a, vec &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
//另外一个获取两个向量的内积
double dot_product(vector<double> &a, vector<double> &b)
{
    double x=0;
    for (int i = 0; i < a.size(); ++i)
    {
        x += a[i] * b[i];
    }
    return x;
}
//获取矩阵的行列式值
double detJ(double a[3][3])
{
    double t;
    t = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[0][2] * a[1][0] * a[2][1];
    t = t - a[0][2] * a[1][1] * a[2][0] - a[0][1] * a[1][0] * a[2][2] - a[0][0] * a[1][2] * a[2][1];
    if (t < 0)
    {
        t = -t;
    }
    else
    {
        // cout<<t<<endl;
    }
    return t;
}
//获取四面体体积，已知四个顶点的坐标
double sgn(double a)
{
    if (a > 0)
    {
        return 1;
    }
    else if (a < 0)
    {
        return -1;
    }
    return 0;
}
double get_tet_volume(vector<Point> &v1)
{
    double mat[3][3];
    mat[0][0] = v1[0].x - v1[3].x;
    mat[0][1] = v1[1].x - v1[3].x;
    mat[0][2] = v1[2].x - v1[3].x;
    mat[1][0] = v1[0].y - v1[3].y;
    mat[1][1] = v1[1].y - v1[3].y;
    mat[1][2] = v1[2].y - v1[3].y;
    mat[2][0] = v1[0].z - v1[3].z;
    mat[2][1] = v1[1].z - v1[3].z;
    mat[2][2] = v1[2].z - v1[3].z;

    double d = detJ(mat) / 6;
    return d;
};
//获取三角形面积
double get_tri_area(vector<Point> &tri)
{
    double a = get_distance(tri[0], tri[1]);
    double b = get_distance(tri[1], tri[2]);
    double c = get_distance(tri[0], tri[2]);
    double p = (a + b + c) / 2;
    return sqrt(p * (p - a) * (p - b) * (p - c));
};
//两向量相减，返回一个新的向量
vec minus(Point &a, Point &b)
{
    vec c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}
//获取两向量各个分量的差值的绝对值的最大值
double margin(vector<double> &a, vector<double> &b)
{
    double res = 0;
    for (int i = 0; i < a.size(); i++)
    {
        res = max(abs(a[i] - b[i]), res);
    }
    return res;
}
//获取高斯积分点的坐标
vector<double> GaussIntegrationPoints(double a, double b, int N)
{
    N = N - 1;
    int N1 = N + 1;
    int N2 = N + 2;
    vector<double> gauss, y, y0;
    auto *Lp = new double[N1];
    auto **L = new double *[N1];
    for (int i = 0; i < N1; i++)
    {
        L[i] = new double[N2];
    }
    for (int i = 0; i < N1; i++)
    {
        L[i][0] = 1;
    }
    for (int i = 0; i < N1; i++)
    {
        y.push_back(cos((2 * (N - i) + 1) * PI / (2 * N + 2)));
    }
    for (int i = 0; i < N1; i++)
    {
        y0.push_back(2);
    }
    while (margin(y, y0) > Eps)
    {
        for (int i = 0; i < N1; i++)
        {
            L[i][1] = y[i];
        }
        for (int j = 2; j < N1 + 1; j++)
        {
            for (int i = 0; i < N1; i++)
            {
                L[i][j] = ((2 * j - 1) * y[i] * L[i][j - 1] - (j - 1) * L[i][j - 2]) / j;
            }
        }
        for (int i = 0; i < N1; i++)
        {
            Lp[i] = (N2) * (L[i][N1 - 1] - y[i] * L[i][N2 - 1]) / (1 - y[i] * y[i]);
            y0[i] = y[i];
            y[i] = y0[i] - L[i][N2 - 1] / Lp[i];
        }
    }
    for (int i = 0; i < N1; i++)
    {
        gauss.push_back((a * (1 - y[i]) + b * (1 + y[i])) / 2);
        gauss.push_back((b - a) / ((1 - y[i] * y[i]) * Lp[i] * Lp[i]) * N2 / N1 * N2 / N1);
    }
    for (int i = 0; i < N1; i++)
    {
        delete[] L[i];
    }
    delete[] L;
    delete[] Lp;
    return gauss;
}
//获取两个向量的向量积
vec cross_product(vec &p1, vec &p2)
{
    vec p3;
    p3.x = p1.y * p2.z - p2.y * p1.z;
    p3.y = p2.x * p1.z - p1.x * p2.z;
    p3.z = p1.x * p2.y - p2.x * p1.y;
    return p3;
}
int parseLine(char *line)
{
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char *p = line;
    while (*p < '0' || *p > '9')
        p++;
    line[i - 3] = '\0';
    i = atoi(p);
    return i;
}

size_t get_host_memory()
{
    FILE *file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];
    while (fgets(line, 128, file) != nullptr)
    {
        if (strncmp(line, "VmRSS:", 6) == 0)
        {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

vec Get_gradient_2D(vector<vector<double>> &J, vector<double> &d)
{
    vec gradient;
    if (J.size() == 2)
    {
        double Jni[2][2];
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        Jni[0][0] = J[1][1] / det;
        Jni[0][1] = -J[0][1] / det;
        Jni[1][1] = J[0][0] / det;
        Jni[1][0] = -J[1][0] / det;
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                double sum = 0;
                for (int k = 0; k < 2; ++k)
                {

                    sum += Jni[i][k] * J[k][j];
                }
                // cout<<sum<<endl;
            }
        }
        gradient.x = Jni[0][0] * d[0] + Jni[0][1] * d[1];
        gradient.y = Jni[1][0] * d[0] + Jni[1][1] * d[1];
        gradient.z = 0;
    }
    else
    {
        vector<vector<double>> JT;
        JT.resize(J[0].size());
        for (int i = 0; i < J[0].size(); ++i)
        {
            JT[i].resize(J.size());
        }
        for (int i = 0; i < J.size(); ++i)
        {
            for (int j = 0; j < J[i].size(); ++j)
            {
                JT[j][i] = J[i][j];
            }
        }
        vector<vector<double>> Jmul;
        Jmul.resize(JT.size());
        for (int i = 0; i < JT.size(); ++i)
        {
            Jmul[i].resize(J.size());
        }
        for (int i = 0; i < JT.size(); ++i)
        {
            for (int j = 0; j < J[0].size(); ++j)
            {
                Jmul[i][j] = 0;
                for (int k = 0; k < JT[i].size(); ++k)
                {
                    Jmul[i][j] += JT[i][k] * J[k][j];
                }
            }
        }
        double Jni[2][2];
        double det = Jmul[0][0] * Jmul[1][1] - Jmul[0][1] * Jmul[1][0];
        Jni[0][0] = Jmul[1][1] / det;
        Jni[0][1] = -Jmul[0][1] / det;
        Jni[1][1] = Jmul[0][0] / det;
        Jni[1][0] = -Jmul[1][0] / det;
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                double sum = 0;
                for (int k = 0; k < 2; ++k)
                {

                    sum += Jni[i][k] * Jmul[k][j];
                }
                // cout<<sum<<endl;
            }
        }

        vector<vector<double>> Jfin;
        Jfin.resize(2);
        for (int i = 0; i < 2; ++i)
        {
            Jfin[i].resize(d.size(), 0);
            for (int j = 0; j < d.size(); ++j)
            {
                for (int k = 0; k < 2; ++k)
                {
                    Jfin[i][j] += Jni[i][k] * JT[k][j];
                }
            }
        }
        gradient.x = 0;
        gradient.y = 0;
        gradient.z = 0;
        for (int i = 0; i < d.size(); ++i)
        {
            gradient.x += Jfin[0][i] * d[i];
            gradient.y += Jfin[1][i] * d[i];
        }
    }

    return gradient;
}
vec Get_gradient_3D(vector<vector<double>> &J, vector<double> &d)
{
    double detJ1 = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2]) + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    vec gradient;
    vector<vector<double>> JT;
    JT.resize(J[0].size());
    for (int i = 0; i < J[0].size(); ++i)
    {
        JT[i].resize(J.size());
    }
    for (int i = 0; i < J.size(); ++i)
    {
        for (int j = 0; j < J[i].size(); ++j)
        {
            JT[j][i] = J[i][j];
        }
    }
    vector<vector<double>> Jmul;
    Jmul.resize(JT.size());
    for (int i = 0; i < JT.size(); ++i)
    {
        Jmul[i].resize(J.size());
    }
    for (int i = 0; i < JT.size(); ++i)
    {
        for (int j = 0; j < J[0].size(); ++j)
        {
            Jmul[i][j] = 0;
            for (int k = 0; k < JT[i].size(); ++k)
            {
                Jmul[i][j] += JT[i][k] * J[k][j];
            }
        }
    }
    double Jni[3][3];

    double detJ = Jmul[0][0] * (Jmul[1][1] * Jmul[2][2] - Jmul[1][2] * Jmul[2][1]) + Jmul[0][1] * (Jmul[1][2] * Jmul[2][0] - Jmul[1][0] * Jmul[2][2]) + Jmul[0][2] * (Jmul[1][0] * Jmul[2][1] - Jmul[1][1] * Jmul[2][0]);
    Jni[0][0] = 1.0 / detJ * (Jmul[1][1] * Jmul[2][2] - Jmul[1][2] * Jmul[2][1]);
    Jni[1][0] = 1.0 / detJ * (Jmul[0][2] * Jmul[2][1] - Jmul[0][1] * Jmul[2][2]);
    Jni[2][0] = 1.0 / detJ * (Jmul[0][1] * Jmul[1][2] - Jmul[0][2] * Jmul[1][1]);
    Jni[0][1] = 1.0 / detJ * (Jmul[1][2] * Jmul[2][0] - Jmul[1][0] * Jmul[2][2]);
    Jni[1][1] = 1.0 / detJ * (Jmul[0][0] * Jmul[2][2] - Jmul[0][2] * Jmul[2][0]);
    Jni[2][1] = 1.0 / detJ * (Jmul[1][0] * Jmul[0][2] - Jmul[0][0] * Jmul[1][2]);
    Jni[0][2] = 1.0 / detJ * (Jmul[1][0] * Jmul[2][1] - Jmul[1][1] * Jmul[2][0]);
    Jni[1][2] = 1.0 / detJ * (Jmul[0][1] * Jmul[2][0] - Jmul[0][0] * Jmul[2][1]);
    Jni[2][2] = 1.0 / detJ * (Jmul[0][0] * Jmul[1][1] - Jmul[1][0] * Jmul[0][1]);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double sum = 0;
            for (int k = 0; k < 3; ++k)
            {

                sum += Jni[i][k] * Jmul[k][j];
            }
            // cout<<sum<<endl;
        }
    }

    vector<vector<double>> Jfin;
    Jfin.resize(3);
    for (int i = 0; i < 3; ++i)
    {
        Jfin[i].resize(d.size(), 0);
        for (int j = 0; j < d.size(); ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                Jfin[i][j] += Jni[i][k] * JT[k][j];
            }
        }
    }
    gradient.x = 0;
    gradient.y = 0;
    gradient.z = 0;
    for (int i = 0; i < d.size(); ++i)
    {
        gradient.x += Jfin[0][i] * d[i];
        gradient.y += Jfin[1][i] * d[i];
        gradient.z += Jfin[2][i] * d[i];
    }
    return gradient;
}

void Get_inverse_2D(double **Jfin, vector<vector<double>> &J)
{
    vec gradient;
    if (J.size() == 2)
    {
        double Jni[2][2];
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        Jni[0][0] = J[1][1] / det;
        Jni[0][1] = -J[0][1] / det;
        Jni[1][1] = J[0][0] / det;
        Jni[1][0] = -J[1][0] / det;
    }
    else
    {
        vector<vector<double>> JT;
        JT.resize(J[0].size());
        for (int i = 0; i < J[0].size(); ++i)
        {
            JT[i].resize(J.size());
        }
        for (int i = 0; i < J.size(); ++i)
        {
            for (int j = 0; j < J[i].size(); ++j)
            {
                JT[j][i] = J[i][j];
            }
        }
        vector<vector<double>> Jmul;
        Jmul.resize(JT.size());
        for (int i = 0; i < JT.size(); ++i)
        {
            Jmul[i].resize(J.size());
        }
        for (int i = 0; i < JT.size(); ++i)
        {
            for (int j = 0; j < J[0].size(); ++j)
            {
                Jmul[i][j] = 0;
                for (int k = 0; k < JT[i].size(); ++k)
                {
                    Jmul[i][j] += JT[i][k] * J[k][j];
                }
            }
        }
        double Jni[2][2];
        double det = Jmul[0][0] * Jmul[1][1] - Jmul[0][1] * Jmul[1][0];
        Jni[0][0] = Jmul[1][1] / det;
        Jni[0][1] = -Jmul[0][1] / det;
        Jni[1][1] = Jmul[0][0] / det;
        Jni[1][0] = -Jmul[1][0] / det;
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                double sum = 0;
                for (int k = 0; k < 2; ++k)
                {

                    sum += Jni[i][k] * Jmul[k][j];
                }
                // cout<<sum<<endl;
            }
        }

        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < J.size(); ++j)
            {
                Jfin[i][j]=0;
            }
        }
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < J.size(); ++j)
            {
                for (int k = 0; k < 2; ++k)
                {
                    Jfin[i][j] += Jni[i][k] * JT[k][j];
                    //cout<<Jfin[i][j] <<endl;
                }
            }
        }
        //cout<<endl;
    }
}


void Get_inverse_3D(double **Jfin,vector<vector<double>> &J)
{
    double detJ1 = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2]) + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    vec gradient;
    vector<vector<double>> JT;
    JT.resize(J[0].size());
    for (int i = 0; i < J[0].size(); ++i)
    {
        JT[i].resize(J.size());
    }
    for (int i = 0; i < J.size(); ++i)
    {
        for (int j = 0; j < J[i].size(); ++j)
        {
            JT[j][i] = J[i][j];
        }
    }
    vector<vector<double>> Jmul;
    Jmul.resize(JT.size());
    for (int i = 0; i < JT.size(); ++i)
    {
        Jmul[i].resize(J.size());
    }
    for (int i = 0; i < JT.size(); ++i)
    {
        for (int j = 0; j < J[0].size(); ++j)
        {
            Jmul[i][j] = 0;
            for (int k = 0; k < JT[i].size(); ++k)
            {
                Jmul[i][j] += JT[i][k] * J[k][j];
            }
        }
    }
    double Jni[3][3];

    double detJ = Jmul[0][0] * (Jmul[1][1] * Jmul[2][2] - Jmul[1][2] * Jmul[2][1]) + Jmul[0][1] * (Jmul[1][2] * Jmul[2][0] - Jmul[1][0] * Jmul[2][2]) + Jmul[0][2] * (Jmul[1][0] * Jmul[2][1] - Jmul[1][1] * Jmul[2][0]);
    Jni[0][0] = 1.0 / detJ * (Jmul[1][1] * Jmul[2][2] - Jmul[1][2] * Jmul[2][1]);
    Jni[1][0] = 1.0 / detJ * (Jmul[0][2] * Jmul[2][1] - Jmul[0][1] * Jmul[2][2]);
    Jni[2][0] = 1.0 / detJ * (Jmul[0][1] * Jmul[1][2] - Jmul[0][2] * Jmul[1][1]);
    Jni[0][1] = 1.0 / detJ * (Jmul[1][2] * Jmul[2][0] - Jmul[1][0] * Jmul[2][2]);
    Jni[1][1] = 1.0 / detJ * (Jmul[0][0] * Jmul[2][2] - Jmul[0][2] * Jmul[2][0]);
    Jni[2][1] = 1.0 / detJ * (Jmul[1][0] * Jmul[0][2] - Jmul[0][0] * Jmul[1][2]);
    Jni[0][2] = 1.0 / detJ * (Jmul[1][0] * Jmul[2][1] - Jmul[1][1] * Jmul[2][0]);
    Jni[1][2] = 1.0 / detJ * (Jmul[0][1] * Jmul[2][0] - Jmul[0][0] * Jmul[2][1]);
    Jni[2][2] = 1.0 / detJ * (Jmul[0][0] * Jmul[1][1] - Jmul[1][0] * Jmul[0][1]);

    for (int i = 0; i < 3; ++i)
    {

        for (int j = 0; j < J.size(); ++j)
        {
            Jfin[i][j]=0;
            for (int k = 0; k < 3; ++k)
            {
                Jfin[i][j] += Jni[i][k] * JT[k][j];
            }
        }
    }


}
#ifndef USE_GPU
// #else


#endif