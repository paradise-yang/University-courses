#include <vector>
using namespace std;

typedef unsigned index;
typedef double T;
typedef int Status;
typedef std::vector<std::vector<T>> matrix;

//相关运算
//向量点乘向量
T vectorMultiplyvector(const vector<T> &u, const vector<T> &v)
{
    index n1 = u.size();
    index n2 = v.size();
    T s = 0;
    if (n1 != n2) {
        printf("vectorMultiplyvector:向量点乘尺寸异常\n");
        return s;
    }
    for (index i = 0; i < n1; i++) {
        s = s + u[i] * v[i];
    }
    return s;
}

//向量加法重载
vector<T> operator+(const vector<T> &u, const vector<T> &v)
{
    index n1 = u.size();
    index n2 = v.size();
    if (n1 != n2) {
        printf("operator+:向量相加尺寸不同\n");
        return u;
    }

    vector<T> w(n1);
    for (index i = 0; i < n1; i++)
        w[i] = u[i] + v[i];

    return w;
}

//向量减法重载
vector<T> operator-(const vector<T> &u, const vector<T> &v)
{
    index n1 = u.size();
    index n2 = v.size();
    if (n1 != n2) {
        printf("operator-:向量相减尺寸不同\n");
        return u;
    }

    vector<T> w(n1);
    for (index i = 0; i < n1; i++)
        w[i] = u[i] - v[i];

    return w;
}

//向量数乘重载
vector<T> operator*(T x, const vector<T> &u)
{
    index n = u.size();
    vector<T> w(n);
    for (index i = 0; i < n; i++) {
        w[i] = x * u[i];
    }
    return w;
}


vector<T> quick_multiply_3_line(const matrix &A, const vector<T> &b)
{
    index n = A.size();
    if (n != A[0].size() || n != b.size()) {
        printf("error\n");
        return b;
    }
    vector<T> result(n,0);

    for (index i = 0; i < n - 1; i++) {
        result[i] += (A[i][i] * b[i]);
        result[i] += (A[i][i + 1] * b[i + 1]);
        result[i + 1] += (A[i + 1][i] * b[i]);
    }
    result[n - 1] += (A[n - 1][n - 1] * b[n - 1]);

    return result;
}

vector<T> Conjugate_Gradient(const matrix &A, const vector<T> &b, const vector<T> &start, T precision, index num_max)
{
    if (A.size() != A[0].size() || A[0].size() != b.size() || b.size() != start.size()) {
        printf("共轭梯度法:矩阵与向量尺寸有问题\n");
        return b;
    }

    vector<T> x = start;
    vector<T> r(b), p(b), w(b), temp(b);

    index k = 0;
    index n = b.size();
    index k_max = 10 * n;

    if (num_max > 0 && k_max < num_max) k_max = num_max; //妥善处理最大迭代次数的参数

    T rou = 0, pre_rou = 0;
    T beta = 0, alpha = 0;
    T temp_T = 0;
    T epsilon = 1e-8;
    if (precision > 0) epsilon = precision; //妥善处理精度参数

    T up = epsilon * epsilon * (vectorMultiplyvector(b,b)); //误差上限平方

    k = 0;
    // temp = A * x;
    temp = quick_multiply_3_line(A, x);
    r = b - temp;
    rou = vectorMultiplyvector(r,r);
    while (rou > up && k < k_max) {
        k = k + 1;
        if (k == 1) {
            p = r;
        }
        else {
            beta = rou / pre_rou;
            p = r + beta * p;
        }
        // w = A * p; //每次迭代只有这里有一次矩阵乘以向量，是最大的计算量
        w = quick_multiply_3_line(A, p);
        temp_T = vectorMultiplyvector(p,w);
        alpha = rou / temp_T;
        x = x + alpha * p;
        r = r - alpha * w;

        pre_rou = rou;
        rou = vectorMultiplyvector(r,r);
    }
    if (k >= k_max) {
        printf("共轭梯度法:超过最大迭代次数%d\n", k_max);
    }

    return x;
}
