#include "Lab1.h"

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Решение линейной алгебраической системы, нахождение числа обусловленности\n";
    cout << "Вариант 7\n";
    vector<vector<double>> A = {{-401.43, 200.19}, {1201.14, -601.62}};
    vector<vector<double>> b = { {200}, {-600} };
    vector<vector<double>> b_ = { {199}, {-601} }; /// b 
    
    cout << "A:" << endl;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "\n";

    double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0]; /// count reverse to A matrix
    vector<vector<double>> M = { {A[1][1], -A[0][1]}, { -A[1][0], A[0][0]} };
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            M[i][j] = M[i][j] * 1.0 / detA; 
        }
    }
    
    cout << "Решение системы с точной правой частью x:" << endl;
    vector<vector<double>> x = countX(M, b);

    cout << "Решение системы с возмущённой правой частью x_:" << endl;
    vector<vector<double>> x_ = countX(M, b_);/// M = A^(-1) ; Ax = b ; x = A^(-1) * B

    double cond = norm(M) * norm(A); /// count cond(A)
    cout << "Число обусловленности cond: " << cond << endl;

    vector<vector<double>> x_x = { {x[0][0] - x_[0][0]}, {x[1][0] - x_[1][0]} };
    double deltaX = norm(x_x) / norm(x);

    /// <summary>
    /// count norm of matrices |A - A_|, |b - b_|
    /// </summary>
    cout << "Фактическая относительная погрешность deltaX: " << deltaX << endl;
    vector<vector<double>> dA = { {0, 0}, {0, 0} };
    vector<vector<double>> db = countDelta(b, b_);

    double R = cond * (norm(db) / norm(b) + norm(dA) / norm(A))
        / (1 - cond * norm(dA) / norm(A));
    cout << "Оценка погрешности R: " << R << endl;
}

vector<vector<double>> countDelta (vector<vector<double>> A, vector<vector<double>> A_)
{
    vector<vector<double>> dA = A_;
    for (int i = 0; i < dA.size(); ++i)
    {
        for (int j = 0; j < dA[0].size(); ++j)
        {
            dA[i][j] = dA[i][j] - A[i][j];
        }
    }
    return dA;
}

/// <summary>
/// multiply matrix to vector
/// </summary>
vector<vector<double>> countX(vector<vector<double>> M, vector<vector<double>> b) 
{
    vector<vector<double>> x = { {M[0][0] * b[0][0] + M[0][1] * b[1][0]}, {M[1][0] * b[0][0] + M[1][1] * b[1][0]} };
    cout << x[0][0] << endl;
    cout << x[1][0] << endl;
    return x;
}

/// <summary>
/// count matrix norm (p = 2)
/// </summary>
double norm(vector<vector<double>> A)
{
    double sum = 0;
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A[0].size(); ++j)
        {
            sum += pow(A[i][j], 2);
        }
    }
    return pow(sum, 0.5);
}