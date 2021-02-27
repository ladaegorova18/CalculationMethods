#include "Lab1.2.h"

const double epsilon = 0.00001;

int main()
{
    setlocale(LC_ALL, "Russian");

    cout << "Решение линейной алгебраической системы, нахождение числа обусловленности\n";
    cout << "Вариант 3\n";
    vector<vector<double>> A = { {8.29381, 0.995516, -0.560617, 0.766522}, {0.995516, 6.298198, 0.595772, 3.84422}, {-0.560617, 0.595772, 4.997407, 5.239231} };
    cout << "A:" << endl;
    printMatrix(A);

    vector<vector<double>> C = A;
    C[0][0] = A[0][0] * pow(10, -8);
    cout << "C:" << endl;
    printMatrix(C);

    cout << "Метод Гаусса с выбором главного элемента для матриц A и С:\n";
    
    vector<double> resultA = gaussWithMainElement(A);
    cout << "Матрица A:" << endl;
    printResult(resultA, "X");
    printPrecision(A, resultA);

    vector<double> resultC = gaussWithMainElement(C);
    cout << "Матрица C:" << endl;
    printResult(resultC, "X");
    printPrecision(C, resultC);

    cout << "Метод LU-разложения для нахождения определителя матрицы А:\n";
    double det = LUdecomposition(A.size(), A);
    cout << "det A: " << setprecision(15) << det << endl;

    cout << "Обратная матрица для A:\n";
    auto reverseA = reverseMatrix(A, A.size());
    printMatrix(reverseA);
}

void printPrecision(vector<vector<double>> A, vector<double> X)
{
    vector<double> R;
    size_t n = A.size();
    vector<double> Ax(n);
    vector<double> b(n);

    for (size_t i = 0; i < A.size(); ++i)
    {
        double sum = 0;
        for (size_t j = 0; j < A.size(); ++j)
        {
            sum += A[i][j] * X[j];
        }
        Ax[i] = sum;
        b[i] = A[i][n];
    }
    for (size_t j = 0; j < n; ++j)
    {
        R.push_back(Ax[j] - b[j]);
    }
    cout << "Невязка:" << endl;
    printResult(R, "R");
}

void printResult(vector<double> X, string name)
{
    cout << name << ":" << endl;
    for (size_t i = 0; i < X.size(); ++i)
    {
        cout << setprecision(15) << X[i] << endl;
    }
    cout << endl;
}

void printMatrix(vector<vector<double>> A)
{
    for (int i = 0; i < A.size(); i++)
        for (int j = 0; j < A[i].size(); j++)
        {
            bool isRect = j == A[i].size() - 1;
            cout << " |"[isRect && A.size() + 1 == A[0].size()];
            cout << setprecision(15) << setw(20) << A[i][j] << " \n"[isRect];
        }
    cout << endl;
}

vector<double> GaussMethod(vector<vector<double>> A)
{
    A = ForvardMove(A, false);
    auto X = ReverseMove(A);
    return X;
}

vector<double> gaussWithMainElement(vector<vector<double>> A)
{
    A = ForvardMove(A, true);
    auto X = ReverseMove(A);
    return X;
}

vector<vector<double>> ForvardMove(vector<vector<double>> A, bool chooseMainElement)
{
    size_t n = A.size();
    for (size_t k = 0; k < n; ++k)
    {
        double tmp = abs(A[k][k]);
        if (tmp < epsilon)
        {
            if (chooseMainElement)
            {
                A = chooseMain(A, k);
            }
            else
            {
                cout << "Слишком маленький ведущий элемент " << setprecision(15) 
                    << A[k][k] << endl;
            }
        }
        for (size_t j = k; j < n + 1; ++j)
        {
            A[k][j] = A[k][j] / tmp;
        }
        for (size_t i = k + 1; i < n; ++i)
        {
            tmp = A[i][k];
            for (size_t j = k; j < n + 1; ++j)
            {
                A[i][j] = A[i][j] - A[k][j] * tmp;
            }
        }
    }
    return A;
}

vector<double> ReverseMove(vector<vector<double>> A)
{
    size_t n = A[0].size() - 1;
    vector<double> X(A.size());
    for (int i = n - 1; i > -1; i--)
    {
        double sum = 0;
        for (size_t j = i + 1; j < n; ++j)
        {
            sum += A[i][j] * X[j];
        }
        X[i] = A[i][n] - sum;
    }
    return X;
}

double LUdecomposition(size_t n, vector<vector<double>> A)
{
    vector<vector<double>> L(n, vector<double>(n));
    vector<vector<double>> U(n, vector<double>(n));

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i; j < n; ++j)
        {
            double sum = 0;
            for (size_t k = 0; k < i; ++k)
            {
                sum += L[j][k] * U[k][i];
            }
            L[j][i] = A[j][i] - sum;

            sum = 0;
            for (size_t k = 0; k < i; ++k)
            {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - sum) / L[i][i];
        }
    }
    double det = 1;
    for (size_t i = 0; i < n; ++i)
    {
        det *= L[i][i];
    }
    return det;
}

vector<vector<double>> Jordan(vector<vector<double>> A)
{
    size_t n = A.size();
    for (size_t k = 0; k < n; ++k)
    {
        double tmp = abs(A[k][k]);
        if (tmp < epsilon)
        {
            cout << "Слишком маленький ведущий элемент " << setprecision(15)
                << A[k][k] << endl;
        }
        for (size_t j = k; j < A[0].size(); ++j)
        {
            A[k][j] = A[k][j] / tmp;
        }
        for (size_t i = 0; i < n; ++i)
        {
            tmp = A[i][k];
            if (i != k)
            {
                for (size_t j = k; j < A[0].size(); ++j)
                {
                    A[i][j] = A[i][j] - A[k][j] * tmp;
                }
            }
        }
    }
    return A;
}

vector<vector<double>> reverseMatrix(vector<vector<double>> A, size_t n)
{
    vector<vector<double>> reverseA (n, vector<double>((A[0].size() - 1) * 2));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            reverseA[i][j] = A[i][j];
            if (i == j)
            {
                reverseA[i][j + n] = 1;
            }
        }
    }
    reverseA = Jordan(reverseA);
    vector<vector<double>> result(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            result[i][j] = reverseA[i][j + n];
        }
    }

    return result;
}

vector<vector<double>> chooseMain(vector<vector<double>> A, int k)
{
    int p = 0;
    double max = abs(A[k][k]);
    for (size_t i = k; i < A.size(); ++i)
    {
        if (abs(A[i][k]) > max)
        {
            max = abs(A[i][k]);
            p = i;
        }
    }
    auto tmp = A[k];
    A[k] = A[p];
    A[p] = tmp;
    return A;
}