#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class ColumnVector {
private:
    int size;
    vector<double> data;
public:
    ColumnVector(int m) {
        size = m;
        data.resize(m);
    }

    double& operator [] (int index) {
        return data[index];
    }

    const double& operator [] (int index) const {
        return data[index];
    }

    friend istream &operator >> (istream& input, ColumnVector& vector) {
        for (int i = 0; i < vector.size; i++) {
            input >> vector.data[i];
        }
        return input;
    }

    friend ostream& operator << (ostream& output, ColumnVector& vector) {
        output << fixed << setprecision(4);
        for (int i = 0; i < vector.size; i++) {
            output << (abs(vector[i]) < 1e-10 ? 0.00 : vector[i]) << endl;
        }
        return output;
    }

    ColumnVector operator + (ColumnVector& other_vector) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = data[i] + other_vector[i];
        }
        return result;
    }

    ColumnVector operator - (ColumnVector& other_vector) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = data[i] - other_vector[i];
        }
        return result;
    }

    double operator * (ColumnVector& other_vector) {
        double result = 0;
        for (int i = 0; i < size; i++) {
            result += data[i] * other_vector[i];
        }
        return result;
    }

    ColumnVector operator * (double scalar) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    double norm() const {
        double result = 0;
        for (int i = 0; i < size; i++) {
            result += data[i] * data[i];
        }
        return sqrt(result);
    }

};

class Matrix {
protected:
    int n, m; // n - rows, m - columns
    vector<vector<double>> matrix;
public:
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        matrix.resize(n);
        for (int i = 0; i < n; i++) {
            matrix[i].resize(m);
        }
    }

    int& getN() {
        return n;
    }

    int& getM() {
        return m;
    }

    vector<double>& operator [] (int index) {
        return matrix[index];
    }

    const vector<double>& operator [] (int index) const {
        return matrix[index];
    }

    friend istream &operator >> (istream& input, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                input >> matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator << (ostream& output, Matrix &matrix) {
        output << fixed << setprecision(4);
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                if (j == matrix.m - 1) {
                    output << (abs(matrix[i][j]) < 1e-10 ? 0.00 : matrix[i][j]);
                } else {
                    output << (abs(matrix[i][j]) < 1e-10 ? 0.00 : matrix[i][j]) << " ";
                }
            }
            output << endl;
        }
        return output;
    }

    Matrix operator = (const Matrix &matrix) {
        if (this != &matrix) {
            this->n = matrix.n;
            this->m = matrix.m;
            this->matrix = matrix.matrix;
        }
        return *this;
    }

    Matrix operator + (Matrix &second) {
        if (n != second.n || m != second.m) {
            return Matrix(0, 0);
        }
        Matrix result(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = matrix[i][j] + second[i][j];
            }
        }
        return result;
    }

    Matrix operator - (Matrix &second) {
        if (n != second.n || m != second.m) {
            return Matrix(0, 0);
        }
        Matrix result(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = matrix[i][j] - second[i][j];
            }
        }
        return result;
    }

    Matrix operator * (Matrix &second) {
        if (m != second.n) {
            return Matrix(0, 0);
        }
        Matrix result(n, second.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < second.m; j++) {
                double sum = 0;
                for (int k = 0; k < m; k++) {
                    sum += matrix[i][k] * second[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    ColumnVector operator * (ColumnVector& vector) {
        ColumnVector result(n);
        for (int i = 0; i < n; i++) {
            result[i] = 0;
            for (int j = 0; j < m; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    Matrix transpose() {
        Matrix result(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {}
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1;
        }
    }
};

class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(SquareMatrix& A, int n, int j, int i) : SquareMatrix(n) {
        double ratio = A[j][i]/ A[i][i];
        for (int k = 0; k < n; k++) {
            matrix[k][k] = 1;
        }
        matrix[j][i] = - ratio;
    }
};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n, int i, int j) : SquareMatrix(n) {
        for (int k = 0; k < n; k++) {
            if (k == i) {
                matrix[k][j] = 1;
            } else if (k == j) {
                matrix[k][i] = 1;
            } else {
                matrix[k][k] = 1;
            }
        }
    }
};

class AugmentedMatrix : public Matrix {
public:
    AugmentedMatrix(Matrix& A, Matrix& I, int n, int m) : Matrix(n, m) {
        for (int i = 0; i < n; i++) {
            matrix[i] = A[i];
            matrix[i].insert(matrix[i].end(), I[i].begin(), I[i].end());
        }
    }
};

class MatrixA : public  Matrix {
public:
    MatrixA(int n, int m, ColumnVector& t) : Matrix(n, m) {
        for (int i = 0; i < n; i++) {
            matrix[i][0] = 1;
            for (int j = 1; j < m; j++) {
                matrix[i][j] = pow(t[i], j);
            }
        }
    }
};

Matrix findInverseMatrix(SquareMatrix& A) {
    int n = A.getN();
    IdentityMatrix I(n);
    Matrix* pA = &A;
    Matrix M = *pA;
    Matrix* pI = &I;
    // perform permutation as a first step if necessary
    int pivot_row = 0;
    for (int j = 1; j < n; j++) {
        if (abs(M[j][0]) > abs(M[pivot_row][0])) {
            pivot_row = j;
        }
    }
    if (pivot_row != 0) {
        PermutationMatrix P(n, 0, pivot_row);
        M = P * M;
        *pI = P * *pI;
    }
    for (int i = 0; i < n; i++) {
        SquareMatrix N = *(SquareMatrix*)(&M);
        // elimination step
        for (int j = i + 1; j < n; j++) {
            if (abs(N[j][i]) != 0) {
                EliminationMatrix E(N, n, j, i);
                M = E * M;
                *pI = E * *pI;
            }
        }
        for (int k = i; k < n; k++) {
            // find pivot row
            pivot_row = k;
            for (int j = k + 1; j < n; j++) {
                if (abs(M[j][k]) > abs(M[pivot_row][k])) {
                    pivot_row = j;
                }
            }
            // perform permutation if necessary
            if (pivot_row != k) {
                PermutationMatrix P(n, k, pivot_row);
                M = P * M;
                *pI = P * *pI;
            }
        }
    }
    SquareMatrix N = *(SquareMatrix*)(&M);
    // elimination step
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (abs(N[j][i]) != 0) {
                EliminationMatrix E(N, n, j, i);
                M = E * M;
                *pI = E * *pI;
            }
        }
    }
    // diagonal normalization
    IdentityMatrix I1 = *(IdentityMatrix*)(&(*pI));
    for (int i = 0; i < n; i++) {
        double diagonal = M[i][i];
        for (int j = 0; j < n; j++) {
            M[i][j] /= diagonal;
            I1[i][j] /= diagonal;
        }
    }
    *pI = I1;
    return *pI;
}

void approximation(ColumnVector& t, ColumnVector& b, int m, int n) {

    // print matrix A
    MatrixA A(m, n + 1, t);
    cout << "A:" << endl;
    cout << A;

    Matrix* pA = &A;

    // compute A_T
    Matrix aT = pA->transpose();

    // compute and print A_T*A
    Matrix RES1 = aT * *pA;
    cout << "A_T*A:" << endl;
    cout << RES1;

    // compute and print (A_T*A)^-1
    SquareMatrix tempRES2 = *(SquareMatrix*)(&RES1);
    Matrix RES2 = findInverseMatrix(tempRES2);
    cout << "(A_T*A)^-1:" << endl;
    cout << RES2;

    // compute and print A_T*b
    ColumnVector RES3 = aT * b;
    cout << "A_T*b:" << endl;
    cout << RES3;

    // compute and print x~ = (A_T*A)^-1 * A_T*b
    ColumnVector x = RES2 * RES3;
    cout << "x~:" << endl;
    cout << x;
}

int main() {

    int m;
    cin >> m;

    ColumnVector t(m);
    ColumnVector b(m);

    for (int i = 0; i < m; i++) {
        cin >> t[i] >> b[i];
    }

    int n;
    cin >> n;

    approximation(t, b, m, n);

    return 0;
}
