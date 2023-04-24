// Demyan Zverev
// d.zverev@innopolis.university
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

class Matrix {
public:
    int row, column;
    vector<vector<double>>matrix;

    Matrix() {
        row = 0, column = 0;
    }

    Matrix(int n, int m) {
        row = n;
        column = m;
        matrix.resize(n, vector<double>(m, 0));
    }

    friend istream& operator >> (istream& input, Matrix& cur_matrix) {
        input >> cur_matrix.row;
        cur_matrix.row = 2;
        cur_matrix.matrix.resize(cur_matrix.row, vector<double>(cur_matrix.column));
        for (int i = 0; i < cur_matrix.row; i++)
            for (int j = 0; j < cur_matrix.column; j++)
                input >> cur_matrix.matrix[i][j];
        return input;
    }

    friend ostream& operator << (ostream& output, Matrix& cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++) {
            for (int j = 0; j < cur_matrix.column; j++)
                output << fixed << setprecision(4) << cur_matrix.matrix[i][j] << " ";
            output << endl;
        }
        return output;
    }

    Matrix operator = (Matrix cur_matrix) {
        column = cur_matrix.column;
        row = cur_matrix.row;
        for (int i = 0; i < column; i++)
            for (int j = 0; j < row; j++)
                matrix[i][j] = cur_matrix.matrix[i][j];
        return *this;
    }

    friend Matrix operator + (Matrix first, Matrix second) {
        Matrix ans(first.row, first.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = first.matrix[i][j] + second.matrix[i][j];
        return ans;
    }

    friend Matrix operator - (Matrix first, Matrix second) {
        Matrix ans(first.row, first.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = first.matrix[i][j] - second.matrix[i][j];
        return ans;
    }

    friend Matrix operator * (Matrix first, Matrix second) {
        Matrix ans(first.row, second.column);
        for (int i = 0; i < first.row; i++)
            for (int j = 0; j < second.column; j++)
                for (int l = 0; l < first.matrix[i].size(); l++)
                    ans.matrix[i][j] += first.matrix[i][l] * second.matrix[l][j];
        return ans;
    }

    Matrix transpose() {
        Matrix ans(this->column, this->row);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < this->column; j++)
                ans.matrix[j][i] = this->matrix[i][j];
        return ans;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int x) {
        row = column = x;
        matrix.resize(x, vector<double>(x));
    }
    SquareMatrix() {
        column = row = 0;
    }

    friend istream& operator >> (istream& input, SquareMatrix& cur_matrix) {
        input >> cur_matrix.row;
        cur_matrix.column = cur_matrix.row;
        cur_matrix.matrix.resize(cur_matrix.row, vector<double>(cur_matrix.row));
        for (int i = 0; i < cur_matrix.row; i++)
            for (int j = 0; j < cur_matrix.row; j++)
                input >> cur_matrix.matrix[i][j];
        return input;
    }

    friend ostream& operator << (ostream& output, SquareMatrix& cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++) {
            for (int j = 0; j < cur_matrix.row; j++)
                output << fixed << setprecision(2) << (cur_matrix.matrix[i][j] == 0 ? 0 : cur_matrix.matrix[i][j]) << " ";
            output << endl;
        }
        return output;
    }

    friend SquareMatrix operator +(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir + *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    friend SquareMatrix operator-(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir - *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    friend SquareMatrix operator*(SquareMatrix first, SquareMatrix second) {
        Matrix* fir = &first, * sec = &second;
        Matrix ans = *fir * *sec;
        SquareMatrix* an = (SquareMatrix*)&ans;
        return *an;
    }

    SquareMatrix transpose() {
        Matrix fir = (*this);
        fir = (fir).transpose();
        SquareMatrix* ans = (SquareMatrix*)&fir;
        return *ans;
    }

    int max_in_column(int i) {
        int ans = i;
        for (int j = i; j < matrix.size(); j++)
            if (abs(matrix[ans][i]) < abs(matrix[j][i]))
                ans = j;
        return ans;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() {
        row = column = 0;
    }
    IdentityMatrix(int n) {
        row = column = n;
        matrix.resize(n, vector<double>(n, 0));
        for (int i = 0; i < n; i++)
            matrix[i][i] = 1;
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix() {}
    EliminationMatrix(SquareMatrix x, int i, int j) : IdentityMatrix(x.column) {
        matrix[i - 1][j - 1] = (-x.matrix[i - 1][j - 1]);
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(SquareMatrix x, int i, int j) : IdentityMatrix(x.column) {
        swap(matrix[i - 1], matrix[j - 1]);
    }
};

Matrix Inverse(SquareMatrix a) {
    int n = a.column;
    SquareMatrix id = IdentityMatrix(n);
    int k = 1;
    for (int l = 0; l < n - 1; l++) {
        int cur = a.max_in_column(l);
        if (cur != l) {
            PermutationMatrix p(a, cur + 1, l + 1);
            a = p * a;
            id = p * id;
        }
        for (int i = l + 1; i < n; i++) {
            if (a.matrix[i][l] != 0) {
                double now = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j++) {
                    a.matrix[i][j] -= now * a.matrix[l][j];
                    id.matrix[i][j] -= now * id.matrix[l][j];
                }
            }
        }
    }
    for (int l = n - 1; l >= 0; l--) {
        for (int i = l - 1; i >= 0; i--) {
            if (a.matrix[i][l] != 0) {
                double now = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j++) {
                    a.matrix[i][j] -= now * a.matrix[l][j];
                    id.matrix[i][j] -= now * id.matrix[l][j];
                }
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double now = a.matrix[i][i];
        for (int j = 0; j < n; j++)
            id.matrix[i][j] /= now;
        a.matrix[i][i] = 1;
    }
    return id;
}

double pow(double x, int p) {
    double ans = 1;
    for (int i = 0; i < p; i++)
        ans *= x;
    return ans;
}

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#endif

int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#endif
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);

    int k; cin >> k;
    Matrix b(k, 1);
    vector<double>t(k);
    for (int i = 0; i < k; i++) {
        cin >> t[i] >> b.matrix[i][0];
    }
    int n; cin >> n;
    Matrix a(k, n + 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n + 1; j++)
            if (j >= 1)
                a.matrix[i][j] = pow(t[i], j);
            else
                a.matrix[i][j] = 1;
    }
    cout << "A:\n" << a;
    Matrix s = a.transpose();
    Matrix c = s * a;
    cout << "A_T*A:\n" << c;
    SquareMatrix* w = (SquareMatrix*)&c;
    c = Inverse(*w);
    cout << "(A_T*A)^-1:\n" << c;
    Matrix ww = s * b;
    cout << "A_T*b:\n" << ww;
    cout << "x~:\n";
    Matrix ans = c * ww;
    cout << ans;
    fprintf(pipe, "plot [-30 : 30] [-30 : 30] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", ans.matrix[3][0], ans.matrix[2][0], ans.matrix[1][0], ans.matrix[0][0]);
    fprintf(pipe, "plot '-' w p ls 1, '-' w p ls 2, '-' w p ls 3, '-' w p ls 4, '-' w p ls 5, '-' w p ls 6, '-' w p ls 7, '-' w p ls 8, '-' w p ls 9, '-' w p ls 10\n");
    for (int i = 0; i < k; i++) {
        fprintf(pipe, "%f\t%f\n", t[i], b.matrix[i][0]);
    }
#ifdef WIN32
    _pclose(pipe);
#else
    _pclose(pipe);
#endif
}