

#ifndef LABWORK1_MATRIX_WORK_H
#define LABWORK1_MATRIX_WORK_H

#include <iostream>
#include <fstream>
#include "vector"
#include "cmath"

using namespace std;

template <class T>
class matrix {
private:
    int k, l;
    T** m_data;
public:
    matrix(int n)
    {
        k = n;
        l = n;
        m_data = new T * [k];
        for (int i = 0; i < k; ++i) {
            m_data[i] = new T[l];
        }
    }
    matrix(int n, int m) {
        k = n;
        l = m;
        m_data = new T * [k];
        for (int i = 0; i < k; ++i) {
            m_data[i] = new T[l];
        }
    }
    matrix(int param, char* link) {
        ifstream fin(link);
        k = 0;
        l = 0;
        fin >> k;
        if (param == 2) {
            fin >> l;
        }
        else {
            l = k;
        }
        m_data = new T * [k];
        for (int i = 0; i < k; ++i) {
            m_data[i] = new T[l];
            for (int j = 0; j < l; ++j) {
                fin >> m_data[i][j];
            }
        }
    }
    matrix(matrix<T> const& input_m) {
        k = input_m.k;
        l = input_m.l;
        m_data = new T * [k];
        for (int i = 0; i < k; ++i) {
            m_data[i] = new T[l];
            for (int j = 0; j < l; ++j) {
                m_data[i][j] = input_m[i][j];
            }
        }
    }
    T*& operator[](int n) const {
        return m_data[n];
    }
    void init_zero() {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < l; ++j) {
                m_data[i][j] = 0;
            }
        }
    }
    void init_unit() {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < l; ++j) {
                if (i == j) {
                    m_data[i][j] = 1;
                }
                else {
                    m_data[i][j] = 0;
                }
            }
        }
    }
    void console_output() const {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < l; ++j) {
                cout << m_data[i][j] << " ";
            }
            cout << std::endl;
        }
    }
    void output_f(ofstream& out) const {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < l; ++j) {
                out << m_data[i][j] << " ";
            }
            out << std::endl;
        }
    }
    int lin() const {
        return k;
    }

    int col() const {
        return l;
    }
    matrix transpose() {
        matrix<T> output(l, k);
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < k; ++j) {
                output[i][j] = m_data[j][i];
            }
        }
        return output;
    }
    matrix operator*(matrix<T> const& right) {
        matrix<T> result(k, right.col());
        if (l != right.lin()) {
            throw "Ошибка размеров при матричном умножении";
        }
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < right.col(); ++j) {
                result[i][j] = 0;
                for (int m = 0; m < l; ++m) {
                    result[i][j] += m_data[i][m] * right[m][j];
                }
            }
        }
        return result;
    }
    matrix operator=(matrix<T> const& right) {
        for (int i = 0; i < k; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
        k = right.k;
        l = right.l;
        m_data = new T * [k];
        for (int i = 0; i < k; ++i) {
            m_data[i] = new T[l];
            for (int j = 0; j < l; ++j) {
                m_data[i][j] = right[i][j];
            }
        }
        return *this;
    }
    T det() const {
        if (k != l) {
            throw "Попытка найти определитель неквадратной матрицы";
        }
        if (k == 1) return m_data[0][0];
        T res = 0;
        T sign = 1;
        for (int i = 0; i < k; ++i) {
            matrix<T> h(k - 1, k - 1);
            for (int j = 1; j < k; ++j) {
                for (int m = 0; m < i; ++m) {
                    h[j - 1][m] = m_data[j][m];
                }
                for (int m = i + 1; m < k; ++m) {
                    h[j - 1][m - 1] = m_data[j][m];
                }
            }
            res += m_data[0][i] * sign * h.det();
            sign = sign * (-1);
        }
        return res;
    }
    vector<T> operator*(vector<T> const& right) const {
        if (l != right.size()) {
            throw "Ошибка размеров при умножении вектора на матрицу";
        }
        vector<T> result(right.size());
        for (int i = 0; i < k; ++i) {
            result[i] = 0;
            for (int j = 0; j < l; ++j) {
                result[i] += m_data[i][j] * right[j];
            }
        }
        return result;
    }
    T abs_1() {
        T result = 0;
        for (int j = 0; j < l; ++j) {
            T cand = 0;
            for (int i = 0; i < k; ++i) {
                cand += abs(m_data[i][j]);
            }
            if (cand > result) result = cand;
        }
        return result;
    }
    T abs_inf() {
        T result = 0;
        for (int i = 0; i < l; ++i) {
            T cand = 0;
            for (int j = 0; j < k; ++j) {
                cand += abs(m_data[i][j]);
            }
            if (cand > result) result = cand;
        }
        return result;
    }
    matrix<T> inverse() {
        if (k != l) {
            throw "Поиск обратной у неквадратной матрицы";
        }
        if (this->det() == 0) {
            throw "Невозможно найти обратную матрицу";
        }
        matrix<T> result(k);
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < l; ++j) {
                matrix<T> h(k - 1, k - 1);
                for (int u = 0; u < i; ++u) {
                    for (int m = 0; m < j; ++m) {
                        h[u][m] = m_data[u][m];
                    }
                    for (int m = j + 1; m < l; ++m) {
                        h[u][m - 1] = m_data[u][m];
                    }
                }
                for (int u = i + 1; u < k; ++u) {
                    for (int m = 0; m < j; ++m) {
                        h[u - 1][m] = m_data[u][m];
                    }
                    for (int m = j + 1; m < l; ++m) {
                        h[u - 1][m - 1] = m_data[u][m];
                    }
                }
                result[j][i] = pow(-1, i + j) * h.det() / this->det();
            }
        }
        return result;
    }
    ~matrix() {
        for (int i = 0; i < l; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }
};


#endif //LABWORK1_MATRIX_WORK_H
