
#include "Header1.h"
#define not_square 2
#define square 1
 //Данные макросы передаются при считывании матрицы из файла
 //Параметр 1 означает, что матрица квадратная, параметр 2 - что не квадратная
#include "vector"
#include <iostream>
#include <fstream>
#include "cmath"
#define link_data "/home/qw/Рабочий стол/laba1/"
// Данный макрос используется для указания пути к файлам с исходными данными

using namespace std;

template <typename float_type>
vector<float_type> read_vector(char* link) {
    ifstream in(link);
    int n;
    in >> n;
    vector<float_type> read_vec(n);
    for (int i = 0; i < n; ++i) {
        double a;
        in >> a;
        read_vec [i] = a;
    }
    return read_vec;
}

template <typename T>
int tVectorOut(vector<T> const& V, ofstream& out)
//Вывод вектора
{
    for (int i = 0; i < V.size(); ++i) {
        out << V[i] << " ";
    }
    out << endl;
    return 0;
}

template <typename float_type>
float_type abs_1(vector<float_type> vec)
//Норма вектора ||*||_1
{
    float_type result = 0;
    for (int i = 0; i < vec.size(); ++i) {
        result += abs(vec[i]);
    }
    return result;
}

template <typename float_type>
float_type abs_inf(vector<float_type> vec)
//Норма вектора ||*||_{/inf}
{
    float_type result = 0;
    for (int i = 0; i < vec.size(); ++i) {
        if (abs(vec[i]) > result) result = abs(vec[i]);
    }
    return result;
}

template <typename float_type>
vector<float_type> gauss_method(matrix<float_type> const& A, vector<float_type> const& B, ofstream& out)
//Метод Гаусса.
{
    out << "Метод Гаусса" << endl;
    matrix<float_type> AH(A);
    vector<float_type> BH(B);
    for (int i = 0; i < AH.lin() - 1; ++i) {
        float_type max = AH[i][i];
        int imax = i;
        for (int j = i + 1; j < AH.lin(); ++j) {
            if (AH[j][i] > max) {
                max = AH[j][i];
                imax = j;
            }
        }
        float_type* c = AH[i];
        AH[i] = AH[imax];
        AH[imax] = c;
        c = nullptr;
        if (abs(AH[i][i]) < 1e-15) {
            throw "Решение СЛАУ не единственное";
        }
        float_type cb = BH[i];
        BH[i] = BH[imax];
        BH[imax] = cb;
        for (int j = i + 1; j < AH.lin(); ++j) {
            float_type k = AH[j][i] / AH[i][i];
            for (int l = i; l < AH.col(); ++l) {
                AH[j][l] = AH[j][l] - k * AH[i][l];
            }
            BH[j] = BH[j] - k * BH[i];
        }
    }
    if (abs(AH[AH.lin() - 1][AH.col() - 1]) < 1e-15) {
        throw "Решение СЛАУ не единственное";
    }
    out << "Приведение к ступенчатому виду" << endl;
    AH.output_f(out);
    out << endl << "Правая часть" << endl;
    tVectorOut(BH, out);
    vector<float_type> X(BH);
    for (int i = X.size() - 1; i >= 0; --i) {
        for (int j = i + 1; j < X.size(); ++j) {
            X[i] = X[i] - X[j] * AH[i][j];
        }
        X[i] = X[i] / AH[i][i];
    }
    return X;
}

template <typename float_type>
vector<float_type> QR_method(matrix<float_type> const& A, vector<float_type> const& B, ofstream& out)
//Алгоритм QR-разложения
{
    matrix<float_type> AH(A);
    vector<float_type> BH(B);
    matrix<float_type> Tm(A.col());
    Tm.init_unit();
    for (int i = 0; i < AH.lin() - 1; ++i) {
        for (int j = i + 1; j < AH.lin(); ++j) {
            if (abs(pow(AH[i][i], 2) + pow(AH[j][i], 2)) > 1e-15) {

                float_type c = AH[i][i] / sqrt(pow(AH[i][i], 2) + pow(AH[j][i], 2));
                float_type s = AH[j][i] / sqrt(pow(AH[i][i], 2) + pow(AH[j][i], 2));
                matrix<float_type> Th(A.col());
                for (int k = 0; k < AH.col(); ++k) {
                    float_type k1 = AH[i][k];
                    float_type k2 = AH[j][k];
                    AH[i][k] = c * k1 + s * k2;
                    AH[j][k] = -s * k1 + c * k2;
                }
                Th.init_unit();
                Th[i][i] = c;
                Th[i][j] = s;
                Th[j][i] = -s;
                Th[j][j] = c;
                Tm = Th * Tm;
            }
        }
    }


    matrix<float_type> Rm(Tm * A);
    matrix<float_type> Qm(Tm.transpose());
    out << "Матрица Q:" << endl;
    Qm.output_f(out);
    out << endl;
    if (abs(Rm[Rm.lin() - 1][Rm.col() - 1]) < 1e-15) {
        throw "Решение СЛАУ не единственное";
    }
    out << "Матрица R:" << endl;
    Rm.output_f(out);
    out << endl;
    BH = Tm * BH;
    out << "Вектор b*:" << endl;
    tVectorOut(BH, out);
    out << endl;
    vector<float_type> X(BH);
    for (int i = X.size() - 1; i >= 0; --i) {
        for (int j = i + 1; j < X.size(); ++j) {
            X[i] = X[i] - X[j] * Rm[i][j];
        }
        X[i] = X[i] / Rm[i][i];
    }
    return X;
}

template <typename float_type>
float_type vec_abs(vector<float_type> vec)
//Норма вектора ||*||_2
{
    float_type result = 0;
    for (int i = 0; i < vec.size(); ++i) {
        result += pow(vec[i], 2);
    }
    return sqrt(result);
}

template <typename float_type>
vector<float_type> operator-(vector<float_type> const& left, vector<float_type> const& right)
//Оператор разности векторов
{
    if (left.size() != right.size()) {
        throw "Ошибка размеров векторов при операции разности";
    }
    vector<float_type> result(left.size());
    for (int i = 0; i < left.size(); ++i) {
        result[i] = left[i] - right[i];
    }
    return result;
}

template <typename float_type>
vector<float_type> residual(matrix<float_type> const& A, vector<float_type> const& B, vector<float_type> const& X)
//Вычисления вектора невязки
{
    vector<float_type> B1 = A * X;
    vector<float_type> v_residual = B - B1;
    return v_residual;
}

template <typename float_type>
vector<float_type> rage(vector<float_type> const& vec, vector<float_type>& delta)
//Внесение возмущения в вектор
{
    vector<float_type> result(vec);
    delta = vec;
    for (int i = 0; i < vec.size(); ++i) {
        int pw = rand() % 2;
        delta[i] = 0.01 * pow(-1, pw);
        result[i] = vec[i] + 0.01 * delta[i];
    }
    return result;
}

template <typename float_type>
int search_component(matrix<float_type> const& A, vector<float_type> const& B, ofstream& crash)
//Поиск наиболее влияющей на решение компоненты вектора столбца.
{
    vector<float_type> X = gauss_method(A, B, crash);
    int imax = 0;
    float_type max = -1;
    for (int i = 0; i < B.size(); ++i) {
        vector<float_type> Bh(B);
        Bh[i] += 0.01;
        vector<float_type> Xh = gauss_method(A, Bh, crash);
        float_type num = vec_abs(Xh - X);
        if (num > max) {
            max = num;
            imax = i;
        }
    }
    return imax;
}

void lab1(char* const input1, char* const input2, char* const output_name)
{
    ofstream n_out;
    n_out.open(output_name);
    try {
        n_out << "Результаты работы программы:" << endl << endl;
        matrix<double> A(square, input1);
        n_out << "Исходная матрица:" << endl;
        A.output_f(n_out);
        n_out << endl;
        vector<double> B(read_vector<double>(input2));
        n_out << "Исходная правая часть СЛАУ:" << endl;
        tVectorOut(B, n_out);
        n_out << endl;
        vector<double> X1 = gauss_method(A, B, n_out);
        n_out << endl << "Решение по методу Гаусса" << endl;
        tVectorOut(X1, n_out);
        n_out << endl << endl << "Метод QR-разложения" << endl;
        vector<double> X2 = QR_method(A, B, n_out);
        n_out << "Результат QR-разложения" << endl;
        tVectorOut(X2, n_out);
        ofstream crash;
        crash.open(link_data"crash.txt");
        matrix<float> Af(square, input1);
        vector<float> Bf(read_vector<float>(input2));
        vector<float> X1f = gauss_method(Af, Bf, crash);
        vector<float> X2f = QR_method(Af, Bf, crash);
        n_out << endl << "Вектора невязки:" << endl;
        n_out << "1. Для типа float" << endl;
        n_out << "Для метода Гаусса: " << endl;
        tVectorOut(residual(Af, Bf, X1f), n_out);
        n_out << "Для QR-разложения: " << endl;
        tVectorOut(residual(Af, Bf, X2f), n_out);
        n_out << "2. Для типа double" << endl;
        n_out << "Для метода Гаусса: " << endl;
        tVectorOut(residual(A, B, X1), n_out);
        n_out << "Для QR-разложения: " << endl;
        tVectorOut(residual(A, B, X2), n_out);
        n_out << endl << "Модули вектора невязки:" << endl;
        n_out << "1. Для типа float" << endl;
        n_out << "Для метода Гаусса: " << vec_abs(residual(Af, Bf, X1f)) << ";  Для QR-разложения: " << vec_abs(residual(Af, Bf, X2f)) << endl;
        n_out << "2. Для типа double" << endl;
        n_out << "Для метода Гаусса: " << vec_abs(residual(A, B, X1)) << ";  Для QR-разложения: " << vec_abs(residual(A, B, X2)) << endl;
        n_out << endl << endl << "Число обусловленности по заданным нормам:" << endl;
        n_out << A.abs_1() * A.inverse().abs_1() << ";   " << A.abs_inf() * A.inverse().abs_inf() << endl;
        n_out << endl << "Внесение возмущения в B:" << endl;
        vector<double> deltaB(B);
        vector<double> Br = rage(B, deltaB);
        tVectorOut(Br, n_out);
        n_out << "Решение для возмущённой правой части:" << endl;
        vector<double> Xr = gauss_method(A, Br, crash);
        tVectorOut(Xr, n_out);
        vector<double> delta = Xr - X1;
        n_out << endl << "Реакция на возмущение:" << endl;
        tVectorOut(delta, n_out);
        n_out << endl << "Оценка снизу для числа обусловленности по заданным нормам:" << endl;
        double delta1_1 = abs_1(delta) / abs_1(X1);
        double delta2_1 = abs_1(deltaB) / abs_1(B);
        double delta1_inf = abs_inf(delta) / abs_inf(X1);
        double delta2_inf = abs_inf(deltaB) / abs_inf(B);
        n_out << delta1_1 / delta2_1 << ";  " << delta1_inf / delta2_inf << endl;
        n_out << endl << "Hомер компоненты вектора правой части, наиболее сильно влияющей на решение:" << endl;
        n_out << search_component(A, B, crash) + 1 << endl;
        n_out << endl << "Матрица A*A^(-1):" << endl;
        matrix<double> E(A * A.inverse());
        E.output_f(n_out);
    }
    catch (char const*& s)
        //Обработка ошибок выводит в файл название ошибки
    {
        n_out << endl << "Ошибка: " << s;
    }
}

int main() {
   
   lab1(link_data"matrix13_1.dat", link_data"vector13_1.dat", link_data"Kiseleva_13_1.txt");
    lab1(link_data"matrix13_2.dat", link_data"vector13_2.dat", link_data"Kiseleva_13_2.txt");
   lab1(link_data"matrix26_1.dat", link_data"vector26_1.dat", link_data"Chistyakova_26_1.txt");
    lab1(link_data"matrix26_2.dat", link_data"vector26_2.dat", link_data"Chistyakova_26_2.txt");
   
    return 0;
}
