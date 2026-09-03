#include "subdivision_oracle.hpp"

#include <vector>

void get_subdivision_matrices_oracle(Matrix &mat,
                                     Matrix &subMat1,
                                     Matrix &subMat2,
                                     Matrix &subMat3,
                                     Matrix &subMat4)
{
    const int N = 6;
    const double w = 3.0 / 8.0 / N; // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
    const int N1 = 5;
    const double w1 = 3.0 / 8.0 / N1; // w1=1/N1*(5/8-(3/8+1/4*cos(2*pi/N1))^2);
    const double a = 3.0 / 8.0;
    const double b = 1.0 / 8.0;
    std::vector<std::vector<double>> Mtmp{{a, b, a, b, 0, 0, 0, 0, 0, 0, 0},
                                {b, a, a, 0, 0, b, 0, 0, 0, 0, 0},
                                {w1, w1, 1.0 - N1 * w1, w1, 0, w1, w1, 0, 0, 0, 0},
                                {b, 0, a, a, 0, 0, b, 0, 0, 0, 0},
                                {0, a, b, 0, b, a, 0, 0, 0, 0, 0},
                                {0, b, a, 0, 0, a, b, 0, 0, 0, 0},
                                {0, 0, a, b, 0, b, a, 0, 0, 0, 0},
                                {0, 0, b, a, 0, 0, a, b, 0, 0, 0},
                                {0, b, 0, 0, a, a, 0, 0, b, 0, 0},
                                {0, w, w, 0, w, 1.0 - N * w, w, 0, w, w, 0},
                                {0, 0, b, 0, 0, a, a, 0, 0, b, 0},
                                {0, 0, w, w, 0, w, 1.0 - N * w, w, 0, w, w},
                                {0, 0, 0, b, 0, 0, a, a, 0, 0, b},
                                {0, 0, 0, 0, b, a, 0, 0, a, b, 0},
                                {0, 0, 0, 0, 0, a, b, 0, b, a, 0},
                                {0, 0, 0, 0, 0, b, a, 0, 0, a, b},
                                {0, 0, 0, 0, 0, 0, a, b, 0, b, a}};

    mat = Matrix(Mtmp);

    std::vector<std::vector<double>> SM1(12, std::vector<double>(17, 0.0));
    std::vector<int> element1{2, 3, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16};
    for (int i = 0; i < 12; i++)
    {
        SM1[i][element1[i]] = 1.0;
    }
    subMat1 = Matrix(SM1);

    std::vector<std::vector<double>> SM2(12, std::vector<double>(17, 0.0));
    std::vector<int> element2{4, 1, 9, 5, 2, 14, 10, 6, 3, 15, 11, 7};
    for (int i = 0; i < 12; i++)
    {
        SM2[i][element2[i]] = 1.0;
    }
    subMat2 = Matrix(SM2);

    std::vector<std::vector<double>> SM3(12, std::vector<double>(17, 0.0));
    std::vector<int> element3{1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15};
    for (int i = 0; i < 12; i++)
    {
        SM3[i][element3[i]] = 1.0;
    }
    subMat3 = Matrix(SM3);

    std::vector<std::vector<double>> SM4(11, std::vector<double>(17, 0.0));
    std::vector<int> element4{0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11};
    for (int i = 0; i < 11; i++)
    {
        SM4[i][element4[i]] = 1.0;
    }
    subMat4 = Matrix(SM4);
}
