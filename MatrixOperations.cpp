#include "MatrixOperations.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <array>
#include <chrono>
#include <random>

using namespace std;
using clock1 = chrono::high_resolution_clock;

// template<typename T>;

void test2()
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr1(10000, 11000);
    std::uniform_int_distribution<> distr2(10000, 15000);
    std::uniform_int_distribution<> distr3(10000, 13000);

    int iter = 0;

    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            int r1 = distr1(gen);
            int c1, r2;
            c1 = r2 = distr2(gen);
            int c2 = distr3(gen);
            MatrixDefinitionD A, B, C;
            A.fillRandom(r1, c1);
            A.printMatrixCSV(iter, "A");
            B.fillRandom(r2, c2);
            B.printMatrixCSV(iter, "B");
            C.fillzeros(r1, c2);

            chrono::time_point<chrono::system_clock> start = clock1::now();
            C = A * B;

            chrono::duration<float> execution_time = clock1::now() - start;
            string m = "C_calculated_in" + to_string((float)(execution_time.count() * pow(10, 9))) + "_nanoseconds";
            C.printMatrixCSV(iter, m);
            cout << "Time for Multiplication of matrices of size " << r1 << "x" << c1 << " and " << r2 << "x" << c2 << '\t' << execution_time.count() << endl;

            chrono::time_point<chrono::system_clock> start1 = clock1::now();
            C.transpose();
            C.transpose().printMatrixCSV(iter, "Transpose");
            chrono::duration<float> execution_time2 = clock1::now() - start1;
            cout << "Time for transpose of matrices of size " << r1 << "x" << c1 << " and " << r2 << "x" << c2 << '\t' << execution_time2.count() << endl;
            iter++;
        }
    }
}

void test1()
{
    const double arr1[] = {3, 5, 2,
                           8, 4, 8,
                           2, 4, 7};
    // MatrixDefinitionD test_d1(3, 3, vector<double>(arr1, arr1 + sizeof(arr1) / sizeof(arr1[0])));
    MatrixDefinitionD test_d1(6, 6);
    test_d1 = test_d1.fillMatrix(6, 6, 6);
    test_d1.size();

    MatrixDefinitionD d1 = test_d1.transpose();
    d1.size();
    // MatrixDefinitionD c1(6,6) ;

    d1 = test_d1 * d1;
    // std::cout << d1;
    MatrixDefinitionD c1;
    d1 += test_d1;
    // c1 = c1.StrassenMultiplication(test_d1,d1);
    d1.printMatrix();

    MatrixDefinitionD test_d11(6, 6);
    test_d11= test_d11.fillMatrix(6, 6, 6);
    MatrixDefinitionD test_d2(6, 6);
    test_d2 = test_d2.fillMatrix(6, 6, 16);
    MatrixDefinitionD r1, r2, r3;

    r1 = test_d11 + test_d2;
    r2 = test_d11 - test_d2;
    r3 = test_d11 / 65;

    r1.printMatrix();
    r2.printMatrix();
    r3.printMatrix();
}

int main()
{
    cout.precision(12);
    test1();
    // test2();
    return 0;
}