#ifndef MATRIX_OPERATIONS_HPP
#define MATRIX_OPERATIONS_HPP

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <set>
#include <cmath>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;

template <
    typename T,
    typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class MatrixDefinition
{
private:
    /* data */
    std::vector<T> data;
    size_t num_rows;
    size_t num_cols;
    void checkIndexValidity(size_t rows, size_t cols) const
    {
        if (rows < 0 or rows >= num_rows)
            throw invalid_argument("Invalid number of rows");
        if (cols < 0 or cols >= num_cols)
            throw invalid_argument("Invalid number of columns");
    };
    // static pair<MatrixDefinition,MatrixDefinition> eigenSorting(MatrixDefinition eigenval,MatrixDefinition eigenvec);
public:
    enum Axis
    {
        All,
        Rows,
        Cols
    };
    size_t numCols() const { return num_cols; }
    size_t numRows() const { return num_rows; }

    // Initialize empty matrix
    MatrixDefinition()
    {
        num_rows = num_cols = 0;
    }
    ~MatrixDefinition(){

    };

    // Initialize a square matrix of Dimension d
    MatrixDefinition(size_t d)
    {
        MatrixDefinition(d, d);
    }

    MatrixDefinition(size_t r, size_t c) : num_rows(r), num_cols(c), data(r * c)
    {
    }

    MatrixDefinition(size_t r, size_t c, const vector<T> &MatrixData) : num_rows(r), num_cols(c)
    {
        if (r * c != MatrixData.size())
            throw invalid_argument("Please initialize the matrix with size compatible to vector data");
        data = MatrixData;
    }

    template <std::size_t N>
    MatrixDefinition(size_t r, size_t c, T (&data)[N])
    {
        if (N != r * c)
            throw invalid_argument("Matrix dimension incompatible with its initiallizing vector.");
        vector<T> v(data, data + N);
        MatrixDefinition(r, c, v);
    }

    // Scalar Arithmetic Operations (POST Matrix)

    friend MatrixDefinition operator+(const MatrixDefinition &matrix, double scalar_value)
    {
        MatrixDefinition resultingMatrix(matrix.num_rows, matrix.num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < matrix.num_rows; i++)
        {
            for (size_t j = 0; j < matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = scalar_value + matrix(i, j);
            }
        }

        return resultingMatrix;
    }

    friend MatrixDefinition operator-(const MatrixDefinition &matrix, double scalar_value)
    {
        return matrix + (-scalar_value);
    }

    friend MatrixDefinition operator*(const MatrixDefinition &matrix, double scalar_value)
    {
        MatrixDefinition resultingMatrix(matrix.num_rows, matrix.num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < matrix.num_rows; i++)
        {
            for (size_t j = 0; j < matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = scalar_value * matrix(i, j);
            }
        }

        return resultingMatrix;
    }

    friend MatrixDefinition operator/(const MatrixDefinition &matrix, double scalar_value)
    {
        MatrixDefinition resultingMatrix(matrix.num_rows, matrix.num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < matrix.num_rows; i++)
        {
            for (size_t j = 0; j < matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = matrix(i, j) / scalar_value;
            }
        }

        return resultingMatrix;
    }

    // Scalar Arithmetic Operations (PRE Matrix)

    friend MatrixDefinition operator+(double scalar_value, const MatrixDefinition matrix)
    {
        return matrix + scalar_value;
    }

    friend MatrixDefinition operator-(double scalar_value, const MatrixDefinition matrix)
    {
        return matrix - scalar_value;
    }

    friend MatrixDefinition operator*(double scalar_value, const MatrixDefinition matrix)
    {
        return matrix * scalar_value;
    }

    friend MatrixDefinition operator/(double scalar_value, const MatrixDefinition &matrix)
    {
        MatrixDefinition resultingMatrix(matrix.num_rows, matrix.num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < matrix.num_rows; i++)
        {
            for (size_t j = 0; j < matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = scalar_value / matrix(i, j);
            }
        }

        return resultingMatrix;
    }

    // Scalar Assignment Operations

    MatrixDefinition operator+=(double scalar_value)
    {
#pragma omp parallel for
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] += scalar_value;
        }
        return *this;
    }

    MatrixDefinition operator-=(double scalar_value)
    {
#pragma omp parallel for
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] -= scalar_value;
        }
        return *this;
    }

    MatrixDefinition operator*=(double scalar_value)
    {
#pragma omp parallel for
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] *= scalar_value;
        }
        return *this;
    }

    MatrixDefinition operator/=(double scalar_value)
    {
#pragma omp parallel for
        for (size_t i = 0; i < data.size(); i++)
        {
            data[i] /= scalar_value;
        }
        return *this;
    }

    // Matrix Operations

    // Matrix Addition
    MatrixDefinition operator+(const MatrixDefinition &new_matrix) const
    {
        if (num_rows != new_matrix.num_rows || num_cols != new_matrix.num_cols)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }

        MatrixDefinition resultingMatrix(num_rows, num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t j = 0; j < new_matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = operator()(i, j) + new_matrix(i, j);
            }
        }
        return resultingMatrix;
    }

    // Matrix Subtraction
    MatrixDefinition operator-(const MatrixDefinition &new_matrix) const
    {
        if (num_rows != new_matrix.num_rows || num_cols != new_matrix.num_cols)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }

        MatrixDefinition resultingMatrix(num_rows, num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t j = 0; j < new_matrix.num_cols; j++)
            {
                resultingMatrix(i, j) = operator()(i, j) - new_matrix(i, j);
            }
        }
        return resultingMatrix;
    }

    // Matrix Multiplication
    MatrixDefinition operator*(const MatrixDefinition &new_matrix) const
    {
        if (num_cols != new_matrix.num_rows)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }

        MatrixDefinition resultingMatrix = fillzeros(num_rows, new_matrix.num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t k = 0; k < new_matrix.num_cols; k++)
            {
                double temp = operator()(i, k);
                for (size_t j = 0; j < resultingMatrix.num_cols; j++)
                    resultingMatrix(i, j) = temp * new_matrix(k, j);
            }
        }
        return resultingMatrix;
    }

    // Matrix Assignment Operations
    MatrixDefinition &operator+=(const MatrixDefinition &new_matrix)
    {
        if (num_rows != new_matrix.num_rows || num_cols != new_matrix.num_cols)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t j = 0; j < new_matrix.num_cols; j++)
            {
                operator()(i, j) += new_matrix(i, j);
            }
        }

        return *this;
    }

    MatrixDefinition &operator-=(const MatrixDefinition &new_matrix)
    {
        if (num_rows != new_matrix.num_rows || num_cols != new_matrix.num_cols)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t j = 0; j < new_matrix.num_cols; j++)
            {
                operator()(i, j) -= new_matrix(i, j);
            }
        }

        return *this;
    }

    MatrixDefinition &operator*=(const MatrixDefinition &new_matrix)
    {
        if (num_cols != new_matrix.num_rows)
        {
            throw invalid_argument("Cannot operate on these Matrices because their sizes are LeftMatrix:" + to_string(num_rows) + "rows by " + to_string(num_cols) + "columns and RightMatrix:" + to_string(new_matrix.num_rows) + "rows by " + to_string(new_matrix.num_cols) + "columns.");
        }

        MatrixDefinition resultingMatrix(num_rows, num_cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < new_matrix.num_rows; i++)
        {
            for (size_t k = 0; k < new_matrix.num_cols; k++)
            {
                resultingMatrix(i, k) = 0;
                for (size_t j = 0; j < resultingMatrix.num_cols; j++)
                    resultingMatrix(i, k) = operator()(i, j) * new_matrix(j, k);
            }
        }

        num_rows = resultingMatrix.num_rows;
        num_cols = resultingMatrix.num_cols;
        data = resultingMatrix.data;
        return *this;
    }

    T &operator()(size_t i, size_t j)
    {
        checkIndexValidity(i, j);
        return data[i * num_cols + j];
    }

    T operator()(size_t i, size_t j) const
    {
        checkIndexValidity(i, j);
        return data[i * num_cols + j];
    }

    //  Helper Functions

    static MatrixDefinition fillMatrix(size_t r, size_t c, double value)
    {
        MatrixDefinition result(r, c, vector<T>(r * c, value));
        return result;
    }

    static MatrixDefinition fillRandom(size_t r, size_t c)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr1(1, 100);
        vector<T> V1(r * c, 0);
        for (int i=0;i<V1.size();i++) 
            V1.at(i) = distr1(gen);
        MatrixDefinition result(r, c, V1);
        return result;
    }

    static MatrixDefinition fillzeros(size_t r, size_t c)
    {
        return fillMatrix(r, c, 0);
    }

    bool isSquare() const
    {
        return num_cols == num_rows;
    }

    MatrixDefinition subMatrix(size_t r, size_t c) const
    {
        MatrixDefinition resultingMatrix(num_rows - 1, num_cols - 1);

        size_t si = 0;
#pragma omp parallel for
        for (size_t i = 0; i < num_rows; i++)
        {
            size_t subj = 0;
            if (i == r)
                continue;
            for (size_t j = 0; j < num_cols; j++)
            {
                if (j == c)
                    continue;
                resultingMatrix(si, subj) = operator()(i, j);
                subj++;
            }
            si++;
        }
    }

    MatrixDefinition transpose() const
    {
        MatrixDefinition resultingMatrix(num_cols, num_rows);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < resultingMatrix.num_rows; i++)
        {
            for (size_t j = 0; j < resultingMatrix.num_cols; j++)
            {
                resultingMatrix(j, i) = operator()(i, j);
            }
        }

        return resultingMatrix;
    }

    double determinant() const
    {
        if (!isSquare())
        {
            throw runtime_error("Not a Square matrix cannot calculate determinant");
        }

        size_t num = num_rows;
        double det = 0;
        if (num == 2)
        {
            return (operator()(0, 0) * operator()(1, 1) - operator()(1, 0) * operator()(0, 1));
        }
        else
        {
#pragma omp parallel for reduction(+ \
                                   : det)
            for (size_t i = 0; i < num; i++)
            {
                det += pow(-1, i) * operator()(0, i) * subMatrix(0, i).determinant();
            }
        }
        return det;
    }

    void size()
    {
        cout << "rows" << num_rows << ", columns" << num_cols << endl;
        // return 0;
    }

    void printMatrix()
    {
        

        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << to_string(operator()(i, j)) << ",";
            }
            cout << endl;
        }
        
    }

    void printMatrixCSV(int i, string c)
    {
        
        ostringstream filename;
        filename << "Matrix_"<<c<<"for_iteration_"<<to_string(i)<<".csv";
        std::ofstream out;
        out.open(filename.str());

        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                out << to_string(operator()(i, j)) << ",";
            }
            out << endl;
        }
    }

    
};

typedef MatrixDefinition<double> MatrixDefinitionD;
typedef MatrixDefinition<int> MatrixDefinitionI;

#endif