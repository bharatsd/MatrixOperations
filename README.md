#How to use MatrixOperations in your project

> Initializing options for matrix of different diensions

```cpp
#include <iostream>
#include "MatrixOperations.hpp"

using namespace std;

int main (void) {
  MatrixDefinitionD test_1; // initialize Matrix

  // code in the examples goes here...

  return 0;
}
```

> Replace `// code in the examples goes here...` with your own after going through the header, save this file, `Matrix_Code.cpp` and compile it like this: `g++ -pipe -O2 -std=c++11 Matrix_Code.cpp `, where `Matrix` is the path where you have the header files which you obtained after cloning this repo.

##Constructing Matrix

There are several ways to initialize a _Matrix_ object.

- Overloaded with number of rows and columns

  ```cpp
  MatrixDefinitionD test_2(6,16);
  ```

- Square matrix of dimension dxd

  ```cpp
  MatrixDefinitionD test_3(9);
  ```

- Uses data from vector to fill matrix

  ```cpp
  MatrixDefinitionD test_4(3,3,vector<double>(9);
  ```

The most complete example is `MatrixOperations.cpp`. Here is a list of all mathematical operations currently available and examples of usage.

- multiplication

  - `matrix * matrix` and `matrix * number` (scalar operation interchangeable)

  ```cpp
    MatrixDefinitionD A, B, C;
    A.fillRandom(10, 10);
    B.fillRandom(10, 10);

  C = A* B;

  cout << c * 3.1 << endl;

  ```

- transpose

- `matrix'` , helps compute the transpose of the Matrix

```cpp
    A1.fillRandom(10, 10);
    

  A1.transpose();

  ```


- Additional functions
-Matrix Addition, Matrix subtraction, scalar division,scalar addition, scalar subtraction
- `matrix + matrix`, `matrix - matrix` nad `matrix / scalar`

```cpp
    

MatrixDefinitionD test_d1(6, 6);
MatrixDefinitionD test_d2(6, 6);
MatrixDefinitionD r1,r2,r3,r4,r5,r6,r7,r8,r9;

r1 = test_d1+test_d2;
r2 = test_d1-test_d2;
r3 = test_d1/65;
r4 = test_d2-12;
r5 = test_d2 +242;
```
-Assignment Operations with both Scalars and Matrices
- `matrix +=`, `matrix -=` nad `matrix*=`

```cpp
test_d1+=r1; //Matrix test_d1 = test_d1 + r1 
test_d2+=9; // Add 9 to all the roe and column elements

test_d1-=r1; //Matrix test_d1 = test_d1 - r1 
test_d2-=9; // Subtract 9 to all the roe and column elements

test_d1*=r1; //Matrix test_d1 = test_d1 * r1 
test_d2*=9; // Multiplies 9 to all the roe and column elements
 
test_d2/=9; // Divides 9 to all the roe and column elements
```

- Indexing
```cpp

MatrixDefinitionD test_d1(6, 6);
cout<<test_d1(m,n)<<endl; //outputs the mth row nth coulumn element

```

-Helper Functions

```cpp

MatrixDefinitionD test_d1,test_d2,test_d3,test_d4,test_d5;
test_d1.fillMatrix(6,6,26); // creates a 6x6 matrix and fills all elements to be 26
test_d2.fillzeros(7,8);  // creates a 7x8 matrix and fills all elements to be 0
test_d3.fillRandom(8,4); //Random numbers as elements of 8x4 matrix
test_d1.determinant(); //determinant of a square matrix
test_d2.printMatrix(); // prints to console
test_d3.printMatrixCSV(1,string); //prints to csv file with a unique name made of the arguements

```

- Use `MatrixOperations.cpp` to run the test cases for multiplication and transpose of two matrices, and compile it like this: `g++ -pipe -O2 -std=c++11 MatrixOperations.cpp ` and complete by typing `./a.out` in the console
