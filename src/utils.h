#ifndef UTILS_H    
#define UTILS_H

#include <stdio.h>
#include <string>
#include <list>
#include <complex>
#include <vector>
#include <map>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;


using Complex = std::complex<double>;
using Matrix = std::vector<vector<Complex>>;
using Psi = std::vector<Complex>;
using SpinOp = std::map<string, vector<Matrix>>;
using Ematrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;


Psi operator *( const Matrix& A, const Psi& psi);
void normalize(Psi& psi);
Complex vdot(Psi& psi1, Psi& psi2);
void printPsi(Psi psi);

Matrix product( const Matrix& A, const Matrix& B);
Matrix kroneckerProduct(const Matrix& matrix1, const Matrix& matrix2);
void printMatrix(const Matrix& A);
Matrix add( const Matrix& A, const Matrix& B);
Matrix operator * (const Complex& c, Matrix& A);
Matrix exp(Matrix& A);

//Matrix& operator = (const Matrix&);
SpinOp build_ops(const int nq);
void printSpinOp(const SpinOp& spin_op_list);

list<int> findall_digit(const string& text);
list<string> findall_char(const string& text);



#endif
