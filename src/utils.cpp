#include "utils.h"
#include <iostream>
#include <regex>
#include <string>
#include <list>
#include <complex>
#include <vector>
#include <map>

using namespace std;



// Define the complex type for convenience
//using Complex = std::complex<double>;

//using Matrix = std::vector< vector<Complex>>;
// Define the Pauli matrices as 2x2 matrices

Matrix sigmaX { {0.0, 1.0}, {1.0, 0.0} }; 
Matrix sigmaY { {0.0, -Complex(0.0, 1.0)}, {Complex(0.0, 1.0), 0.0} }; 
Matrix sigmaZ { {1.0, 0.0}, {0.0, -1.0} }; 
Matrix sigmaI { {1.0, 0.0}, {0.0, 1.0} }; 


Matrix add( const Matrix& A, const Matrix& B){

    // Define matrix product
    int acol = A[0].size();
    int arow = A.size();

    Matrix C(arow,vector<Complex>(acol));

    for (int i=0; i< arow; i++){
        for(int j=0; j<acol; j++){

            C[i][j] = A[i][j] + B[i][j];

        }

    }
    return C;

}


//Matrix& operator = (const Matrix&);

Matrix operator * (const Complex& c, Matrix& A){
    /*
    Matrix multiplication by a scalar
    Arguments:
    c: Complex number scalar
    A: Matrix to be multiplied by c
    */

    int row = A.size();
    int col = A[0].size();
    //cout<<row<<""<<col<<endl;
    Matrix C ( row, vector<Complex> (col) );
    for(int i = 0; i<row; i++){
        for(int j = 0; j<col; j++){
            C[i][j] = c*A[i][j];
        }    
    }
    //printMatrix(C);
    return C;
} 


Matrix product( const Matrix& A, const Matrix& B){

    // Define matrix product
    int acol = A[0].size();
    int arow = A.size();
    int brow = B.size();
    int bcol = B[0].size();
    if (acol != brow) throw std::invalid_argument( "Column and row mismatch!" );

    Matrix C(arow,vector<Complex>(bcol));
    for(int i = 0; i< arow; i++){
        for(int j = 0; j<bcol; j++){
            C[i][j] = 0.0;
            for(int k = 0; k<acol; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}


// Function to compute the Kronecker product of two matrices
Matrix kroneckerProduct(const Matrix& matrix1, const Matrix& matrix2) {
    // Get the dimensions of the input matrices
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();

    // Create the resulting matrix with the appropriate dimensions
    Matrix result(rows1 * rows2, std::vector<Complex>(cols1 * cols2));

    // Compute the Kronecker product
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols1; ++j) {
            for (int k = 0; k < rows2; ++k) {
                for (int l = 0; l < cols2; ++l) {
                    result[i * rows2 + k][j * cols2 + l] = matrix1[i][j] * matrix2[k][l];
                }
            }
        }
    }

    return result;
}

Matrix exp(Matrix& A){
    /*
    Calculates the exponential of a square matrix
    argument: 
    A: input square matrix
    returns e^A
    */

    if (A.size() != A[0].size()) throw std::invalid_argument( "Matrix must be square to be exponentiated!" );
    //step 1: convert A to an eigen matrix A
    int n = A.size();

    Ematrix A1(n,n);
    //Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> A(n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            A1(i,j) = A[i][j];
        }
    }

    //Eigen::MatrixXd expA = A1.exp();
    //Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> expA = A1.exp(); 
    Ematrix expA = A1.exp();


    Matrix C(n,vector<Complex>(n));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            C[i][j] = expA(i,j);
        }
    }

    return C;

}


void printMatrix(const Matrix& A){

    int nrow = A.size();
    int ncol = A[0].size();
    for (int i=0; i< nrow; i++){
        for(int j=0; j<ncol;j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}


Psi operator *( const Matrix& A, const Psi& psi){

    /*
    Performs matrix-vector multiplication: C = Ax
    Parameters:
    psi: input statevector
    A: matrix operator

    return A*psi
    */
    int acol = A[0].size();
    int vrow = psi.size();

    Psi psi_(vrow,0); 

    if (acol != vrow) throw std::invalid_argument( "Column and row mismatch!" );

    for (int i=0; i<vrow; i++){
        for (int j=0; j<acol; j++){
            psi_[i] += A[i][j]*psi[j]; 
        }
    }

    return psi_;

} 

void normalize(Psi& psi){
    /*
    Normalize psi
    */

    Complex norm = 0;
    for(Complex x:psi){
        norm += std::conj(x)*x;
    }
    
    for(Complex &x:psi){
        x /= sqrt(norm);
    }

}

Complex vdot(Psi& psi1, Psi& psi2){
    /*
    dot product of two vectors
    arguments:
    |psi1>: first vector
    |psi2>: second vector
    returns <psi1 | psi2>
    */

    if (psi1.size() != psi2.size()) throw std::invalid_argument( "dimension mismatch for two vectors in function <vdot>!" );
    int n = psi1.size();

    Complex c = 0;
    for(int i=0; i<n ;i++){

        c += std::conj(psi1[i])*psi2[i];

    }
    return c;
}

void printPsi(Psi psi){
    cout<<"|"<<" ";
    for(Complex i: psi){
        cout<<i<<" ";
    }
    cout<<">"<<endl;
};

SpinOp build_ops(const int nq){

    /*
    Build the Pauli tensor product operators

    parameters:
    nq: number of qubits
    */

    int ndim = pow(2,nq);


    Matrix eye( ndim, vector<Complex>(ndim) );
    for (int i=0; i< ndim; i++){
        eye[i][i] = 1.0;
    }

    std::vector<Matrix> Oplist;
    
    for(int i = 0; i<nq; i++){
        Oplist.push_back(sigmaI);
    }
    
    // Initialize SX matrix

    std::vector<Matrix> sx_list, sy_list, sz_list;
    for (int i = 0; i<nq; i++){

        // Create SX
        Oplist[i] = sigmaX;
        Matrix init = Oplist[0];
        for (int j = 1; j < nq; j++){
            Matrix tmp;
            tmp = kroneckerProduct( init, Oplist[j] );
            init = tmp;
        }
        sx_list.push_back(init);


        // Create SY
        Oplist[i] = sigmaY;
        init = Oplist[0];
        for (int j = 1; j < nq; j++){
            Matrix tmp;
            tmp = kroneckerProduct( init, Oplist[j] );
            init = tmp;
        }
        sy_list.push_back(init);


        // Create SZ
        Oplist[i] = sigmaZ;
        init = Oplist[0];
        for (int j = 1; j < nq; j++){
            Matrix tmp;
            tmp = kroneckerProduct( init, Oplist[j] );
            init = tmp;
        }
        sz_list.push_back(init);

    }

    SpinOp spin_op_list { {"X",sx_list}, {"Y",sz_list}, {"Z",sz_list} };

    return spin_op_list;

}

void printSpinOp(const SpinOp& spin_op_list){

    for (const auto& item : spin_op_list){
        cout<< item.first <<endl;
        for (Matrix a:item.second){
            cout<<"=="<<endl;
            printMatrix(a);
        }
    } 

}



list<int> findall_digit(const string& text) {
    //std::string text = "X0Y1Z2Y4";
    std::regex pattern("\\d+"); // Regular expression pattern for matching words

    std::sregex_iterator it(text.begin(), text.end(), pattern);
    std::sregex_iterator end;
    std::list<int> out;

    while (it != end) {
        std::smatch match = *it;
        //std::cout << match.str() << std::endl;
        out.push_back(stoi(match.str()));
        ++it;
    }

    return out;
}


list<string> findall_char(const string& text) {
    std::regex pattern("[A-Z]"); // Regular expression pattern for matching words

    std::sregex_iterator it(text.begin(), text.end(), pattern);
    std::sregex_iterator end;
    //std::list<int> out;
    std::list<string> out;

    while (it != end) {
        std::smatch match = *it;
        //std::cout << match.str() << std::endl;
        out.push_back(match.str());
        ++it;
    }

    return out;
}


/*
int main(){
    //list<int> out;
    //list<string> out;
    //string str = "X0Y1Z2Y4";
    //out = findall_digit(str);
    //out = findall_char(str);
    //for(int i:out){
    //    cout<<i<<endl;
    //}
    //cout<<out<<endl;
    for (int i = 0; i < sigmaX.size(); i++){
        for (int j = 0; j < sigmaX[i].size(); j++){
        cout<<sigmaX[i][j];
        }
        cout<<endl;
    }

    
}
*/
