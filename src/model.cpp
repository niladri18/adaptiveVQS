#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include "utils.h"
#include "model.h"

using namespace std;
/*
class Model{

    public:

    int nq {};
    std::string model_filename {};

    Model(int nqb = 0, std::string model_fname = "") : nq(nqb), model_filename(model_fname) {}

    //copy constructor
    Model(const Model& rhs) : nq(rhs.nq), model_filename(rhs.model_filename) {} 

    // Copy assignment operator
    Model& operator = (const Model& rhs){
        if (this != &rhs){
            nq = rhs.nq;
            model_filename = rhs.model_filename;
        }
        return *this;

    }
    //destructor
    ~Model(){}

    Matrix load_operator(string filename);
    Matrix Hamiltonian = this->load_operator( this->model_filename);

    
};
*/

Matrix Model::load_operator(string filename){

    SpinOp spin_op_list = build_ops(this->nq);
    //printSpinOp(spin_op_list);
    int ndim = pow(2,this->nq);// Hilbert space dimension

    //std::ifstream inputFile(this->model_filename);  // Open the file
    std::ifstream inputFile(filename);  // Open the file

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        //return 1;
    }

    std::string line;
    std::list<int> qb_list;
    std::list<string> pa_list;
    //std::list<string> out;
    string out, pstring;

    Matrix Ham( ndim, vector<Complex>(ndim) );

    while (std::getline(inputFile, line)) {
        // Process the line here
        std::cout << line << std::endl;
        out = line.substr(0, line.find("*"));
        pstring = line.substr(out.size()+1, line.size());
        double coeff = std::stof(out);
    
        

        //cout<<line.substr(out.size()+1, line.size())<<endl;
        //cout<<coeff<<"\t"<<pstring<<endl;
        
        qb_list = findall_digit(pstring);
        pa_list = findall_char(pstring);
        auto it1 = qb_list.begin();
        auto it2 = pa_list.begin();

        Matrix op;

        Matrix tmp( ndim, vector<Complex>(ndim) );
        for (int i=0; i< ndim; i++){
            tmp[i][i] = 1.0;
        }

        while (it1 != qb_list.end() && it2 != pa_list.end()) {


            //std::cout << *it1 << " " << *it2 << std::endl;

            op = spin_op_list[*it2][*it1];
            tmp = product(tmp, op);
            

            ++it1;
            ++it2;
        }

        //cout<<tmp.size()<<"\t"<<tmp[0].size()<<endl;

        Ham = add(Ham,coeff*tmp);



    }

    inputFile.close();  // Close the file
    //return 0;

    

    return Ham;

    
    
}


/*
int main(){
    
    int nq = 4;
    Model model = Model(nq, "N4_J1.0_h0.8.dat"); 
    //model.load_hamiltonian();
    
    printMatrix(model.Hamiltonian);
}
*/
