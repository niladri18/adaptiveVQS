#ifndef MODEL_H    
#define MODEL_H

#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include "utils.h"


using namespace std;

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
    ~Model(){}
    Matrix load_operator(string filename);
    Matrix Hamiltonian = this->load_operator( this->model_filename);
    Matrix h2 = product(Hamiltonian, Hamiltonian);

    
};


#endif
