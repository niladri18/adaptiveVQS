#ifndef ANSATZ_H    
#define ANSATZ_H


#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <map>
#include <tuple>
#include "utils.h"

using namespace std;
using OpPool = std::map<string, Matrix>;
using AnsatzOp = std::map<std::pair<long,std::string>, Matrix>; // store the ansatz 

class Ansatz{
    public:
    int nq {};
    std::string pool_filename {};
    int nparam = 0;

    Ansatz(int nqb = 0, std::string fname=""): nq(nqb),  pool_filename(fname){}
    Ansatz(const Ansatz& rhs):nq(rhs.nq), pool_filename(rhs.pool_filename){}
    Ansatz& operator = (const Ansatz& rhs){
        if (this != &rhs){
            nq = rhs.nq; 
            pool_filename = rhs.pool_filename;
        }
        return *this;
        }
    ~Ansatz(){}
    //OpPool load_operator_pool();
    OpPool load_operator_pool();
    OpPool operator_pool = this->load_operator_pool();
    std::tuple<long, long> adaptStep( Matrix& hamiltonian, Matrix& h2, double dt, Psi& psi);
    AnsatzOp ansatz;
    //OpPool load_operator_pool_parallel();
    //OpPool operator_pool2 = this->load_operator_pool_parallel();

};


#endif
