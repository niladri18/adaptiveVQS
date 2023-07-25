#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <map>
#include <thread>
#include <future>
#include <tuple>
#include "utils.h"
#include "model.h"
#include "avqds.h"
#include "ansatz.h"

using namespace std;
using OpPool = std::map<string, Matrix>;

/*
std::promise<Matrix> buildOp(std::string& pstring, SpinOp spin_op_list, const int nq){
    int ndim = pow(2,nq);
    std::string line;
    std::string pauli;

    Matrix tmp( ndim, vector<Complex>(ndim) );
    for (int i=0; i< ndim; i++){
        tmp[i][i] = 1.0;
    }

    int n = pstring.size();
    for(int i = 0; i<n; i++){

        

        pauli = pstring[i];
        if (pauli == "I") continue;
        Matrix op = spin_op_list[ pauli ][i];
        tmp = product(op,tmp);
        
    }
        //op_pool[line] = tmp;
        return tmp;

}
*/

OpPool Ansatz::load_operator_pool(){

    /*
    Builds the operator pool by reading from a text file
    */
    
    SpinOp spin_op_list = build_ops(this->nq);
    //printSpinOp(spin_op_list);
    int ndim = pow(2,this->nq);// Hilbert space dimension


    std::ifstream inputFile(this->pool_filename);  // Open the file

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        //return 1;
    }

    std::string line;
    std::string pauli;

    OpPool op_pool;


    std::vector<std::string> pool_list;
    while (std::getline(inputFile, line)) {
        // Process the line here
        std::cout << line << std::endl;
        pool_list.push_back(line);
    }

    for (std::string line: pool_list){

        int n = line.size();

        Matrix tmp( ndim, vector<Complex>(ndim) );
        for (int i=0; i< ndim; i++){
            tmp[i][i] = 1.0;
        }

        for(int i = 0; i<n; i++){

            

            pauli = line[i];
            if (pauli == "I") continue;
            Matrix op = spin_op_list[ pauli ][i];
            tmp = product(op,tmp);
            
        }
        op_pool[line] = tmp;


    }

    return op_pool;
}


std::tuple<long, long> Ansatz::adaptStep( Matrix& hamiltonian, Matrix& h2, double dt, Psi& psi){
    /*
    Runs one Trotter step: checks McLachlan's distance and runs adaptive procedure
    Params:
    hamiltonian: model Hamiltonian
    h2: Hamiltonian^2
    dt: time stepsize
    psi: wavefunction at time t
    */

    Psi psi1 = hamiltonian*psi;
    Psi psi2 = h2*psi;
    Complex energy = vdot(psi, psi1);
    Complex h2exp = vdot(psi, psi2);
    Complex var = h2exp - (energy*energy);
    //cout<<"Var: "<<var<<std::real(var)<<endl;
    
    return make_tuple(std::real(energy), std::real(var));
}


/*
std::map<std::string, std::vector<Psi>>Ansatz::compute_dpsi(Psi& psi){




    for(const auto& item: this->ansatz){

        std::pair<long,std::string> key = item.first;
        Matrix op = item.second;

        

    }

}
*/
    


/*
OpPool Ansatz::load_operator_pool_parallel(){

    
    SpinOp spin_op_list = build_ops(this->nq);
    //printSpinOp(spin_op_list);
    int nq = this->nq;


    std::ifstream inputFile(this->pool_filename);  // Open the file

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        //return 1;
    }

    std::string line;
    std::string pauli;

    OpPool op_pool;
    std::promise<Matrix> tmp1;
    std::future<Matrix> tmp2 = tmp1.get_future();


    std::vector<std::thread> threads;
    std::vector<std::string> pool_list;
    while (std::getline(inputFile, line)) {
        // Process the line here
        std::cout << line << std::endl;
        pool_list.push_back(line);
    }


    for (std::string line: pool_list){

        auto tmp = std::async(buildOp, line, spin_op_list, nq);
        Matrix ops = tmp.get(); 

    }

    // Join the threads to wait for them to finish
    for (auto& thread : threads) {
        thread.join();
    }
    //for (auto& thread: threads){
    //    cout<<thread<<endl;
    //}

    return op_pool;
}
*/

