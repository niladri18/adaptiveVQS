#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include "utils.h"
#include "model.h"
#include "avqds.h"
#include "ansatz.h"


using namespace std;


int main(){

    // Define the number of qubits
    int nq = 4; 
    double T = 1;
    double dt = 0.1;
    double rcut = 1e-6;

    // Hamiltonian file name
    std::string model_filename = "N4_J1.0_h0.8.dat";

    // load Hamiltonian
    Model ising = Model(nq, model_filename);

    //Operator pool file name
    std::string operator_filename = "N4hampool.dat";
    Ansatz ansatz = Ansatz(nq, operator_filename);

    std::cout<<"Operator pool created"<<endl;

    //simulator
    AdaptSim simulator = AdaptSim(ising, ansatz, T, dt, rcut);


    // print Hamiltonain
    //printMatrix(ising.Hamiltonian);
    //printMatrix(simulator.model.Hamiltonian);

    simulator.run();
}
