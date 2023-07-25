#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include <cstdio>
#include <tuple>
#include "utils.h"
#include "model.h"
#include "avqds.h"
#include "ansatz.h"



void AdaptSim::run(){

    // dimenstion of the Hilbert space
    //int ndim = pow(2,this->model.nq);
    int ndim = pow(2,this->ansatz.nq);

    // initialize the wavefunction
    Psi psi(ndim,1);
    normalize(psi);

    double dt = this->dt;

    const Complex im (0,1);
    Matrix Ham = this->model.Hamiltonian;
    Matrix Ham2 = this->model.h2;
    Matrix itHam = (-im*dt)*Ham; // -iHdt
    Matrix EvolutionOp = exp(itHam); // e^{-iHdt}
    printPsi(psi);
    normalize(psi);

    double t = 0.;
    Psi psi1;

    long energy;
    long var;

    while (t<this->T){

        psi1 = Ham*psi;
        //energy = vdot(psi1 , psi);
        tie(energy, var) = ansatz.adaptStep(Ham, Ham2, dt, psi);
        printf("t: %f energy: %ld\t var: %ld\n",t, energy,var );
        psi = EvolutionOp*psi;
        printPsi(psi);


        t += dt;

    }
    
    
}
