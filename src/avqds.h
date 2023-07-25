#ifndef AVQDS_H    
#define AVQDS_H

#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <list>
#include "utils.h"
#include "model.h"
#include "ansatz.h"


class AdaptSim{

    public:

    Model model {};
    Ansatz ansatz {};
    double T;
    double dt;
    double rcut;

    //constructor
    //AdaptSim( Ansatz& a={}, const Model& m = {}, const double T = 0, const double rcut = 1e-4): T(T), ansatz(a), model(m), T(T), rcut(rcut) {} 
    AdaptSim( const Model& m = {}, \
            const Ansatz& a={}, \
            const double T = 0, \
            const double dt = 0.1, \
            const double rcut = 1e-4): \
            model(m), ansatz(a), T(T), dt(dt), rcut(rcut) {} 

    //copy constructor
    AdaptSim(const AdaptSim& rhs) : \
        model(rhs.model), \
        ansatz(rhs.ansatz), \
        T(rhs.T), \
        dt(rhs.dt),\
        rcut(rhs.rcut) \
        {}

    //copy assignment operator
    AdaptSim& operator = (const AdaptSim& rhs){
        if(this != &rhs){
            model = rhs.model;
            ansatz = rhs.ansatz;
            T = rhs.T;
            dt = rhs.dt;
            rcut = rhs.rcut;
        }
        return *this;
    }
    ~AdaptSim(){}


    void run();

};



#endif
