#pragma once

#include<string>
#include"../Types/Types.h"
#include<hdf5.h>
#include<cmath>
#include<blaze/Math.h>
#include<random>

namespace Parameter_Space
{

// ===============================================================
// ============== HEADER FOR PARAMETER SPACE CLASS ===============
// ===============================================================

class ParameterSpace
{
    
    public:

    // CONSTRUCTORS
    ParameterSpace() = default;
    ParameterSpace( const int argC, char* const argV[] );

    // PUBLIC MEMBERS
    // Data storage
    std::string couplings_filename;
    std::string project_name;
    
    // Model parameters
    uint num_Spins;
    uint HilbertSpaceDimension;
    uint CoordinationNumber;
    CouplingMatrix couplings; // Coupling matrix
    std::string spin_model;
    RealType beta; // Inverse temperature
    RealType rescale;

    // Numerical parameters
    uint num_TimePoints;
    uint num_Vectors;
    RealType Gauss_covariance; // Covariance of the distribution used to draw states
    RealType dt;


    // PUBLIC METHODS
    void read_SpinSystem();
    
};

}