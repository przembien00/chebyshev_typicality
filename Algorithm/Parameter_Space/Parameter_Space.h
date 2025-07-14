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
    ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank );

    // PUBLIC MEMBERS
    // Data storage
    std::string couplings_filename;
    std::string project_name;
    std::string filename_extension;
    int my_rank;
    int world_size;

    // Model parameters
    uint num_Spins;
    uint HilbertSpaceDimension;
    uint CoordinationNumber;
    CouplingMatrix couplings; // Coupling matrix
    std::string spin_model;
    RealType beta; // Inverse temperature
    RealType rescale;

    // Numerical parameters
    uint Chebyshev_cutoff;
    std::string seed;
    uint num_TimePoints;
    uint num_Vectors_Per_Core;
    RealType Gauss_covariance; // Covariance of the distribution used to draw states
    RealType dt;


    // PUBLIC METHODS
    void read_SpinSystem();
    
};

}