#pragma once

#include<string>
#include"../Types/Types.h"
#include"../cpp_libs/Print_Routines.h"
#include<hdf5.h>
#include<cmath>
#include<blaze/Math.h>
#include<random>

namespace print = Print_Routines;

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
    CouplingMatrix couplings; // Coupling matrix
    std::string spin_model;
    RealType beta; // Inverse temperature
    RealType rescale;

    // Numerical parameters
    uint Chebyshev_cutoff;
    std::string seed;
    uint num_TimePoints;
    uint num_TimeSteps_therm;
    uint num_Vectors_Per_Core;
    RealType Gauss_covariance; // Covariance of the distribution used to draw states
    RealType dt;
    RealType CET_rescale;
    RealType dt_therm;
    uint num_PrintDigits;

    // PUBLIC METHODS
    void read_SpinSystem();
    std::string create_essentials_string() const;    
};

}