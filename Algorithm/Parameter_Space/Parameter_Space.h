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
    RealType h_z;
    char symmetry_type;
    std::string evol_type;
    uint spin_site;

    // Numerical parameters
    RealType CET_therm_error; // Error threshold for thermalization
    RealType CET_evol_error; // Error threshold
    std::string seed;
    uint num_TimePoints;
    uint num_Vectors_Per_Core;
    RealType Gauss_covariance; // Covariance of the distribution used to draw states
    RealType dt;
    RealType CET_rescale;
    uint num_PrintDigits;
    RealType Tmax;
    bool determine_bandwidth;
    RealType E_max;
    RealType E_min;

    // PUBLIC METHODS
    void read_SpinSystem();
    std::string create_essentials_string() const;    
};

}