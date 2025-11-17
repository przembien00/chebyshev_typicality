#include"Parameter_Space.h"

#include<iostream>
#include<string>
#include<hdf5.h>
#include<blaze/Math.h>
#include<boost/program_options.hpp>
#include"../Types/Types.h"


namespace Parameter_Space
{

namespace bpo = boost::program_options;

// ===============================================================
// ==================== PARAMETER SPACE CLASS ====================
// ===============================================================
// constructor : build parameter space
ParameterSpace::ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank ):
    my_rank( my_rank ),
    world_size( world_size )
{
    // 1a.) Define Options
    bpo::options_description description("All Options:");
    bpo::options_description description_help("Help options:");
    bpo::options_description description_physics("General options concerning physics:");
    bpo::options_description description_numerics("General options concerning numerics:");
    bpo::options_description description_storing("General options concerning storing and naming:");

    description_help.add_options()
    (
    "help", "show all options"
    )(
    "helpNum", "show allowed options concerning numerics"
    );

    // ========== model and physical parameters ==========
    description_physics.add_options()
    (
    "srcfile", bpo::value<std::string>()->default_value("-"),
    "specify the file from which the couplings should be read; leave out file-ending hdf5"
    )(
    "spinmodel", bpo::value<std::string>()->default_value("ISO"),
    "set the spin model || options are : \
    ISO = Isotropic Heisenberg Model, \
    XXZ = Anisotropic Heisenberg Model"
    )(
    "beta", bpo::value<RealType>()->default_value(RealType{0.1}),
    "set the value of the inverse temperature"
    )(
    "h_z", bpo::value<RealType>()->default_value(RealType{0.}),
    "set the value of the magnetic field in the z direction"
    )(
    "lambda", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the anisotropy parameter lambda for the XXZ model"
    )(
    "symm_type", bpo::value<char>()->default_value('D'),
    "set the correlation symmetry type || options are: \
    A = (gab=0, gxx=gyy=gzz), \
    B = (gab=0, gxx=gyy), \
    C = (gaz=0, gxx=gyy), \
    D = (no constraints)"
    )(
    "evol_type", bpo::value<std::string>()->default_value("imaginary"),
    "set the type of time evolution || options are : \
    imaginary = imaginary time evolution, \
    real = real time evolution"
    )(
    "spin_site", bpo::value<uint>()->default_value(uint{0}),
    "set the spin with which the correlations are computed \
    (they look like < S_0^a(t) S_i^b(0) >, where i is spin_site and 0 the first spin in the coupling file)"
    )(
    "numSpins", bpo::value<uint>()->default_value(uint{16}),
    "set the number of spins in the system"
    )(
    "Tmax", bpo::value<RealType>()->default_value(RealType{5.0}),
    "set the maximum time for real time evolution"
    );

    // ========== general numerical parameters ==========
    description_numerics.add_options()
    (
    "seed", bpo::value<std::string>()->default_value("random"),
    "set the seed for drawing the states"
    )(
    "numCouplingConfigs", bpo::value<uint>()->default_value(uint{1}),
    "set the number of different coupling configurations to be drawn"
    )(
    "numTimePoints", bpo::value<uint>()->default_value(uint{100}),
    "set the number of time points for the equidistant time discretization"
    )(
    "fulldiag", "perform a loop over the whole Hilbert space for exact diagonalization"
    )(
    "GaussCovariance", bpo::value<RealType>()->default_value(RealType{1.0}), 
    "set the covariance of the distribution used to draw random states"
    )(
    "numVectorsPerCore", bpo::value<uint>()->default_value(uint{1}),
    "set the number of vectors drawn and averaged over"
    )(
    "CET_therm_error", bpo::value<RealType>()->default_value(RealType{1e-10}),
    "set the error threshold for thermalization, used to determine the Chebyshev expansion depth"
    )(
    "CET_evol_error", bpo::value<RealType>()->default_value(RealType{1e-6}),
    "set the error threshold for time evolution, used to determine the Chebyshev expansion depth"
    )(
    "determine_bandwidth", bpo::value<bool>()->default_value(true),
    "if true, the bandwidth of the Hamiltonian is determined via a Lanczos procedure \
    if false, the values E_min and E_max are used"
    )(
    "E_max", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the maximum eigenvalue of the Hamiltonian, needed for the Chebyshev expansion"    
    )(
    "E_min", bpo::value<RealType>()->default_value(RealType{-1.0}),
    "set the minimum eigenvalue of the Hamiltonian, needed for the Chebyshev expansion"
    );

    // ========== storing and naming ==========
    description_storing.add_options()
    (
    "project", bpo::value<std::string>()->default_value(""),
    "sort the data into a project folder"
    )(
    "fileext", bpo::value<std::string>()->default_value(""),
    "Define an extension to the filename; it will be appended according to : filename__fileext"
    )(
    "numPrintDigits", bpo::value<size_t>()->default_value(size_t{4}),
    "set the value precision for printing to the terminal"
    );

    description.add( description_help );
    description.add( description_physics );
    description.add( description_numerics );
    description.add( description_storing );
    bpo::variables_map vm;
    bpo::store( bpo::parse_command_line(argC, argV, description), vm );
    bpo::notify( vm );

    // 2.) Output of help descriptions and termination of the program
    if( vm.count("help") || vm.count("helpNum") )
    {
        if( vm.count("help") ){
            std::cout << description << "\n";
        }
        exit(0);
    }

    // 3.) Store bp options in parameter space 

    // ========== physical parameters ==========
    num_Spins = vm["numSpins"].as<uint>();
    HilbertSpaceDimension = 1 << num_Spins; // 2^num_Spins
    beta = vm["beta"].as<RealType>();
    spin_model = vm["spinmodel"].as<std::string>();
    h_z = vm["h_z"].as<RealType>();
    lambda = vm["lambda"].as<RealType>();
    symmetry_type = vm["symm_type"].as<char>();
    evol_type = vm["evol_type"].as<std::string>();
    spin_site = vm["spin_site"].as<uint>();

    // ========== general numerical parameters ==========
    num_Coupling_Configs = vm["numCouplingConfigs"].as<uint>();
    seed = vm["seed"].as<std::string>();
    CET_therm_error = vm["CET_therm_error"].as<RealType>();
    CET_evol_error = vm["CET_evol_error"].as<RealType>();
    num_TimePoints = vm["numTimePoints"].as<uint>();
    if( evol_type == "real" )
    {
        Tmax = vm["Tmax"].as<RealType>();
    }
    else if( evol_type == "imaginary" )
    {
        Tmax = beta * RealType{0.5};
    }
    dt = Tmax / static_cast<RealType>(num_TimePoints - 1);
    num_Vectors_Per_Core = vm["numVectorsPerCore"].as<uint>();
    Gauss_covariance = vm["GaussCovariance"].as<RealType>();
    determine_bandwidth = vm["determine_bandwidth"].as<bool>();
    if( vm.count("fulldiag") )
    {
        full_diagonalization = true;
        num_Vectors_Per_Core = HilbertSpaceDimension / static_cast<uint>(world_size) + 1;
    }
    E_max = vm["E_max"].as<RealType>();
    E_min = vm["E_min"].as<RealType>();

    // ========== saving and naming ==========
    project_name = vm["project"].as<std::string>();
    filename_extension = vm["fileext"].as<std::string>();
    num_PrintDigits = vm["numPrintDigits"].as<size_t>();
}

void ParameterSpace::draw_couplings( uint seed, uint config )
{
    std::mt19937 gen{ static_cast<uint>( seed + 1000 * config ) };
    std::normal_distribution<RealType> d{0., 1.0};
    RealType JQ = RealType{0.};
    couplings.resize( num_Spins );
    for( uint i = 0; i < num_Spins; ++i )
    {
        for( uint j = i+1; j < num_Spins; ++j )
        {
            couplings(i,j) = d(gen);
        }
        JQ += std::pow(couplings(0,i),2);
    }
    JQ = std::sqrt(JQ);
    couplings /= JQ; // Normalize couplings such that JQ=1
}

// printing method : return the essential parameters string
std::string ParameterSpace::create_essentials_string() const
{   
    size_t pre_colon_space = 35;
    std::stringstream ss{};

    ss
    << print::quantity_to_output_line( pre_colon_space, "spin_model"     , spin_model )
    << print::quantity_to_output_line( pre_colon_space, "num_Spins"     , std::to_string(num_Spins) )
    << print::quantity_to_output_line( pre_colon_space, "beta"          , print::remove_zeros(print::round_value_to_string(beta,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "h_z"          , print::remove_zeros(print::round_value_to_string(h_z,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "num_TimePoints" , std::to_string(num_TimePoints) ) 
    << print::quantity_to_output_line( pre_colon_space, "dt"       , print::remove_zeros(print::round_value_to_string(dt,num_PrintDigits)) ) 
    << print::quantity_to_output_line( pre_colon_space, "num_Vectors"   , std::to_string(num_Vectors_Per_Core*world_size) );
    if( spin_model == "XXZ" )
    {
        ss << print::quantity_to_output_line( pre_colon_space, "lambda"          , print::remove_zeros(print::round_value_to_string(lambda,num_PrintDigits)) );
    }
    
    return ss.str();
}

}