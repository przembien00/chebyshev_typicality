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
    ISO = Isotropic Heisenberg Model, "
    )(
    "beta", bpo::value<RealType>()->default_value(RealType{0.1}),
    "set the value of the inverse temperature"
    )(
    "rescale", bpo::value<RealType>()->default_value(RealType{1.0}),
    "the factor rescale is multiplied to the couplings (for example for them to be in units of JQ)"
    );

    // ========== general numerical parameters ==========
    description_numerics.add_options()
    (
    "seed", bpo::value<std::string>()->default_value("random"),
    "set the seed for drawing the states"
    )(
    "numTimePoints", bpo::value<uint>()->default_value(uint{100}),
    "set the number of time points for the equidistant time discretization"
    )(
    "GaussCovariance", bpo::value<RealType>()->default_value(RealType{1.0}), 
    "set the covariance of the distribution used to draw random states"
    )(
    "numVectorsPerCore", bpo::value<uint>()->default_value(uint{1}),
    "set the number of vectors drawn and averaged over"
    )(
    "ChebyshevRescale", bpo::value<RealType>()->default_value(RealType{4}),
    "set the rescaling of the Hamiltonian for CET"
    )(
    "dt", bpo::value<RealType>()->default_value(RealType{0.05}),
    "set the timestep, RK4 error scales like dt^5"
    )(
    "CET_therm_error", bpo::value<RealType>()->default_value(RealType{1e-10}),
    "set the error threshold for thermalization, used to determine the Chebyshev expansion depth"
    )(
    "CET_evol_error", bpo::value<RealType>()->default_value(RealType{1e-7}),
    "set the error threshold for time evolution, used to determine the Chebyshev expansion depth"
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
    couplings_filename      = vm["srcfile"].as<std::string>();
    rescale = vm["rescale"].as<RealType>();
    read_SpinSystem(); // couplings and num_Spins are set
    HilbertSpaceDimension = 1 << num_Spins; // 2^num_Spins
    beta = vm["beta"].as<RealType>();
    spin_model = vm["spinmodel"].as<std::string>();

    // ========== general numerical parameters ==========
    seed = vm["seed"].as<std::string>();
    CET_therm_error = vm["CET_therm_error"].as<RealType>();
    CET_evol_error = vm["CET_evol_error"].as<RealType>();
    num_TimePoints = vm["numTimePoints"].as<uint>();
    dt = beta * RealType{0.5} / static_cast<RealType>(num_TimePoints - 1);
    num_Vectors_Per_Core = vm["numVectorsPerCore"].as<uint>();
    Gauss_covariance = vm["GaussCovariance"].as<RealType>();
    CET_rescale = vm["ChebyshevRescale"].as<RealType>();

    // ========== saving and naming ==========
    project_name = vm["project"].as<std::string>();
    filename_extension = vm["fileext"].as<std::string>();
    num_PrintDigits = vm["numPrintDigits"].as<size_t>();
}

void ParameterSpace::read_SpinSystem()
{
    // 1) open file and group:
    std::string total_filename = "Couplings/" + couplings_filename + ".hdf5";
    hid_t file_id = H5Fopen( total_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    hid_t group_id = H5Gopen( file_id, "/all", H5P_DEFAULT );

    // 2) read the number of spins:
    // a) open attribute and read type
    hid_t attr_id = H5Aopen( group_id, "num_Spins", H5P_DEFAULT );
    hid_t datatype_id = H5Aget_type( attr_id ); // datatype of the dataset -> uint

    // b) read attribute data
    H5Aread( attr_id, datatype_id, &num_Spins );

    // c) close attribute
    H5Tclose( datatype_id );
    H5Aclose( attr_id );

    // 3) read the coupling data:
    // a) open dataset and read type
    hid_t dataset = H5Dopen2( group_id, "/all/J_ij", H5P_DEFAULT );
    hid_t datatype = H5Dget_type( dataset ); // datatype of the dataset -> double

    // b) read linearized data from the dataset
    std::vector<double> linearized_data( num_Spins * num_Spins );
    H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized_data.data() );

    // c) close resources
    H5Tclose( datatype );
    H5Dclose( dataset );

    // 4) close resources
    H5Gclose( group_id );
    H5Fclose( file_id );

    // 5) tensorize the data:
    couplings.resize( num_Spins );
    for( size_t i = 0; i < num_Spins; i++ ) 
    {
        for( size_t j = 0; j < num_Spins; j++ ) 
        {
            couplings(i,j) = linearized_data[i*num_Spins + j] * rescale;
        }
    }
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
    << print::quantity_to_output_line( pre_colon_space, "rescale"            , print::remove_zeros(print::round_value_to_string(rescale,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "num_TimePoints" , std::to_string(num_TimePoints) ) 
    << print::quantity_to_output_line( pre_colon_space, "dt"       , print::remove_zeros(print::round_value_to_string(dt,num_PrintDigits)) ) 
    << print::quantity_to_output_line( pre_colon_space, "num_Vectors"   , std::to_string(num_Vectors_Per_Core*world_size) );
    return ss.str();
}

}