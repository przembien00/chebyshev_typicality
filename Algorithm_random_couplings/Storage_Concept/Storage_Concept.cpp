#include"Storage_Concept.h"

#include<fstream>
#include<iostream>

#include"../cpp_libs/File_Management.h"
namespace fm = File_Management;

#include"../cpp_libs/STOC_Error_Handling.h"
namespace error = Storage_Concept::Error_Handling;

#include"../cpp_libs/Print_Routines.h"
namespace print = Print_Routines;

namespace Storage_Concept
{

namespace ps = Parameter_Space;

namespace hdf5r = HDF5_Routines;

// ================= CLASS IMPLEMENTATIONS =================
// HDF5 STORAGE CLASS
// constructor : create folder tree and data file
HDF5_Storage::HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace ):
    m_storing_permission( my_rank == 0 ),
    m_fname_max_length( 200 ),
    m_num_TriesToBuildFile( 5 )
{
    create_folder_branch( pspace );
    create_file( pspace );
}
// create the folder branch in which the data will be stored
void HDF5_Storage::create_folder_branch( const ps::ParameterSpace& pspace )
{
    if( !m_storing_permission ){ return; } // permission request

    // determine folder branch list
    std::vector<std::string> folder_branch_list{};
    folder_branch_list.push_back( "Data" );
    if( pspace.project_name != "" )
    {
        folder_branch_list.push_back( pspace.project_name );
    }

    // create folder branch:
    m_filename = fm::create_folder_tree( folder_branch_list, m_fname_max_length );
}

// create the file in which the data will be stored
void HDF5_Storage::create_file( const ps::ParameterSpace& pspace )
{
    if( !m_storing_permission ){ return; } // permission request

    // create filename:
    std::string filename = pspace.spin_model; // spin model info
    filename += "__Random";
    std::stringstream params_stream;
    if( pspace.spin_site != 0 )
    {
        params_stream << "__site=" << pspace.spin_site;
    }
    params_stream << "__beta=" << pspace.beta;
    if (pspace.lambda != RealType{1.0})
    {
        params_stream << "__lambda=" << pspace.lambda;
    }
    if( pspace.h_z != RealType{0.} )
    {
        params_stream << "__h_z=" << pspace.h_z;
    }
    params_stream << "__numConfigs=" << pspace.num_Coupling_Configs;
    filename += params_stream.str(); 
    if( pspace.filename_extension != "" )
    {
        filename += "__" + pspace.filename_extension;
    }

    // create the file:
    print::cut_if_too_large( filename, m_fname_max_length );
    m_filename += filename;
    size_t count = 0;
    do{
        std::string tmp = m_filename + ".hdf5";
        std::ifstream f( tmp.c_str() );
        if( !f.good() ) // then the file doesn't exist yet
        {          
            m_file_id = H5Fcreate( tmp.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            if( count > 0 ) // inform about the difficulties -> creating file didn't work in the first try
            {
                std::cout << "\033[1;31mdata file already exists!\nsuccessfully created " << m_filename << " instead\033[0m\n";
            } 
            break;
        }
        else // file exists already 
        {
            m_filename += "X"; // adds an X to the filename and retries
        }
        f.close();
        count++;
    }while( count < m_num_TriesToBuildFile );
    m_filename += ".hdf5";
    if( count == m_num_TriesToBuildFile )
    {
        error::CREATE_FILE( m_filename, __PRETTY_FUNCTION__ );
    }
}


void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, const CorrelationTensor& corr_Re, const CorrelationTensor& corr_Im )
{
    if( !m_storing_permission ){ return; } // permission request

    // ===== store parameters =====
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    hdf5r::store_scalar( ps_group_id, "num_Spins",                    pspace.num_Spins ); 
    hdf5r::store_scalar( ps_group_id, "num_HilbertSpaceDimension",    pspace.HilbertSpaceDimension ); 
    hdf5r::store_string( ps_group_id, "spin_model",                 pspace.spin_model );
    hdf5r::store_scalar( ps_group_id, "spin_site",                  pspace.spin_site );
    hdf5r::store_string( ps_group_id, "evol_type",                  pspace.evol_type );
    hdf5r::store_scalar( ps_group_id, "Tmax",                       pspace.Tmax );

    hdf5r::store_scalar( ps_group_id, "beta", pspace.beta );
    hdf5r::store_scalar( ps_group_id, "h_z", pspace.h_z );
    hdf5r::store_scalar( ps_group_id, "num_TimePoints",               pspace.num_TimePoints ); 
    hdf5r::store_scalar( ps_group_id, "delta_t",                  pspace.dt );
    hdf5r::store_scalar( ps_group_id, "num_Cores", pspace.world_size);
    hdf5r::store_scalar( ps_group_id, "E_max", pspace.E_max );
    hdf5r::store_scalar( ps_group_id, "E_min", pspace.E_min );

    hdf5r::store_string( ps_group_id, "original project_name",       pspace.project_name );

    H5Gclose( ps_group_id );

    // ===== store results =====
    m_results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    
    store_correlation_tensor( corr_Re, m_results_group_id, "Re_correlation", "Real part of correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );
    store_correlation_tensor( corr_Im, m_results_group_id, "Im_correlation", "Imaginary part of correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );

    H5Gclose( m_results_group_id );
}

void HDF5_Storage::store_correlation_tensor( const CorrelationTensor& CT, const hid_t group_id, const std::string dataset_name, const std::string dataset_info )
{   
    if( !m_storing_permission ){ return; } // permission request

    hdf5r::store_2D_tensor<RealType>(group_id, dataset_name, H5_REAL_TYPE, CT, dataset_info);
}

void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}



};
