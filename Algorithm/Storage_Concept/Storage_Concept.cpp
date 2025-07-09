#include"Storage_Concept.h"

#include<fstream>
#include<iostream>

#include"File_Management.h"
namespace fm = File_Management;

#include"STOC_Error_Handling.h"
namespace error = Storage_Concept::Error_Handling;

#include"Print_Routines.h"
namespace print = Print_Routines;

namespace Storage_Concept
{

namespace ps = Parameter_Space;

namespace hdf5r = HDF5_Routines;

// ================= CLASS IMPLEMENTATIONS =================
// HDF5 STORAGE CLASS
// constructor : create folder tree and data file
HDF5_Storage::HDF5_Storage( const ps::ParameterSpace& pspace ):
    m_storing_permission( true ),
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
    filename += "__" + pspace.couplings_filename;
    std::stringstream params_stream;
    params_stream << "__beta=" << pspace.beta << "__rescale=" << pspace.rescale;
    filename += params_stream.str(); 

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


void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, const CorrTen& corr )
{
    if( !m_storing_permission ){ return; } // permission request

    // ===== store parameters =====
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    hdf5r::store_scalar( ps_group_id, "num_Spins",                    pspace.num_Spins ); 
    hdf5r::store_scalar( ps_group_id, "num_HilbertSpaceDimension",    pspace.HilbertSpaceDimension ); 
    std::string tmp{ pspace.couplings_filename + ".hdf5" };
    hdf5r::store_string( ps_group_id, "src_file",                   tmp );
    hdf5r::store_string( ps_group_id, "spin_model",                 pspace.spin_model );

    //hdf5r::add_char_to_group( ps_group_id, "mean_symmetry_type",           pspace.spin_model.m_means._symmetry_type );
    //hdf5r::add_char_to_group( ps_group_id, "correlation_symmetry_type",    pspace.spin_model.m_correlations._symmetry_type );

    hdf5r::store_scalar( ps_group_id, "beta", pspace.beta );
    hdf5r::store_scalar( ps_group_id, "num_TimePoints",               pspace.num_TimePoints ); 
    hdf5r::store_scalar( ps_group_id, "delta_t",                  pspace.dt );

    hdf5r::store_string( ps_group_id, "original project_name",       pspace.project_name );

    H5Gclose( ps_group_id );

    // ===== store results =====
    m_results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    
    hdf5r::store_list( m_results_group_id, "g_zz", corr );

    H5Gclose( m_results_group_id );
}


void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}



};
