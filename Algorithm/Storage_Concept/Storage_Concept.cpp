#include"Storage_Concept.h"

#include<fstream>
#include<iostream>
#include<sstream>
#include<stdexcept>

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

namespace
{

std::string spin_sites_to_filename_fragment( const std::vector<uint>& sites )
{
    std::stringstream ss;
    for( size_t i = 0; i < sites.size(); ++i )
    {
        if( i > 0 )
        {
            ss << "-";
        }
        ss << sites[i];
    }
    return ss.str();
}

std::string site_group_name( const uint site )
{
    return std::to_string( site ) + "-0";
}

}

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
    filename += "__" + pspace.couplings_filename;
    std::stringstream params_stream;
    if( !( pspace.spin_sites.size() == 1 && pspace.spin_sites.front() == 0 ) )
    {
        params_stream << "__sites=" << spin_sites_to_filename_fragment( pspace.spin_sites );
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
    if( pspace.rescale != RealType{1.} )
    {
        params_stream << "__rescale=" << pspace.rescale;
    }
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


void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, const std::vector<CorrelationTensor>& corr_Re, const std::vector<CorrelationTensor>& corr_Im, const std::vector<CorrelationTensor>& stds_Re, const std::vector<CorrelationTensor>& stds_Im )
{
    if( !m_storing_permission ){ return; } // permission request
    if( corr_Re.size() != pspace.spin_sites.size() || corr_Im.size() != pspace.spin_sites.size() || stds_Re.size() != pspace.spin_sites.size() || stds_Im.size() != pspace.spin_sites.size() )
    {
        throw std::runtime_error( "site result arrays do not match the configured spin_sites list" );
    }

    // ===== store parameters =====
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    hdf5r::store_scalar( ps_group_id, "num_Spins",                    pspace.num_Spins ); 
    hdf5r::store_scalar( ps_group_id, "num_HilbertSpaceDimension",    pspace.HilbertSpaceDimension ); 
    std::string tmp{ pspace.couplings_filename + ".hdf5" };
    hdf5r::store_string( ps_group_id, "src_file",                   tmp );
    hdf5r::store_string( ps_group_id, "spin_model",                 pspace.spin_model );
    hdf5r::store_scalar( ps_group_id, "rescale",                    pspace.rescale );
    std::stringstream spin_sites_stream;
    for( size_t i = 0; i < pspace.spin_sites.size(); ++i )
    {
        if( i > 0 )
        {
            spin_sites_stream << ",";
        }
        spin_sites_stream << pspace.spin_sites[i];
    }
    hdf5r::store_string( ps_group_id, "spin_sites",                 spin_sites_stream.str() );
    hdf5r::store_string( ps_group_id, "evol_type",                  pspace.evol_type );
    hdf5r::store_scalar( ps_group_id, "Tmax",                       pspace.Tmax );

    hdf5r::store_scalar( ps_group_id, "beta", pspace.beta );
    hdf5r::store_scalar( ps_group_id, "h_z", pspace.h_z );
    hdf5r::store_scalar( ps_group_id, "num_TimePoints",               pspace.num_TimePoints ); 
    hdf5r::store_scalar( ps_group_id, "delta_t",                  pspace.dt );
    hdf5r::store_scalar( ps_group_id, "num_Cores", pspace.world_size);
    hdf5r::store_scalar( ps_group_id, "num_Vectors_Per_Core",    pspace.num_Vectors_Per_Core );
    hdf5r::store_scalar( ps_group_id, "E_max", pspace.E_max );
    hdf5r::store_scalar( ps_group_id, "E_min", pspace.E_min );

    hdf5r::store_string( ps_group_id, "original project_name",       pspace.project_name );

    H5Gclose( ps_group_id );

    // ===== store results =====
    m_results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hid_t re_group = H5Gcreate( m_results_group_id, "Re_correlation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t im_group = H5Gcreate( m_results_group_id, "Im_correlation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t re_stds_group = H5Gcreate( m_results_group_id, "Re_stds", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hid_t im_stds_group = H5Gcreate( m_results_group_id, "Im_stds", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    for( size_t site_idx = 0; site_idx < pspace.spin_sites.size(); ++site_idx )
    {
        const std::string group_name = site_group_name( pspace.spin_sites[site_idx] );
        store_correlation_tensor( corr_Re[site_idx], re_group, group_name, "Real part of correlations <S^a_i(t)S^b_0(0)> for a given site" );
        store_correlation_tensor( corr_Im[site_idx], im_group, group_name, "Imaginary part of correlations <S^a_i(t)S^b_0(0)> for a given site" );
        store_correlation_tensor( stds_Re[site_idx], re_stds_group, group_name, "Standard deviations of the real part of correlations for a given site" );
        store_correlation_tensor( stds_Im[site_idx], im_stds_group, group_name, "Standard deviations of the imaginary part of correlations for a given site" );
    }

    H5Gclose( re_group );
    H5Gclose( im_group );
    H5Gclose( re_stds_group );
    H5Gclose( im_stds_group );

    H5Gclose( m_results_group_id );
}

void HDF5_Storage::store_correlation_tensor( const CorrelationTensor& CT, const hid_t group_id, const std::string dataset_name, const std::string dataset_info )
{   
    if( !m_storing_permission ){ return; } // permission request

    hdf5r::store_2D_tensor<RealType>(group_id, dataset_name, H5_REAL_TYPE, CT, dataset_info);
}

void HDF5_Storage::store_runtime( const Time_Measurement::Clock& clock )
{
    if( !m_storing_permission ){ return; }

    auto group_id = H5Gcreate( m_file_id, "runtime_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    for( const auto& [name, secs] : clock.get_measurements() )
    {
        hdf5r::store_scalar( group_id, name + "_s", secs );
    }
    hdf5r::store_scalar( group_id, "total_s", clock.get_total_s() );
    H5Gclose( group_id );
}

void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}



};
