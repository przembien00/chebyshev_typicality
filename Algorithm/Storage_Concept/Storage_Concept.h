#pragma once

#include"../Types/Types.h"
#include"../Types/Tensors.h"
#include"../Types/Correlations.h"
#include"../cpp_libs/HDF5_Routines.h"
#ifdef USE_DOUBLE
static hid_t H5_REAL_TYPE = H5T_IEEE_F64LE;
#endif 
#ifdef USE_FLOAT 
static hid_t H5_REAL_TYPE = H5T_IEEE_F32LE;
#endif
#include"../Parameter_Space/Parameter_Space.h"


namespace Storage_Concept
{

using CorrelationTensor = Tensors::CorrelationTensor<Correlations::CorrelationVector>;

namespace ps = Parameter_Space;


// ================= CLASS DEFINITIONS =================
// HDF5 STORAGE CLASS
class HDF5_Storage
{
 public:
    HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace );

    void store_main( const ps::ParameterSpace& pspace, const CorrelationTensor& corr_Re, const CorrelationTensor& corr_Im );

    void finalize();

 private:
    const bool m_storing_permission{true}; // permission for the class instance to create folders/files and store data
    const size_t m_fname_max_length{};
    const size_t m_num_TriesToBuildFile{};
    hid_t m_file_id;
    hid_t m_results_group_id;
    std::string m_filename{};

    void create_folder_branch( const ps::ParameterSpace& pspace );
    void create_file( const ps::ParameterSpace& pspace );
    void store_correlation_tensor( const CorrelationTensor& CT, const hid_t group_id, const std::string dataset_name, const std::string dataset_info );
};


};