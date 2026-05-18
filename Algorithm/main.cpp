#include"main_header.h"



int main( const int argC, char* const argV[] )
{


// ====== Initialize MPI ======
// world_size is the number of cores, my_rank is the number of "this" core 
int world_size, my_rank;
MPI_Init( nullptr, nullptr );
MPI_Comm_size( MPI_COMM_WORLD, &world_size );
MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

tmm::Clock my_clock( my_rank );

// Initialize
ps::ParameterSpace my_pspace( argC, argV, world_size, my_rank );
print::print_R0( my_rank, my_pspace.create_essentials_string() );
print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n" );
ham::Hamiltonian my_H( my_pspace );
size_t seed = func::generate_seed( my_pspace, my_rank );
const size_t num_sites = my_pspace.spin_sites.size();
std::vector<CorrelationTensor> correlations_R( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> correlations_I( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> correlations_R_sq( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> correlations_I_sq( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> covariances_Z_R( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> covariances_Z_I( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
RealType Z = RealType{0.};
RealType Z_sq = RealType{0.};
auto [depth_beta, depth_dt] = func::determine_CET_depth( my_H, my_pspace );
my_clock.measure("Initialization");

tmm::Simple_Estimator my_estimator( my_rank, my_pspace.num_Vectors_Per_Core, "typicality sampling" );

my_estimator.enter_loop();
// ====== Main loop ======
for( int k=0; k < my_pspace.num_Vectors_Per_Core; k++ )
{
    std::vector<CorrelationTensor> new_correlations_R( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
    std::vector<CorrelationTensor> new_correlations_I( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
    std::vector<CorrelationTensor> new_correlations_R_sq( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
    std::vector<CorrelationTensor> new_correlations_I_sq( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
    // Initialize and thermalize
    State psi_L;
    if( !func::initialize_state( my_pspace, seed, k, psi_L ) )
    {
        continue;
    }
    if(my_pspace.beta != RealType{0.})
    {
        func::CET( my_H, psi_L, -my_pspace.beta * RealType{0.5}, depth_beta, "imaginary" ); // e^(-beta*H/2)|psi_0>
    }
    RealType Z_sample  = std::pow( std::real(blaze::norm(psi_L)) , 2 );
    Z += Z_sample; // Z = <psi_0|e^(-beta*H)|psi_0>
    Z_sq += std::pow( Z_sample , 2 ); // for stds
    // Evolve: loop over times

    States v_psi_R = func::S_i_act( my_pspace, psi_L, 0 ); // S^a_0 e^(-beta*H/2)|psi_0>
    func::compute_correlations_at_sites( 0, my_pspace, psi_L, v_psi_R, my_pspace.spin_sites, new_correlations_R, new_correlations_I, new_correlations_R_sq, new_correlations_I_sq );
    for( uint t_point = 1; t_point < my_pspace.num_TimePoints; t_point++ )
    {
        std::for_each( v_psi_R.begin(), v_psi_R.end(), [&my_H, &my_pspace, &depth_dt]( State& psi_R )
        {
        func::CET( my_H, psi_R, -my_pspace.dt, depth_dt, my_pspace.evol_type ); // e^(-tau H) S^a_0 e^(-beta H/2)|psi_0>
        } );
        func::CET( my_H, psi_L, my_pspace.dt, depth_dt, my_pspace.evol_type ); // e^(tau H) e^(-beta H/2)|psi_0>
        func::compute_correlations_at_sites( t_point, my_pspace, psi_L, v_psi_R, my_pspace.spin_sites, new_correlations_R, new_correlations_I, new_correlations_R_sq, new_correlations_I_sq ); // <psi_0|e^(-beta H/2) e^(tau H) S^a_i e^(-tau H) S^b_0 e^(-beta H/2)|psi_0>
    }
    for( size_t site_idx = 0; site_idx < num_sites; ++site_idx )
    {
        correlations_R[site_idx] += new_correlations_R[site_idx];
        correlations_I[site_idx] += new_correlations_I[site_idx];
        covariances_Z_R[site_idx] += Z_sample * new_correlations_R[site_idx];
        covariances_Z_I[site_idx] += Z_sample * new_correlations_I[site_idx];
        correlations_R_sq[site_idx] += new_correlations_R_sq[site_idx];
        correlations_I_sq[site_idx] += new_correlations_I_sq[site_idx];
    }
    my_estimator.estimate(k);
}
my_estimator.leave_loop();
if( num_sites > 0 )
{
    func::MPI_share_results( Z, Z_sq, correlations_R[0], correlations_I[0], correlations_R_sq[0], correlations_I_sq[0], covariances_Z_R[0], covariances_Z_I[0] );
    for( size_t site_idx = 1; site_idx < num_sites; ++site_idx )
    {
        func::MPI_share_results( correlations_R[site_idx], correlations_I[site_idx], correlations_R_sq[site_idx], correlations_I_sq[site_idx], covariances_Z_R[site_idx], covariances_Z_I[site_idx] );
    }
}
for( size_t site_idx = 0; site_idx < num_sites; ++site_idx )
{
    func::normalize( Z, correlations_R[site_idx] );
    func::normalize( Z, correlations_I[site_idx] );
}
std::vector<CorrelationTensor> stds_R( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
std::vector<CorrelationTensor> stds_I( num_sites, CorrelationTensor{my_pspace.symmetry_type, my_pspace.num_TimePoints} );
for( size_t site_idx = 0; site_idx < num_sites; ++site_idx )
{
    func::compute_stds( my_pspace, Z, Z_sq, correlations_R[site_idx], correlations_I[site_idx], correlations_R_sq[site_idx], correlations_I_sq[site_idx], covariances_Z_R[site_idx], covariances_Z_I[site_idx], stds_R[site_idx], stds_I[site_idx] );
}
my_clock.measure("Correlations");

// Store correlations
stor::HDF5_Storage my_data_storage( my_rank, my_pspace );
my_data_storage.store_main( my_pspace, correlations_R, correlations_I, stds_R, stds_I );
my_clock.measure("Saving");
my_clock.finalize();
my_data_storage.store_runtime( my_clock );
my_data_storage.finalize();
MPI_Finalize();
print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n" );


}
