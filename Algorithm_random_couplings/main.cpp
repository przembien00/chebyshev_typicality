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
size_t seed = func::generate_seed( my_pspace, my_rank );
CorrelationTensor averaged_correlations_R( my_pspace.symmetry_type, my_pspace.num_TimePoints );
CorrelationTensor averaged_correlations_I( my_pspace.symmetry_type, my_pspace.num_TimePoints );
CorrelationTensor sqsum_R( my_pspace.symmetry_type, my_pspace.num_TimePoints );
CorrelationTensor sqsum_I( my_pspace.symmetry_type, my_pspace.num_TimePoints );
my_clock.measure("Initialization");

// Distribute the disorder configurations across the MPI ranks: each rank handles its
// configs entirely on its own (couplings, Lanczos bandwidth, typicality sampling and the
// per-config normalization by Z); the config sums are reduced once after the loop.
// Each config is sampled with num_Vectors_Per_Core typicality vectors.
const uint num_local_configs = ( static_cast<uint>(my_rank) < my_pspace.num_Coupling_Configs ) ?
    ( my_pspace.num_Coupling_Configs - static_cast<uint>(my_rank) + static_cast<uint>(world_size) - 1 ) / static_cast<uint>(world_size) : 0;
tmm::Simple_Estimator my_estimator( my_rank, num_local_configs, "disorder sampling" );
my_estimator.enter_loop();
uint local_config = 0;
for( uint i = static_cast<uint>(my_rank); i < my_pspace.num_Coupling_Configs; i += static_cast<uint>(world_size) )
{
my_pspace.draw_couplings( seed, i );
ham::Hamiltonian my_H( my_pspace );
CorrelationTensor correlations_R{my_pspace.symmetry_type, my_pspace.num_TimePoints};
CorrelationTensor correlations_I{my_pspace.symmetry_type, my_pspace.num_TimePoints};
RealType Z = RealType{0.};
// Precompute the Chebyshev coefficient tables (adaptively truncated); the bandwidth
// a, b changes with each coupling configuration, so the tables are per-config
std::vector<ComplexType> coeffs_therm;
if( my_pspace.beta != RealType{0.} )
{
    coeffs_therm = func::CET_coefficients( -my_pspace.beta * RealType{0.5}, my_H.a, my_H.b, "imaginary", my_pspace.CET_therm_error );
}
const std::vector<ComplexType> coeffs_dt_plus = func::CET_coefficients( my_pspace.dt, my_H.a, my_H.b, my_pspace.evol_type, my_pspace.CET_evol_error );
const std::vector<ComplexType> coeffs_dt_minus = func::CET_coefficients( -my_pspace.dt, my_H.a, my_H.b, my_pspace.evol_type, my_pspace.CET_evol_error );
// mix the config index into the seed so every config is sampled with independent states
const size_t config_seed = func::throw_seed( seed, static_cast<size_t>(i), 0 );

// ====== Main loop ======
for( int k=0; k < my_pspace.num_Vectors_Per_Core; k++ )
{
    CorrelationTensor new_correlations_R(my_pspace.symmetry_type, my_pspace.num_TimePoints);
    CorrelationTensor new_correlations_I(my_pspace.symmetry_type, my_pspace.num_TimePoints);
    // Initialize and thermalize
    State psi_L;
    if( !func::initialize_state( my_pspace, config_seed, k, psi_L ) )
    {
        continue;
    }
    if(my_pspace.beta != RealType{0.})
    {
        func::CET( my_H, psi_L, coeffs_therm ); // e^(-beta*H/2)|psi_0>
    }
    Z += std::pow( std::real(blaze::norm(psi_L)) , 2 ); // Z = <psi_0|e^(-beta*H)|psi_0>
    // Evolve: loop over times

    States v_psi_R = func::S_i_act( my_pspace, psi_L, my_pspace.spin_site ); // S^a_0 e^(-beta*H/2)|psi_0>
    func::compute_correlations_at( 0, 0, my_pspace, psi_L, v_psi_R, new_correlations_R, new_correlations_I );
    for( uint t_point = 1; t_point < my_pspace.num_TimePoints; t_point++ )
    {
        std::for_each( v_psi_R.begin(), v_psi_R.end(), [&my_H, &coeffs_dt_minus]( State& psi_R )
        {
        func::CET( my_H, psi_R, coeffs_dt_minus ); // e^(-tau H) S^a_0 e^(-beta H/2)|psi_0>
        } );
        func::CET( my_H, psi_L, coeffs_dt_plus ); // e^(tau H) e^(-beta H/2)|psi_0>
        func::compute_correlations_at( t_point, 0, my_pspace, psi_L, v_psi_R, new_correlations_R, new_correlations_I ); // <psi_0|e^(-beta H/2) e^(tau H) S^a_0 e^(-tau H) S^b_0 e^(-beta H/2)|psi_0>
    }
    correlations_R += new_correlations_R;
    correlations_I += new_correlations_I;
}
// all typicality vectors of this config live on this rank: normalize locally
func::normalize( Z, correlations_R );
func::normalize( Z, correlations_I );
averaged_correlations_R += correlations_R;
averaged_correlations_I += correlations_I;
func::add_sqs( correlations_R, correlations_I, sqsum_R, sqsum_I );
my_estimator.estimate( local_config++ );
}
my_estimator.leave_loop();
// sum the per-config results over all ranks
func::MPI_share_results( averaged_correlations_R, averaged_correlations_I );
func::MPI_share_results( sqsum_R, sqsum_I );
RealType N = static_cast<RealType>(my_pspace.num_Coupling_Configs);

CorrelationTensor stds_R{my_pspace.symmetry_type, my_pspace.num_TimePoints};
CorrelationTensor stds_I{my_pspace.symmetry_type, my_pspace.num_TimePoints};
func::compute_stds( N, averaged_correlations_R, averaged_correlations_I, sqsum_R, sqsum_I, stds_R, stds_I );
func::normalize( N, averaged_correlations_R );
func::normalize( N, averaged_correlations_I );
my_clock.measure("Correlations");

// Store correlations
stor::HDF5_Storage my_data_storage( my_rank, my_pspace );
my_data_storage.store_main( my_pspace, averaged_correlations_R, averaged_correlations_I, stds_R, stds_I );
my_data_storage.finalize();
my_clock.measure("Saving");
my_clock.finalize();
MPI_Finalize();
print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n" );


}
