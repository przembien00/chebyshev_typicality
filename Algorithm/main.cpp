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
CorrelationTensor correlations_R{my_pspace.symmetry_type, my_pspace.num_TimePoints};
CorrelationTensor correlations_I{my_pspace.symmetry_type, my_pspace.num_TimePoints};
RealType Z = RealType{0.};
auto [depth_beta, depth_dt] = func::determine_CET_depth( my_H, my_pspace );

my_clock.measure("Initialization");

tmm::Simple_Estimator my_estimator( my_rank, my_pspace.num_Vectors_Per_Core, "typicality sampling" );

my_estimator.enter_loop();
// ====== Main loop ======
for( int k=0; k < my_pspace.num_Vectors_Per_Core; k++ )
{
    CorrelationTensor new_correlations_R(my_pspace.symmetry_type, my_pspace.num_TimePoints);
    CorrelationTensor new_correlations_I(my_pspace.symmetry_type, my_pspace.num_TimePoints);
    // Initialize and thermalize
    State psi_L = func::initialize_state( my_pspace, seed, k ); // |psi_0>
    func::CET( my_H, psi_L, my_pspace.beta * RealType{0.5}, depth_beta ); // e^(-beta*H/2)|psi_0>
    Z += std::pow( std::real(blaze::norm(psi_L)) , 2 ); // Z = <psi_0|e^(-beta*H)|psi_0>
    // Evolve: loop over times

    States v_psi_R = func::S_i_act( my_pspace, psi_L, 0 ); // S^a_0 e^(-beta*H/2)|psi_0>
    func::compute_correlations_at( 0, 0, my_pspace, psi_L, v_psi_R, new_correlations_R, new_correlations_I );
    for( uint i = 1; i < my_pspace.num_TimePoints; i++ )
    {
        std::for_each( v_psi_R.begin(), v_psi_R.end(), [&my_H, &my_pspace, &depth_dt]( State& psi_R )
        {
        func::CET( my_H, psi_R, my_pspace.dt, depth_dt ); // e^(-tau H) S^a_0 e^(-beta H/2)|psi_0>
        } );
        func::CET( my_H, psi_L, - my_pspace.dt, depth_dt ); // e^(tau H) e^(-beta H/2)|psi_0>
        func::compute_correlations_at( i, 0, my_pspace, psi_L, v_psi_R, new_correlations_R, new_correlations_I ); // <psi_0|e^(-beta H/2) e^(tau H) S^a_0 e^(-tau H) S^b_0 e^(-beta H/2)|psi_0>
    }
    correlations_R += new_correlations_R;
    correlations_I += new_correlations_I;
    my_estimator.estimate(k);
}
my_estimator.leave_loop();
func::MPI_share_results( Z, correlations_R, correlations_I );
func::normalize( Z, correlations_R );
func::normalize( Z, correlations_I );
my_clock.measure("Correlations");

// Store correlations
stor::HDF5_Storage my_data_storage( my_rank, my_pspace );
my_data_storage.store_main( my_pspace, correlations_R, correlations_I );
my_data_storage.finalize();
my_clock.measure("Saving");
my_clock.finalize();
MPI_Finalize();
print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n" );


}