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
CorrTen Correlations(my_pspace.num_TimePoints);
RealType Z = RealType{0.};

auto [depth_beta, depth_dt] = func::determine_CET_depth( my_H, my_pspace );

my_clock.measure("Initialization");

tmm::Simple_Estimator my_estimator( my_rank, my_pspace.num_Vectors_Per_Core, "typicality sampling" );

my_estimator.enter_loop();
// ====== Main loop ======
for( int k=0; k < my_pspace.num_Vectors_Per_Core; k++ )
{
    State psi_L = func::initialize_state( my_pspace, seed, k ); // |psi_0>

    func::CET( my_H, psi_L, my_pspace.beta * RealType{0.5}, depth_beta ); // e^(-beta*H/2)|psi_0>
    Z += std::pow( std::real(blaze::norm(psi_L)) , 2 ); // Z = <psi_0|e^(-beta*H)|psi_0>

    // Evolve: loop over times

    State psi_R = func::S_z_i_act( psi_L, 0 ); // S^z_0 e^(-beta*H/2)|psi_0>
    Correlations[0] += func::cdot( psi_L, func::S_z_i_act( psi_R, 0 ) );

    for( uint i = 1; i < my_pspace.num_TimePoints; i++ )
    {
        // psi_R = func::RK4( my_H, psi_R, - my_pspace.dt );
        // psi_L = func::RK4( my_H, psi_L, my_pspace.dt );
        func::CET( my_H, psi_R, my_pspace.dt, depth_dt ); // e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
        func::CET( my_H, psi_L, - my_pspace.dt, depth_dt ); // e^(tau H) e^(-beta H/2)|psi_0>
        Correlations[i] += func::cdot( psi_L, func::S_z_i_act(psi_R, 0) ); // <psi_0|e^(-beta H/2) e^(tau H) S^z_0 e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
    }
    my_estimator.estimate(k);
}
my_estimator.leave_loop();

func::MPI_share_results( Z, Correlations );
func::normalize( Z, Correlations );
my_clock.measure("Correlations");

// Store correlations
stor::HDF5_Storage my_data_storage( my_rank, my_pspace );
my_data_storage.store_main( my_pspace, Correlations );
my_data_storage.finalize();
my_clock.measure("Saving");
my_clock.finalize();
MPI_Finalize();
print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++\n" );


}