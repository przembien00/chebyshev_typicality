#include"main_header.h"



int main( const int argC, char* const argV[] )
{
// ====== Initialize MPI ======
// world_size is the number of cores, my_rank is the number of "this" core 
int world_size, my_rank;
MPI_Init( nullptr, nullptr );
MPI_Comm_size( MPI_COMM_WORLD, &world_size );
MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

// Initialize
ps::ParameterSpace my_pspace( argC, argV, world_size, my_rank );
ham::Hamiltonian my_H( my_pspace );
size_t seed = func::generate_seed( my_pspace, my_rank );
print::print_R0( my_rank, my_pspace.create_essentials_string() );

CorrTen Correlations(my_pspace.num_TimePoints);
RealType Z = RealType{0.};

// Estimate the order of expansion needed to minimize the thermalization error.
RealType bound = my_H.a * my_pspace.beta * RealType(0.25) * std::exp(1.0);
uint depth_beta = static_cast<uint>(bound) + 20;
std::stringstream ss;
ss << "Thermalization error = " << std::pow( bound/static_cast<RealType>(depth_beta), static_cast<RealType>(depth_beta) ) << '\n';
// Estimate the order of expansion needed to minimize the evolution error.
bound = my_H.a * my_pspace.dt * RealType(0.5) * std::exp(1.0);
uint depth_dt = static_cast<uint>(bound) + my_pspace.Chebyshev_cutoff;
ss.clear();
ss << "Time evolution error = " << std::pow( bound/static_cast<RealType>(depth_dt), static_cast<RealType>(depth_dt) ) << '\n';
print::print_R0( my_rank, ss.str() );
std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

for( int k=0; k < my_pspace.num_Vectors_Per_Core; k++ )
{
    if( k == 0 && my_rank == 0 )
    {
        begin = std::chrono::steady_clock::now();
    }
    State psi_L = func::initialize_state( my_pspace, seed, k ); // |psi_0>


    // Thermalize. Compute Z.
    // for( uint i=1; i < my_pspace.num_TimeSteps_therm; i++ )
    // {
    //     psi_L = func::RK4( my_H, psi_L, - my_pspace.dt );
    // }

    func::CET( my_H, psi_L, my_pspace.beta * RealType{0.5}, depth_beta ); // e^(-beta*H/2)|psi_0>
    Z += std::pow( std::real(blaze::norm(psi_L)) , 2 ); // Z = <psi_0|e^(-beta*H)|psi_0>

    // Evolve: loop over times

    State psi_R = func::S_z_0_act( psi_L ); // S^z_0 e^(-beta*H/2)|psi_0>
    Correlations[0] += func::cdot( psi_L, func::S_z_0_act( psi_R ) );

    for( uint i = 1; i < my_pspace.num_TimePoints; i++ )
    {
        // psi_R = func::RK4( my_H, psi_R, - my_pspace.dt );
        // psi_L = func::RK4( my_H, psi_L, my_pspace.dt );
        func::CET( my_H, psi_R, my_pspace.dt, depth_dt ); // e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
        func::CET( my_H, psi_L, - my_pspace.dt, depth_dt ); // e^(tau H) e^(-beta H/2)|psi_0>
        Correlations[i] += func::cdot( psi_L, func::S_z_0_act(psi_R) ); // <psi_0|e^(-beta H/2) e^(tau H) S^z_0 e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
    }
    if( k == 0 && my_rank == 0 )
    {
        end = std::chrono::steady_clock::now();
        std::cout << "Single vector duration = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
    }
}

func::MPI_share_results( Z, Correlations );
func::normalize( Z, Correlations );

// Store correlations
stor::HDF5_Storage my_data_storage( my_rank, my_pspace );
my_data_storage.store_main( my_pspace, Correlations );
my_data_storage.finalize();
MPI_Finalize();

}