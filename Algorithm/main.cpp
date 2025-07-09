#include"main_header.h"



int main( const int argC, char* const argV[] )
{

// Initialize
ps::ParameterSpace my_pspace( argC, argV );
ham::Hamiltonian my_H( my_pspace );

CorrTen Correlations(my_pspace.num_TimePoints);
RealType Z = RealType{0.};

// Estimate the order of expansion needed to minimize the thermalization error.
RealType bound = my_H.a * my_pspace.beta * RealType(0.25) * std::exp(1.0);
uint depth_beta = static_cast<uint>(bound) + 20;
std::cout << "beta N_C > " << depth_beta << '\n';
// Estimate the order of expansion needed to minimize the evolution error.
bound = my_H.a * my_pspace.dt * RealType(0.5) * std::exp(1.0);
uint depth_dt = static_cast<uint>(bound) + 3;
std::cout << "dt N_C > " << depth_dt << '\n';

for( int k=0; k < my_pspace.num_Vectors; k++ )
{
    State psi_L = func::initialize_state( my_pspace ); // |psi_0>
    // Thermalize. Compute Z.

    func::CET( my_H, psi_L, my_pspace.beta * RealType{0.5}, depth_beta ); // e^(-beta*H/2)|psi_0>
    Z += std::pow( std::real(blaze::norm(psi_L)) , 2 ); // Z = <psi_0|e^(-beta*H)|psi_0>

    // Evolve: loop over times

    State psi_R = func::S_z_0_act( psi_L ); // S^z_0 e^(-beta*H/2)|psi_0>
    Correlations[0] += func::cdot( psi_L, func::S_z_0_act( psi_R ) );

    for( uint i = 1; i < my_pspace.num_TimePoints; i++ )
    {
        func::CET( my_H, psi_R, my_pspace.dt, depth_dt ); // e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
        func::CET( my_H, psi_L, - my_pspace.dt, depth_dt ); // e^(tau H) e^(-beta H/2)|psi_0>
        Correlations[i] += func::cdot( psi_L, func::S_z_0_act(psi_R) ); // <psi_0|e^(-beta H/2) e^(tau H) S^z_0 e^(-tau H) S^z_0 e^(-beta H/2)|psi_0>
    }
}



for( uint i = 0; i < Correlations.size(); i++ )
{
    Correlations[i] /= Z;
    std::cout << Correlations[i] << '\n';
}

// Store correlations
stor::HDF5_Storage my_data_storage( my_pspace );
my_data_storage.store_main( my_pspace, Correlations );
my_data_storage.finalize();

}

    // std::cout << "a=" << my_H.a << '\n';
    // State psi_test(my_pspace.HilbertSpaceDimension);
    // psi_test[3] = 1;
    // std::cout << "psi test" << '\n';
    // for( uint i= 0; i < psi_test.size(); i++)
    // {
    //     std::cout << psi_test[i] << '\n';
    // }
    // psi_test = my_H.act(psi_test);
    // std::cout << "psi test after H" << '\n';
    // for( uint i= 0; i < psi_test.size(); i++)
    // {
    //     std::cout << psi_test[i] << '\n';
    // }