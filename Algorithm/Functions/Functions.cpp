#include"../Types/Types.h"
#include"../Hamiltonians/Hamiltonians.h"
#include<string>
#include<cmath>
#include <boost/math/special_functions/bessel.hpp>
#include<blaze/Math.h>
#include<random>
#include<iostream>

namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{

RealType cdot( const State& state1, const State& state2 )
// Complex scalar product of two states
{
    return std::real(blaze::dot( blaze::conj(state1), state2 ));
}

size_t generate_seed( const ps::ParameterSpace& pspace, const size_t my_rank )
{
    if( pspace.seed.find("random") != std::string::npos ) // random seed
    {
        size_t seed{};
        if( my_rank == 0 ) // draw seed on rank 0 
        {
            std::random_device my_seed{};
            seed = my_seed();
        }
        MPI_Bcast( &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); // broadcast seed to all the other ranks
        return seed;
    }  
    else // preset seed
    {
        return std::stoi( pspace.seed );
    }    
}

size_t throw_seed( const size_t seed, const size_t my_rank, const size_t sample )
{
    return (size_t) std::abs( static_cast<int>(seed) - static_cast<int>((sample + pow(10,5)) * (my_rank+1)) );
}


State initialize_state( const ps::ParameterSpace& pspace, uint seed, uint sample )
// Get a state with random (Gaussian) complex coefficients. The coeffs are drawn according to a seed
// specified in the parameter space (by deafult random).
{
    std::mt19937 gen{ static_cast<uint>( throw_seed( seed, pspace.my_rank, sample ) ) };
    State state(pspace.HilbertSpaceDimension);
    std::normal_distribution<RealType> d{0., pspace.Gauss_covariance}; 
    
    for( uint i = 0; i < state.size(); ++i )
    {
        RealType a = d(gen);
        RealType b = d(gen);

        state[i] = ComplexType{a, b};
    }
    
    return state;
}

State S_z_0_act( const State& state )
// Act on the state with the operator S_0^z.
{
    State new_state(state.size());

    for( unsigned long ident = 0; ident < state.size(); ++ident )
    {
        if( ident & 1L )
        {
            new_state[ident] = RealType{0.5} * state[ident];
        }
        else
        {
            new_state[ident] = - RealType{0.5} * state[ident];
        }

    }

    return new_state;
}

RealType CET_coeff( int n, RealType t, RealType a, RealType b )
// Compute the nth coefficient of the Chebyshev polynomial expansion.
{
    if( n == 0 )
    {
        return std::exp(b*t) * boost::math::cyl_bessel_i( RealType{0.0}, std::abs(a*t) );
    }
    else
    {
        RealType sign = (a*t>0. && n&1)? RealType{-1.0} : RealType{1.0}; // (-1)^n
        return 2 * sign * std::exp(b*t) * boost::math::cyl_bessel_i( static_cast<RealType>(n), std::abs(a*t) );
    }
}


void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth )
// Apply e^(-tH) to the state using the Chebyshev expansion technique.
// The state is modified in place.
{
    State state_final = CET_coeff( 0, t, H.a, H.b ) * state; // a_0|psi_0>
    State state_aux = H.act( state ); // |psi_1> = H|psi_0>
    state_final += CET_coeff( 1, t, H.a, H.b ) * state_aux; // a_0|psi_0> + a_1|psi_1>
   
    for( uint n=2; n < depth + 1; n++ )
    {
        state = 2 * H.act( state_aux ) - state; // |psi_n> = 2H|psi_n-1> - |psi_n-2>
        state_final += CET_coeff( n, t, H.a, H.b ) * state; // + a_n|psi_n> 
        std::swap( state, state_aux );
    }
    state = state_final;
}

State RK4( ham::Hamiltonian& H, State& state, const RealType dt )
{
    State state_out = state;
    for( uint n = 1; n < 5; n++ )
    {
        state = H.act( state ) * dt / static_cast<RealType>(n);
        state_out += state;
    }
    return state_out;
}

// sum the results of all cores and broadcast the sum to all cores with MPI_Allreduce 
void MPI_share_results( RealType& partition_function, CorrTen& spin_c )
{
    // share correlation results
    std::vector<RealType> rcv_buf( spin_c.size() ); 
    MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    spin_c = rcv_buf;

    // share partition function results
    std::vector<RealType> send_buf = { partition_function };
    std::vector<RealType> receive_buf( 1 );
    MPI_Allreduce( send_buf.data(), receive_buf.data(), 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    partition_function = receive_buf.at(0);
}

void normalize( RealType& partition_function, CorrTen& spin_c )
{
    for( uint i = 0; i < spin_c.size(); i++ )
    {
        spin_c[i] /= partition_function;
    }
}

}