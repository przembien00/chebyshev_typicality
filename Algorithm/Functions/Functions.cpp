#include"../Types/Types.h"
#include"../Hamiltonians/Hamiltonians.h"
#include<cmath>
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

State initialize_state( const ps::ParameterSpace& pspace )
// Get a state with random (Gaussian) complex coefficients.
{
    State state(pspace.HilbertSpaceDimension);
    std::random_device rd{}; 
    std::mt19937 gen{rd()}; 
    std::normal_distribution<float> d{0., pspace.Gauss_covariance}; 
    
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

    for( ulong ident = 0; ident < state.size(); ++ident )
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

RealType CET_coeff( int n, RealType x )
// Compute the nth coefficient of the Chebyshev polynomial expansion.
{
    if( n == 0 )
    {
        return std::cyl_bessel_i( RealType{0.0}, std::abs(x) );
    }
    else
    {
        RealType sign = (x>0. && n&1)? RealType{-1.0} : RealType{1.0}; // (-1)^n
        return 2 * sign * std::cyl_bessel_i( static_cast<RealType>(n), std::abs(x) );
    }
}

void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth )
// Apply e^(-tH) to the state using the Chebyshev expansion technique.
// The state is modified in place.
{
    State state_final = CET_coeff( 0, H.a * t ) * state; // |psi_0>
    State state_aux = H.act( state ); // |psi_1> = H|psi_0>
    state_final += CET_coeff( 1, H.a * t ) * state_aux;
   
    for( uint n=2; n < depth + 1; n++ )
    {
        state = 2 * H.act( state_aux ) - state; // |psi_n> = 2H|psi_n-1> - |psi_n-2>
        state_final += CET_coeff( n, H.a * t ) * state;
        std::swap( state, state_aux );
    }
    state = state_final;
}

}