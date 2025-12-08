#include"../Types/Types.h"
#include"../Types/Tensors.h"
#include"../Types/Correlations.h"
#include"../Hamiltonians/Hamiltonians.h"
#include<string>
#include<cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include<blaze/Math.h>
#include<random>
#include<iostream>



namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{
using CorrelationTensor = Tensors::CorrelationTensor<Correlations::CorrelationVector>;
using CorrelationVector = Correlations::CorrelationVector;

ComplexType cdot( const State& state1, const State& state2 )
// Complex scalar product of two states
{
    return blaze::dot( blaze::conj(state1), state2 );
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
{
    State state(pspace.HilbertSpaceDimension);
    if( pspace.full_diagonalization )
    {
        size_t index = sample + pspace.my_rank * pspace.num_Vectors_Per_Core;
        if( index < pspace.HilbertSpaceDimension )
        {
            state[index] = ComplexType{1.,0.};
        }
    }
    else
    {
        std::mt19937 gen{ static_cast<uint>( throw_seed( seed, pspace.my_rank, sample ) ) };
        std::normal_distribution<RealType> d{0., pspace.Gauss_covariance};
        for( uint i = 0; i < state.size(); ++i )
        {
            RealType a = d(gen);
            RealType b = d(gen);

            state[i] = ComplexType{a, b};
        }
    }
    return state;
}

State S_alpha_i_act( const State& state, const long site, char alpha )
// Act on the state with the operator S_i^alpha, alpha=x,y,z.
{
    State new_state(state.size());
    switch( alpha )
    {
        case 'x':
        {
            for( long ident = 0; ident < state.size(); ++ident )
            {
                long new_ident = ident ^ ( 1L << site );
                new_state[new_ident] += RealType{0.5} * state[ident];
            }
            return new_state;
        }
        case 'y':
        {
            for( long ident = 0; ident < state.size(); ++ident )
            {
                long new_ident = ident ^ ( 1L << site );
                if( ident >> site & 1L )
                {
                    new_state[new_ident] += ComplexType{0.,0.5} * state[ident]; // -1/(2i)
                }
                else
                {
                    new_state[new_ident] += - ComplexType{0.,0.5} * state[ident]; // 1/(2i)
                }
            }
            return new_state;
        }
        case 'z':
        {
            for( long ident = 0; ident < state.size(); ++ident )
            {
                if( ident >> site & 1L )
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
    }
}

void compute_correlations_at( int t, long site, const ps::ParameterSpace& pspace, State& psi_L, States& v_psi_R, CorrelationTensor& corrs_Re, CorrelationTensor& corrs_Im )
{
    ComplexType c;
    switch( pspace.symmetry_type )
    {
        case 'A':
        {
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'z' ) );
            corrs_Re(2,2)[t] = std::real(c);
            corrs_Im(2,2)[t] = std::imag(c);
            break;
        }
        case 'B':
        {
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'x' ) );
            corrs_Re(0,0)[t] = std::real(c);
            corrs_Im(0,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[1], site, 'z' ) );
            corrs_Re(2,2)[t] = std::real(c);
            corrs_Im(2,2)[t] = std::imag(c);
            break;
        }
        case 'C':
        {
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'x' ) );
            corrs_Re(0,0)[t] = std::real(c);
            corrs_Im(0,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'y' ) );
            corrs_Re(1,0)[t] = std::real(c);
            corrs_Im(1,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[1], site, 'x' ) );
            corrs_Re(0,1)[t] = std::real(c);
            corrs_Im(0,1)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[2], site, 'z' ) );
            corrs_Re(2,2)[t] = std::real(c);
            corrs_Im(2,2)[t] = std::imag(c);
            break;
        }
        case 'D':
        {
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'x' ) );
            corrs_Re(0,0)[t] = std::real(c);
            corrs_Im(0,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[1], site, 'x' ) );
            corrs_Re(0,1)[t] = std::real(c);
            corrs_Im(0,1)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[2], site, 'x' ) );
            corrs_Re(0,2)[t] = std::real(c);
            corrs_Im(0,2)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'y' ) );
            corrs_Re(1,0)[t] = std::real(c);
            corrs_Im(1,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[1], site, 'y' ) );
            corrs_Re(1,1)[t] = std::real(c);
            corrs_Im(1,1)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[2], site, 'y' ) );
            corrs_Re(1,2)[t] = std::real(c);
            corrs_Im(1,2)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[0], site, 'z' ) );
            corrs_Re(2,0)[t] = std::real(c);
            corrs_Im(2,0)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[1], site, 'z' ) );
            corrs_Re(2,1)[t] = std::real(c);
            corrs_Im(2,1)[t] = std::imag(c);
            c = cdot( psi_L, S_alpha_i_act( v_psi_R[2], site, 'z' ) );
            corrs_Re(2,2)[t] = std::real(c);
            corrs_Im(2,2)[t] = std::imag(c);
            break;
        }
    }
}

std::tuple<RealType, RealType> correlation( State& psi_L, State& psi_R, long site, char alpha )
{
    ComplexType c = cdot( psi_L, S_alpha_i_act( psi_R, site, alpha ) );
    return std::make_tuple( std::real(c), std::imag(c) );
}

States S_i_act( const ps::ParameterSpace& pspace, State& state, const long site )
{
    States states;
    switch( pspace.symmetry_type )
    {
        case 'A':
        {
            states.emplace_back(S_alpha_i_act(state,site,'z'));     
            return states;
        }
        case 'B':
        {
            states.emplace_back(S_alpha_i_act(state,site,'x'));
            states.emplace_back(S_alpha_i_act(state,site,'z'));
            return states;
        }
        case 'C':
        {
            states.emplace_back(S_alpha_i_act(state,site,'x'));
            states.emplace_back(S_alpha_i_act(state,site,'y'));
            states.emplace_back(S_alpha_i_act(state,site,'z'));
            return states;
        }
        case 'D':
        {
            states.emplace_back(S_alpha_i_act(state,site,'x'));
            states.emplace_back(S_alpha_i_act(state,site,'y'));
            states.emplace_back(S_alpha_i_act(state,site,'z'));
            return states;
        }
    }
}


std::tuple<uint, uint> determine_CET_depth( const ham::Hamiltonian& H, const ps::ParameterSpace& pspace )
// Determine the depth of the Chebyshev expansion needed to minimize the thermalization and evolution errors
{
    RealType factor_therm = H.a * pspace.beta * RealType{0.25} * std::exp(1.0);
    RealType factor_evol = H.a * pspace.dt * RealType{0.5} * std::exp(1.0);
    if( factor_therm == RealType{0.} )
    {
        return std::make_tuple( 0, 
         static_cast<uint>(std::ceil( - std::log( pspace.CET_evol_error ) / boost::math::lambert_w0( - std::log( pspace.CET_evol_error ) / factor_evol ) ) ) );
    }
    else
    {
        return std::make_tuple(
        static_cast<uint>(std::ceil( - std::log( pspace.CET_therm_error ) / boost::math::lambert_w0( - std::log( pspace.CET_therm_error ) / factor_therm ) ) ),
        static_cast<uint>(std::ceil( - std::log( pspace.CET_evol_error ) / boost::math::lambert_w0( - std::log( pspace.CET_evol_error ) / factor_evol ) ) ) ); 
    }
   }

ComplexType CET_coeff( int n, RealType t, RealType a, RealType b, std::string evol_type )
// Compute the nth coefficient of the Chebyshev polynomial expansion.
{   
    if( evol_type == "imaginary" )
    {
    if( n == 0 )
    {
        return std::exp(b*t) * boost::math::cyl_bessel_i( RealType{0.0}, std::abs(a*t) );
    }
    else
    {
        RealType sign = (a*t<0. && n&1)? RealType{-1.0} : RealType{1.0};
        return RealType{2.0} * sign * std::exp(b*t) * boost::math::cyl_bessel_i( static_cast<RealType>(n), std::abs(a*t) );
    }
    }
    else
    {
    if( n == 0 )
    {
        return std::exp(ComplexType{0.,b*std::abs(t)}) * boost::math::cyl_bessel_j( RealType{0.0}, std::abs(a*t) );
    }
    else
    {
        // RealType sign = (a*t<0. && n&1)? RealType{-1.0} : RealType{1.0};
        return RealType{2.0} * pow(ComplexType{0.,1.}, n) * std::exp(ComplexType{0.,b*std::abs(t)}) * boost::math::cyl_bessel_j( static_cast<RealType>(n), std::abs(a*t) );
    }
    }
}
void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth, std::string evol_type )
// Apply e^(tH)/e^(-itH) to the state using the Chebyshev expansion technique.
// The state is modified in place.
{
    State state_final = CET_coeff( 0, t, H.a, H.b, evol_type ) * state; // a_0|psi_0>
    State state_aux = H.act( state ); // |psi_1> = H|psi_0>
    state_final += CET_coeff( 1, t, H.a, H.b, evol_type ) * state_aux; // a_0|psi_0> + a_1|psi_1>
   
    for( uint n=2; n < depth + 1; n++ )
    {
        state = 2 * H.act( state_aux ) - state; // |psi_n> = 2H|psi_n-1> - |psi_n-2>
        state_final += CET_coeff( n, t, H.a, H.b, evol_type ) * state; // + a_n|psi_n> 
        std::swap( state, state_aux );
    }
    state = state_final;
}

void MPI_share_results( RealType& partition_function, CorrelationTensor& correlations_Re, CorrelationTensor& correlations_Im )
// sum the results of all cores and broadcast the sum to all cores with MPI_Allreduce
{
    // share correlation results
    std::for_each( correlations_Re.begin(), correlations_Re.end(), []( CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    std::for_each( correlations_Im.begin(), correlations_Im.end(), []( CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share partition function results
    std::vector<RealType> send_buf = { partition_function };
    std::vector<RealType> receive_buf( 1 );
    MPI_Allreduce( send_buf.data(), receive_buf.data(), 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    partition_function = receive_buf.at(0);
}

void normalize( RealType& partition_function, CorrelationTensor& correlations )
{
    std::for_each( correlations.begin(), correlations.end(), [&partition_function]( CorrelationVector& spin_c ) 
    {
        spin_c *= RealType{1.0}/partition_function; 
    } );
}

void add_sqs( const CorrelationTensor& correlations_Re, const CorrelationTensor& correlations_Im, CorrelationTensor& sqsum_Re, CorrelationTensor& sqsum_Im )
{
    std::transform( correlations_Re.cbegin(), correlations_Re.cend(), sqsum_Re.begin(), sqsum_Re.begin(), [&]( const CorrelationVector& sample_C, CorrelationVector& sample_sqsum_C  )
    {
        CorrelationVector sqsum_C( sample_sqsum_C.size() );
        std::transform( sample_C.cbegin(), sample_C.cend(), sample_sqsum_C.begin(), sqsum_C.begin(), [&]( const auto& sample, auto& sample_sqsum )
        {
            return sample_sqsum + std::pow( sample, 2 );
        } );
        return sqsum_C;
    } );

    std::transform( correlations_Im.cbegin(), correlations_Im.cend(), sqsum_Im.begin(), sqsum_Im.begin(), [&]( const CorrelationVector& sample_C, CorrelationVector& sample_sqsum_C  )
    {
        CorrelationVector sqsum_C( sample_sqsum_C.size() );
        std::transform( sample_C.cbegin(), sample_C.cend(), sample_sqsum_C.cbegin(), sqsum_C.begin(), [&]( const auto& sample, auto& sample_sqsum )
        {
            return sample_sqsum + std::pow( sample, 2 );
        } );
        return sqsum_C;
    } );
}

void compute_stds( RealType M, const CorrelationTensor& correlations_Re, const CorrelationTensor& correlations_Im, const CorrelationTensor& sqsum_Re, const CorrelationTensor& sqsum_Im, CorrelationTensor& stds_Re, CorrelationTensor& stds_Im )
{
    std::transform( correlations_Re.cbegin(), correlations_Re.cend(), sqsum_Re.cbegin(), stds_Re.begin(), [&]( const CorrelationVector& sample_sum_C, const CorrelationVector& sample_sqsum_C )
    {
        CorrelationVector std_C( sample_sum_C.size() );
        std::transform( sample_sum_C.cbegin(), sample_sum_C.cend(), sample_sqsum_C.cbegin(), std_C.begin(), [&]( const auto& sample_sum, const auto& sample_sqsum )
        {
            return std::sqrt(std::abs(  sample_sqsum / M - std::pow( sample_sum / M, 2 )  ));
        } );
        return std_C;
    } );

    std::transform( correlations_Im.cbegin(), correlations_Im.cend(), sqsum_Im.cbegin(), stds_Im.begin(), [&]( const CorrelationVector& sample_sum_C, const CorrelationVector& sample_sqsum_C )
    {
        CorrelationVector std_C( sample_sum_C.size() );
        std::transform( sample_sum_C.cbegin(), sample_sum_C.cend(), sample_sqsum_C.cbegin(), std_C.begin(), [&]( const auto& sample_sum, const auto& sample_sqsum )
        {
            return std::sqrt(std::abs(  sample_sqsum / M - std::pow( sample_sum / M, 2 )  ));
        } );
        return std_C;
    } );
}

}