#pragma once

#include"../Types/Types.h"
#include"../Hamiltonians/Hamiltonians.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{

std::tuple<uint, uint> determine_CET_depth( const ham::Hamiltonian& H, const ps::ParameterSpace& pspace );

ComplexType cdot( const State& state1, const State& state2 );

size_t generate_seed( const ps::ParameterSpace& pspace, const size_t my_rank );

size_t throw_seed( const size_t seed, const size_t my_rank, const size_t sample );

States initialize_states( const ps::ParameterSpace& pspace, uint seed, uint sample );

void thermalize( const ham::Hamiltonian& H, States& states, const RealType beta, const int depth );

void evolve( const ham::Hamiltonian& H, States& states_L, States& states_R, const RealType dt, const int depth );

void compute_correlations_at( int t, States& v_psi_L, States& v_psi_R, CorrelationTensor& corrs );

State S_alpha_i_act( const State& state, const long site, char alpha );

States S_i_act( States& states, const long site );

RealType CET_coeff( int n, RealType x );

void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth );

void MPI_share_results( RealType& partition_function, CorrelationTensor& spin_c );

void normalize( RealType& partition_function, CorrelationTensor& spin_c );

}