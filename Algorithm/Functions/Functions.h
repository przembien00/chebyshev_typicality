#pragma once

#include"../Types/Types.h"
#include"../Types/Tensors.h"
#include"../Types/Correlations.h"
#include"../Hamiltonians/Hamiltonians.h"
#include"../Parameter_Space/Parameter_Space.h"



namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{
using CorrelationTensor = Tensors::CorrelationTensor<Correlations::CorrelationVector>;
using CorrelationVector = Correlations::CorrelationVector;

std::tuple<uint, uint> determine_CET_depth( const ham::Hamiltonian& H, const ps::ParameterSpace& pspace );

ComplexType cdot( const State& state1, const State& state2 );

size_t generate_seed( const ps::ParameterSpace& pspace, const size_t my_rank );

size_t throw_seed( const size_t seed, const size_t my_rank, const size_t sample );

State initialize_state( const ps::ParameterSpace& pspace, uint seed, uint sample );

void compute_correlations_at( int t, long site, const ps::ParameterSpace& pspace, State& psi_L, States& v_psi_R, CorrelationTensor& corrs_Re, CorrelationTensor& corrs_Im );

State S_alpha_i_act( const State& state, const long site, char alpha );

States S_i_act( const ps::ParameterSpace& pspace, State& state, const long site );

void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth, std::string evol_type );

void MPI_share_results( RealType& partition_function, CorrelationTensor& correlations_Re, CorrelationTensor& correlations_Im );

void normalize( RealType& partition_function, CorrelationTensor& spin_c );

}