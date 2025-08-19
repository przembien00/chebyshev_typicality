#pragma once

#include"../Types/Types.h"
#include"../Hamiltonians/Hamiltonians.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{

std::tuple<uint, uint> determine_CET_depth( const ham::Hamiltonian& H, const ps::ParameterSpace& pspace );

RealType cdot( const State& state1, const State& state2 );

size_t generate_seed( const ps::ParameterSpace& pspace, const size_t my_rank );

size_t throw_seed( const size_t seed, const size_t my_rank, const size_t sample );

State initialize_state( const ps::ParameterSpace& pspace, uint seed, uint sample );

State S_z_i_act( const State& state, const unsigned long site );

RealType CET_coeff( int n, RealType x );

void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth );

void MPI_share_results( RealType& partition_function, CorrTen& spin_c );

void normalize( RealType& partition_function, CorrTen& spin_c );

}