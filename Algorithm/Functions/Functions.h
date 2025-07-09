#pragma once

#include"../Types/Types.h"
#include"../Hamiltonians/Hamiltonians.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace ps = Parameter_Space;
namespace ham = Hamiltonians;

namespace Functions
{

RealType cdot( const State& state1, const State& state2 );

State initialize_state( const ps::ParameterSpace& pspace );

State S_z_0_act( const State& state );

RealType CET_coeff( int n, RealType x );

void CET( ham::Hamiltonian& H, State& state, const RealType t, uint depth );

}