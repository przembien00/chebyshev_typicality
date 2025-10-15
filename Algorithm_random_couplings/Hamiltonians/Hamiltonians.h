#pragma once
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"


namespace ps = Parameter_Space;

namespace Hamiltonians
{

struct Parameters
{
    RealType lambda;
    RealType h_z;
};

class Hamiltonian
{
    public:
    long numSpins;
    long dim;
    RealType CET_rescale;
    CouplingMatrix couplings;
    Parameters params;
    RealType E_max;
    RealType E_min;
    RealType a;
    RealType b;
    RealType b_over_a;

    // constructors
    Hamiltonian() = default;
    Hamiltonian( ps::ParameterSpace& pspace);

    std::function<State(const State&)> act;
};

}