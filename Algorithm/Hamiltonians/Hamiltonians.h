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
    long coord_num;
    CouplingMatrix couplings;
    Parameters params;
    RealType a;

    // constructors
    Hamiltonian() = default;
    Hamiltonian(const ps::ParameterSpace& pspace);

    std::function<State(const State&)> act;
};

}