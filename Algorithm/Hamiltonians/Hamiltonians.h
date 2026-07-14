#pragma once
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"


namespace ps = Parameter_Space;

namespace Hamiltonians
{

struct Parameters
{
    RealType lambda{1.};
    RealType h_z{0.};
};

struct Bond
// One nonzero coupling J_ij (i < j), stored as single-spin bit masks and the flip-term prefactor
{
    long bit_i; // 1 << i
    long bit_j; // 1 << j
    RealType half_J; // J/2, prefactor of 0.5*(S_i^+ S_j^- + S_i^- S_j^+)
};

class Hamiltonian
{
    public:
    long numSpins;
    long dim;
    CouplingMatrix couplings;
    Parameters params;
    RealType E_max;
    RealType E_min;
    RealType a;
    RealType b;
    RealType b_over_a;
    std::vector<Bond> bonds; // nonzero off-diagonal couplings
    blaze::DynamicVector<RealType> diag; // diagonal: S^z S^z terms + field - spectrum shift

    // constructors
    Hamiltonian() = default;
    Hamiltonian( ps::ParameterSpace& pspace);

    State act( const State& state ) const;
    void act( const State& in, State& prev_inout ) const; // Chebyshev recurrence step: prev = 2*H*in - prev

    private:
    RealType zz_factor{1.}; // lambda for XXZ, 1 for ISO
    void build_tables( const RealType shift );
};

}
