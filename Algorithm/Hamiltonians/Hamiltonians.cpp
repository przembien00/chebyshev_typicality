#include"Hamiltonians.h"
#include<functional>
#include<algorithm>
#include<cmath>
#include<bit>
#include<string>
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"
#include<iostream>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <lambda_lanczos/lambda_lanczos_tridiagonal_impl.hpp>

namespace ps = Parameter_Space;

namespace Hamiltonians
{

State Hamiltonian::act( const State& state ) const
// Apply the (shifted) Hamiltonian using the precomputed diagonal and bond list.
// The full diagonal (S^z S^z terms, field, spectrum shift) is one fused pass;
// the spin-flip terms are one streaming pass per bond, in gather form.
{
    State new_state( dim );
    for( long ident = 0; ident < dim; ++ident )
    {
        new_state[ident] = diag[ident] * state[ident];
    }
    for( const Bond& bond : bonds )
    {
        const long bi = bond.bit_i;
        const long bj = bond.bit_j;
        const long mask = bi | bj;
        const RealType half_J = bond.half_J;
        // 0.5 * ( S_i^+ * S_j^- + S_i^- * S_j^+ ):
        // iterate directly over the anti-aligned pairs, ident has (bit_i, bit_j) = (1, 0)
        // and its partner ident^mask has (0, 1); no branch needed
        for( long high = 0; high < dim; high += 2*bj )
        {
            for( long mid = high; mid < high + bj; mid += 2*bi )
            {
                for( long ident = mid + bi; ident < mid + 2*bi; ++ident )
                {
                    const long partner = ident ^ mask;
                    new_state[ident] += half_J * state[partner];
                    new_state[partner] += half_J * state[ident];
                }
            }
        }
    }
    return new_state;
}

void Hamiltonian::act( const State& in, State& prev_inout ) const
// Chebyshev recurrence step, fused and allocation-free:
// prev_inout = 2*H*in - prev_inout, i.e. |psi_n> = 2H|psi_n-1> - |psi_n-2> in place.
{
    for( long ident = 0; ident < dim; ++ident )
    {
        prev_inout[ident] = RealType{2.} * diag[ident] * in[ident] - prev_inout[ident];
    }
    for( const Bond& bond : bonds )
    {
        const long bi = bond.bit_i;
        const long bj = bond.bit_j;
        const long mask = bi | bj;
        const RealType J = RealType{2.} * bond.half_J; // recurrence factor 2 absorbed
        for( long high = 0; high < dim; high += 2*bj )
        {
            for( long mid = high; mid < high + bj; mid += 2*bi )
            {
                for( long ident = mid + bi; ident < mid + 2*bi; ++ident )
                {
                    const long partner = ident ^ mask;
                    prev_inout[ident] += J * in[partner];
                    prev_inout[partner] += J * in[ident];
                }
            }
        }
    }
}

void Hamiltonian::build_tables( const RealType shift )
// Precompute the bond list and the full diagonal from the current couplings.
// Called once with shift=0 (unrescaled H) for the bandwidth search and once
// with shift=b/a after the rescaling.
{
    bonds.clear();
    RealType constant = -shift;
    for( long i = 0; i < numSpins; ++i )
    {
        for( long j = i; j < numSpins; ++j )
        {
            RealType J = couplings(i,j);
            if( J == RealType{0.0} ) continue; // skip zero couplings
            if( i == j )
            {
                constant += RealType{0.25} * J * zz_factor; // self-coupling: always aligned
            }
            else
            {
                bonds.push_back( Bond{ 1L << i, 1L << j, RealType{0.5} * J } );
            }
        }
    }

    diag.resize( dim );
    diag = constant;
    for( const Bond& bond : bonds )
    {
        const long mask = bond.bit_i | bond.bit_j;
        const RealType quarter_J_zz = RealType{0.5} * bond.half_J * zz_factor;
        for( long ident = 0; ident < dim; ++ident )
        {
            // S_i^z * S_j^z
            const long overlap = ident & mask;
            diag[ident] += ( overlap == 0 || overlap == mask ) ? quarter_J_zz : -quarter_J_zz;
        }
    }
    if( params.h_z != RealType{0.} )
    {
        for( long ident = 0; ident < dim; ++ident )
        {
            // h_z * sum_i S_i^z = h_z/2 * (n_up - n_down)
            const long n_up = std::popcount( static_cast<unsigned long long>( ident ) );
            diag[ident] += RealType{0.5} * params.h_z * static_cast<RealType>( 2*n_up - numSpins );
        }
    }
}


// Vector initializer for Lanczos
template <typename T>
void vector_initializer(std::vector<T>& v);

template <>
void vector_initializer(std::vector<ComplexType>& v) {
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> rand(-1.0, 1.0);

  size_t n = v.size();
  for (size_t i = 0; i < n; ++i) {
    v[i] = ComplexType(rand(mt), rand(mt));
  }
}

Hamiltonian::Hamiltonian( ps::ParameterSpace& pspace ):
numSpins( pspace.num_Spins ),
dim( pspace.HilbertSpaceDimension ),
couplings( pspace.couplings )
{
    params.lambda = pspace.lambda;
    params.h_z = pspace.h_z;
    if( pspace.spin_model == "XXZ" )
    {
        zz_factor = params.lambda;
    }

    // Find the smallest and largest eigenvalue for the rescaling
    // Define matrix-vector multiplication in a suitable form
    auto mv_mul = [&](const std::vector<ComplexType>& in, std::vector<ComplexType>& out){
        State psi( dim );
        std::copy( in.begin(), in.end(), psi.begin() );
        psi = act( psi );
        out.assign( psi.begin(), psi.end() );
    };
    if( pspace.determine_bandwidth )
    {
    build_tables( RealType{0.} ); // unshifted, unrescaled H for the bandwidth search
    // Apply lanczos. The Chebyshev rescaling only needs the bandwidth to a few digits,
    // so a loose tolerance suffices; the result is padded below. This matters because
    // each Lanczos iteration re-orthogonalizes against (and stores) all previous
    // Lanczos vectors, so the cost and memory grow quickly with the iteration count.
    lambda_lanczos::LambdaLanczos<ComplexType> engine_max(mv_mul, dim, true, 1);
    engine_max.init_vector = vector_initializer<ComplexType>;
    engine_max.eps = RealType{1e-5};
    std::vector<ComplexType> eigenvector;
    engine_max.run(E_max, eigenvector);
    lambda_lanczos::LambdaLanczos<ComplexType> engine_min(mv_mul, dim, false, 1);
    engine_min.init_vector = vector_initializer<ComplexType>;
    engine_min.eps = RealType{1e-5};
    engine_min.run(E_min, eigenvector);
    // Pad the measured bandwidth: the Chebyshev expansion requires the spectrum to lie
    // strictly inside [E_min, E_max] (overestimating is safe, undershooting diverges)
    const RealType margin = RealType{0.02} * ( E_max - E_min );
    E_max += margin;
    E_min -= margin;
    pspace.E_max = E_max;
    pspace.E_min = E_min;
    }
    else
    {
        E_max = pspace.E_max;
        E_min = pspace.E_min;
    }
    a = (E_max - E_min) / RealType{2.};
    b = (E_max + E_min) / RealType{2.};
    b_over_a = b/a;
    couplings = couplings / a;
    params.h_z = params.h_z / a;
    build_tables( b_over_a ); // final tables for the rescaled, shifted H
}

};
