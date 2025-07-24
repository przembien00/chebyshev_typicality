#include"Hamiltonians.h"
#include<functional>
#include<cmath>
#include<string>
#include<map>
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"
#include<iostream>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <lambda_lanczos/lambda_lanczos_tridiagonal_impl.hpp>

namespace ps = Parameter_Space;

namespace Hamiltonians
{

State act_ISO( const State& state, const Hamiltonian& H )
// Apply the Isotropic Heisenberg Hamiltonian to the state vector.
{   
    State new_state = - H.b_over_a * state; // shift to make the spectrum symmetrical

    for( long i = 0; i < H.numSpins; ++i )
    {
        for( long j = i; j < H.numSpins; ++j )
        {
            RealType J = H.couplings(i,j);
            if( J == RealType{0.0} ) continue; // skip zero couplings

            for( long ident = 0; ident < H.dim; ++ident )
            {
                // Apply the Hamiltonian term to the state

                // S_i^z * S_j^z

                if( ( ( ident >> i ) & 1L) == ( (ident >> j)  & 1L ) ) // check if the spins on sites i, j are in the same direction
                {
                    new_state[ident] += RealType{0.25} * J * state[ident];
                }
                else
                {
                    new_state[ident] += - RealType{0.25} * J * state[ident];
                }

                // 0.5 * ( S_i^+ * S_j^- + S_i^- * S_j^+ )

                if( ( ( ident >> i ) & 1L) ^ ( (ident >> j)  & 1L ) )
                {
                    long new_ident = ( ident ^ ( 1L << j ) ) ^ ( 1L << i );
                    new_state[new_ident] += RealType{0.5} * J * state[ident];
                }

            }
        }
    }

    return new_state;
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

Hamiltonian::Hamiltonian( const ps::ParameterSpace& pspace ):
numSpins( pspace.num_Spins ),
CET_rescale( pspace.CET_rescale ),
dim( pspace.HilbertSpaceDimension ),
couplings( pspace.couplings )
{
    
    // Define a mapping from models to actions of the Hamiltonian
    std::map< std::string, std::function< State( const State&, const Hamiltonian& ) > > model_map{
        {"ISO", act_ISO},
    };
    // Define the action of the Hamiltonian based on the model.
    act = [model_map, model = pspace.spin_model, this](const State& state){
        return model_map.at(model)(state, *this);};

    // Find the smallest and largest eigenvalue for the rescaling
    // Define matrix-vector multiplication in a suitable form
    auto mv_mul = [&](const std::vector<ComplexType>& in, std::vector<ComplexType>& out){
        State psi(dim);
        for( int i=0; i<in.size(); i++)
        {
            psi[i] = in[i];
        }
        psi = act( psi );
        for( int i=0; i<psi.size(); i++)
        {
            out[i] = psi[i];
        }
    };
    // Apply lanczos
    lambda_lanczos::LambdaLanczos<ComplexType> engine_max(mv_mul, dim, true, 1);
    engine_max.init_vector = vector_initializer<ComplexType>;
    RealType eigenvalue_max;
    std::vector<ComplexType> eigenvector_max;
    engine_max.run(eigenvalue_max, eigenvector_max);
    lambda_lanczos::LambdaLanczos<ComplexType> engine_min(mv_mul, dim, false, 1);
    engine_min.init_vector = vector_initializer<ComplexType>;
    RealType eigenvalue_min;
    std::vector<ComplexType> eigenvector_min;
    engine_min.run(eigenvalue_min, eigenvector_min);
    a = (eigenvalue_max - eigenvalue_min) / RealType{2.};
    b = (eigenvalue_max + eigenvalue_min) / RealType{2.};
    // b = 0;
    b_over_a = b/a;
    // couplings = couplings / a;
    // Initialize Hamiltonian based on couplings
    // Find the largest abs.val. coupling J, divide by 0.5*J*z*N
    // to make the spectrum of H fit into [-1,1] 
    // RealType J = RealType{0.0};
    // for( int i = 0; i < numSpins; ++i )
    // {
    //     for( int j = 0; j < numSpins; ++j )
    //     {
    //         if( std::abs(couplings(i,j)) > J )
    //         {
    //             J = std::abs(couplings(i,j));
    //         }
    //     }
    // }
    // a = RealType{0.25} * J * static_cast<RealType>(numSpins) * CET_rescale;
    // std::cout << "Old bound = " << a << '\n';
    couplings = couplings / a;
}

};