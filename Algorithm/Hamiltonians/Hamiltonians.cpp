#include"Hamiltonians.h"
#include<functional>
#include<string>
#include<map>
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"
#include<iostream>


namespace ps = Parameter_Space;

namespace Hamiltonians
{

State act_ISO( const State& state, const Hamiltonian& H )
// Apply the Isotropic Heisenberg Hamiltonian to the state vector.
{   
    State new_state(H.dim);

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

Hamiltonian::Hamiltonian( const ps::ParameterSpace& pspace ):
numSpins( pspace.num_Spins ),
CET_rescale( pspace.CET_rescale ),
dim( pspace.HilbertSpaceDimension ),
couplings( pspace.couplings )
{
    // Initialize Hamiltonian based on couplings
    // Find the largest abs.val. coupling J, divide by 0.5*J*z*N
    // to make the spectrum of H fit into [-1,1] 
    RealType J = RealType{0.0};
    for( int i = 0; i < numSpins; ++i )
    {
        for( int j = 0; j < numSpins; ++j )
        {
            if( std::abs(couplings(i,j)) > J )
            {
                J = std::abs(couplings(i,j));
            }
        }
    }
    a = RealType{0.25} * J * static_cast<RealType>(numSpins) * CET_rescale;
    // couplings = couplings / a;
    
    // Define a mapping from models to actions of the Hamiltonian
    std::map< std::string, std::function< State( const State&, const Hamiltonian& ) > > model_map{
        {"ISO", act_ISO},
    };
    // Define the action of the Hamiltonian based on the model.
    act = [model_map, model = pspace.spin_model, this](const State& state){
        return model_map.at(model)(state, *this);};
}

};