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

std::vector<ComplexType> CET_coefficients( const RealType t, const RealType a, const RealType b, const std::string& evol_type, const RealType tolerance );

ComplexType cdot( const State& state1, const State& state2 );

size_t generate_seed( const ps::ParameterSpace& pspace, const size_t my_rank );

size_t throw_seed( const size_t seed, const size_t my_rank, const size_t sample );

bool initialize_state( const ps::ParameterSpace& pspace, size_t seed, size_t sample, State& state );

void compute_correlations_at( int t, long site, const ps::ParameterSpace& pspace, State& psi_L, States& v_psi_R, CorrelationTensor& corrs_Re, CorrelationTensor& corrs_Im );

State S_alpha_i_act( const State& state, const long site, char alpha );

States S_i_act( const ps::ParameterSpace& pspace, State& state, const long site );

void CET( const ham::Hamiltonian& H, State& state, const std::vector<ComplexType>& coeffs );

void MPI_share_results( CorrelationTensor& correlations_Re, CorrelationTensor& correlations_Im );

void normalize( RealType& partition_function, CorrelationTensor& spin_c );

void add_sqs( const CorrelationTensor& correlations_Re, const CorrelationTensor& correlations_Im, CorrelationTensor& sqsum_Re, CorrelationTensor& sqsum_Im );

void compute_stds( RealType M, const CorrelationTensor& correlations_Re, const CorrelationTensor& correlations_Im, const CorrelationTensor& sqsum_Re, const CorrelationTensor& sqsum_Im, CorrelationTensor& stds_Re, CorrelationTensor& stds_Im );

}
