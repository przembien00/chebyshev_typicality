#pragma once

#include<complex>
#include<blaze/Math.h>
#include<mpi.h>


// REAL AND COMPLEX NUMBERS

#ifdef USE_DOUBLE
typedef double RealType;
inline MPI_Datatype MPI_REALTYPE = MPI_DOUBLE;
#endif 

#ifdef USE_FLOAT 
typedef float RealType;
inline MPI_Datatype MPI_REALTYPE = MPI_FLOAT;
#endif



typedef std::complex<RealType> ComplexType;     

// STATE

typedef blaze::DynamicVector<ComplexType> State;
typedef std::vector<State> States;

// SYSTEM COUPLINGS

typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<RealType,blaze::rowMajor>> CouplingMatrix;

// CORRELATION

typedef std::vector<ComplexType> Correlation;