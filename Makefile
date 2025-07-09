cppOptions 	=	-std=c++17 \
				-pipe \
				-Wall \
				-O3 \

# additional flags:
#-ffast-math: quicker but "dangerous", needs to be tested if included

# self-written C++ libraries:
CPP_LIBS = ../cpp_libs

# === TODO ===
# diagonalization flags (use either LAPACK or EIGEN, remove or comment out the remaining flag):
diagonal	= 	-I path/to/eigen \
				-llapack \
				-D EIGEN \
#				-D LAPACK \

# === TODO ===
# hdf5 flags (lhdf5 or lhdf5_serial should work):
hdf5		=	-I /path/to/hdf5 \
				-lhdf5 \

# === TODO ===
# used data type (float or double)
USE_TYPE	= 	FLOAT
#DOUBLE

CPP_LIBS = ../cpp_libs

#DOUBLE
#FLOAT

libOptions 	=	-fopenmp \
				-lboost_program_options \
				-D USE_$(USE_TYPE) \
				$(diagonal) \
				$(hdf5) \
				-I $(CPP_LIBS) \
				-I /data/bieniek \

COMPILER	=	g++ # mpicxx
EXECUTABLE  =	executable_$(USE_TYPE).out

REAL_H		=	Global/RealType_Global.h
MPI_H		=	Global/MPI_Global.h
MATR_H		=	Global/Matrices_Global.h
INLMATR_H	=	Global/InlineMatrices_Global.h
EIGEN_H		=	Global/Eigen_Global.h

OUT		   	=	Executables_$(USE_TYPE)
TMM_OUT		=	$(OUT)/TimeMeasure.out
PSPACE_OUT	=	$(OUT)/Parameter_Space.out
RTD_OUT		=	$(OUT)/Run_Time_Data.out
MOD_OUT		= 	$(OUT)/Models.out
INIT_OUT	=	$(OUT)/Initialization.out
STORE_OUT	=	$(OUT)/Storage_Concept.out
MAIN_OUT	=	$(OUT)/main.out

# move errors to file : 2> output

all: $(EXECUTABLE)

$(EXECUTABLE): $(MAIN_OUT) | build 
	$(COMPILER) $(cppOptions) $(MAIN_OUT) -o $(EXECUTABLE) $(libOptions)

#
$(MAIN_OUT): main.cpp | build
	$(COMPILER) $(cppOptions) -c -o $(MAIN_OUT) main.cpp $(libOptions)


build:
	mkdir -p $(OUT)/

clean:
	rm -rf $(OUT); rm $(EXECUTABLE)

# Models/ISO_Disordered_Blockwise_NN.cpp Models/ISO_Disordered.cpp Models/ISO_Polarized.cpp Models/XXZ_Disordered.cpp Models/DRF_Disordered_Blockwise.cpp Models/DRF_Disordered.cpp