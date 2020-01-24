.SUFIXES: .f90

# Compiler
FC = gfortran

# Flags
# Flags for debugging 
#FFLAGS = -Wall -O0
# Flags for production
FFLAGS = -O3
P1=porE
P2=porE_subgrid
P4=get_PSD

# Default target

all: porE porE_subgrid get_PSD

porE:  $(P1).f90
	$(FC) $(P1).f90 -o $(P1) $(FFLAGS)

porE_subgrid: $(P2).f90
	$(FC) $(P2).f90 -o $(P2) $(FFLAGS)

get_PSD: $(P4).f90
	$(FC) $(P4).f90 -o get_PSD $(FFLAGS) 

clean:
	rm -f porE porE_subgrid get_PSD
