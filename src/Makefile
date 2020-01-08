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
P3=porE_window
P4=get_PSD

# Default target

all: porE porE_subgrid porE_window get_PSD

porE:  $(P1).f90
	$(FC) $(P1).f90 -o $(P1) $(FFLAGS)

porE_subgrid: $(P2).f90
	$(FC) $(P2).f90 -o $(P2) $(FFLAGS)

porE_window: $(P3).f90
	$(FC) $(P3).f90 -o $(P3) $(FFLAGS)

get_PSD: $(P4).f90
	$(FC) $(P4).f90 -o get_PSD $(FFLAGS) 

clean:
	rm -f porE porE_subgrid pore_window get_PSD
