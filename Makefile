# Compiler and flags configuration
# NVidia HPC SDK OpenMPI 
FC = mpifort
FFLAGS = -O3 -Mpreprocess -D_MPI 

# Debugging flags (uncomment when needed)
# FFLAGS = -Mpreprocess -Mbounds

# Main program source file 
MAIN_SRC = SBP.f90

# Other source files 
SRC = random.f90 sorts.f90 normals.f90 newuoa.f90 global_data.f90 arrays.f90 timers.f90 functions.f90 marginals.f90 MPI_communications.f90 rootfindmod.f90

# Combined source files 
ALL_SRC = $(SRC) $(MAIN_SRC)

# Object files 
OBJ = $(ALL_SRC:.f90=.o)

# Module files 
MOD = $(SRC:%.f90=%.mod)

# Executable name 
EXEC = $(basename $(MAIN_SRC)).exe

# Default target
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJ)
	$(FC) $(OBJ) -o $(EXEC) $(FFLAGS)

# Rule to compile Fortran source files into object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJ) $(MOD) $(EXEC)

# Run the program (using 7 processes, change as needed)
run: $(EXEC)
	mpirun -np 7 ./$(EXEC)

# Phony targets (not actual files)
.PHONY: all clean run
