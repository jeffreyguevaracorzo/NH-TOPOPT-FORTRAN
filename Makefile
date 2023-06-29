# Compilation parameters
CC = gfortran
CFLAGS = -O1 -fbounds-check -fbacktrace -fcheck=all -g -Wall -Wextra -Wrealloc-lhs-all
TARGET = ExecuteOptimization
EXT = -llapack -lblas -L/opt/homebrew/Cellar/metis/5.1.0/lib -lmetis 
FILES = FortranModules/HSLModules/sdeps90.f90 \
		FortranModules/HSLModules/ddeps90.f90 \
		FortranModules/HSLModules/common90.f90 \
		FortranModules/HSLModules/hsl_ma87s.f90 \
		FortranModules/HSLModules/hsl_ma87d.f90 \
		FortranModules/HomogenizationModule.f90 \
		FortranModules/TopOptimizationModule.f90 \
		FortranModules/PostProcessingModule.f90 \
		Main.f90
OTHERFILES = DataStructure/AdditionalData/*.txt \
			 ParaviewPostprocessing/*.case \
			 ParaviewPostprocessing/*.geom \
			 ParaviewPostprocessing/*.esca \
			 NumericalResults/*.txt \
			 DataStructure/AdditionalData/SparseSystem/*.txt

# Object files
OBJ1 = ${FILES:.f90=.o}

# Compilation and cleanup command
%.o : %.f90
	${CC} ${CFLAGS} -o $@ -c $<

${TARGET} : ${OBJ1}
	${CC} ${CFLAGS} -o $@ ${OBJ1} $(EXT)

.PHONY : clean
clean :
	@rm -f *.o *.mod ${TARGET} ${OBJ1} $(OTHERFILES)