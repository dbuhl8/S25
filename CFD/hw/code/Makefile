OBJECTS = LinAl.o NumDE.o CFD.o advect.o plm.o fog.o write_data.o cfl.o bc.o\
	advect_init.o  grid_init.o advect_update.o upwind.o downwind.o centered.o \
	minmod.o mc.o vanLeers.o mm_2arg.o mm_3arg.o vl.o
	
MODULES = LinAl.mod NumDE.mod CFD.mod

OUTPUTS = *.png *.mod output.txt

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g -llapack -lblas -fopenmp

.PHONY: clean


main: advect.e
	./advect.e

advect.e: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o advect.e

advect.o: CFD.o NumDE.o LinAl.o advect_update.o

advect_update.o: CFD.o NumDE.o LinAl.o

CFD.o: NumDE.o plm.o fog.o write_data.o cfl.o bc.o advect_init.o \
	grid_init.o upwind.o downwind.o centered.o \
	minmod.o mc.o vanLeers.o

plm.o: NumDE.o upwind.o downwind.o centered.o minmod.o mc.o vanLeers.o
fog.o: NumDE.o

NumDE.o: LinAl.o mm_2arg.o mm_3arg.o vl.o

upwind.o: NumDE.o
downwind.o: NumDE.o
centered.o: NumDE.o
minmod.o: NumDE.o 
mc.o: NumDE.o 
vanLeers.o: NumDE.o 

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(OUTPUTS) $(MODULES) advect.e
