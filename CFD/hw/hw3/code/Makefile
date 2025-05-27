OBJECTS = LinAl.o NumDE.o CFD.o driver.o # Make sure the correct driver is in the objects variable
	
MODULES = LinAl.mod NumDE.mod CFD.mod

OUTPUTS = *.png *.mod output.txt 

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -g -llapack -lblas -fopenmp

.PHONY: clean


main: driver.e
	./driver.e

final: driver.e
	./driver.e
	gnuplot final_plot.p

hw4: driver.e
	./driver.e > output.txt
	gnuplot hw4_plot.p

hw3: driver.e
	./driver.e > output.txt
	gnuplot hw3_plot.p

midterm: driver.e
	./driver.e > output.txt
	gnuplot midterm_plot.p

hw2: driver.e
	./driver.e > output.txt
	gnuplot hw2_plot.p

hw1: driver.e
	./driver.e > output.txt
	gnuplot hw1_plot.p

driver.e: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o driver.e

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(OUTPUTS) $(MODULES) driver.e
