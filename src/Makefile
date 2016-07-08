FC=gfortran
CFLAGS = -Wall -O2 -fbounds-check

#***name of program and objects
PROGRAM = sinzi.exe
OBJS    = types.o \
	  linear_algebra.o\
          inputs_class.o \
          pdb_class.o\
	  calc.o     \
	  align_class.o\
          presto_class.o \
	  analysis_class.o\
          main.o 

#***Making program
$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(CFLAGS) 

types.o: types.f03
	$(FC) -c $<

linear_algebra.o: linear_algebra.f03
	$(FC) -c $<

inputs_class.o: inputs_class.f03
	$(FC) -c $<

pdb_class.o: pdb_class.f03
	$(FC) -c $<

calc.o: calc.f03
	$(FC) -c $<

align_class.o: align_class.f03
	$(FC) -c $<

analysis_class.o: analysis_class.f03
	$(FC) -c $<

presto_class.o: presto_class.f03 
	$(FC) -c $<

main.o: main.f03 
	$(FC) -c $<
#
.PHONY:clean
clean:
	rm $(OBJS) $(PROGRAM) *.mod

