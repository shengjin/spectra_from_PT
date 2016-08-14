ledebut := $(shell echo 'spectra')
#labranche := $(shell basename ${PWD})
lemilieu:= $(shell echo '_')
#lafin := $(shell svn info | grep "Last Changed Rev" | cut -f 2 -d: | sed -e 's/^[ \t]*//')
lafin :=$(shell echo 'X')
codename := ${ledebut}${labranche}${lemilieu}${lafin}

hostname := $(shell hostname )
#hostname := $(shell hostname -s)
#standart compiler
compiler := gfortran
warning := false
#warning := true

DF=
DF90=
DL=

ifeq "$(compiler)" "gfortran"
   FC=gfortran
   F90C=gfortran
   LD=gfortran
   FF=-ffixed-line-length-132 -O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow -Wno-align-commons -c -w -C -fdefault-real-8 -fdefault-double-8 
   F90F=-ffree-line-length-none -O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow -Wall -c -w -fdefault-real-8 -fdefault-double-8     
   LDF=-O3 -Wall -lstdc++ -finit-real=nan
   #additional new optimsation settings: -flto -fstack-arrays
endif


OBJ = readstructure.o getline.o spectra.o readopacity.o natconst.o PTSpetc.o interpolopacity.o calcspectrum.o

all: $(OBJ) $(LIB)
	echo 'Building' ${codename}
	echo 'hostname' ${hostname}
	${LD}  ${LDF} ${DL}-o ${codename} $(MOD) $(OBJ) $(LIB)

%.o: %.f
	${FC} ${FF} ${DF}$<

%.o: %.f90 $(MOD)
	${F90C} ${F90F} ${DF90}$<


clean: clean_program 
#       -rm *.o
#       -rm *.mod
#       cd $(DIREMPS); $(MAKE) clean

clean_program:
	-rm *.o
	-rm *~
