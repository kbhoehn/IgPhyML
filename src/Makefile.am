# Makefile.am for IgPhyML
# Modified by Ken Hoehn 12/1/2017

# Mac compilation options. Uncomment CC line and then replace values 
# after '=' CC requires the llvm clang executable, MACOMP requires 
# the location of omp.h from llvm, and MACLLVM requires the llvm
# library directory. Examples are given below - if you've installed
# llvm through homebrew, yours will probably look similar.

#CC=/usr/local/Cellar/llvm/3.9.1/bin/clang
#MACOMP = /usr/local/Cellar/llvm/3.9.1/lib/clang/3.9.1/include/omp.h
#MACLLVM = /usr/local/Cellar/llvm/3.9.1/lib


DEFS=-DUNIX -D$(PROG) -DDEBUG -D$(OPT)

bin_PROGRAMS = igphyml
PROG = PHYML


igphyml_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
controller.c controller.h\
cl.c  cl.h\
spr.c spr.h\
stats.c stats.h\
nucle2codon.c nucle2codon.h

igphyml_LDADD =
igphyml_LDFLAGS =  


if MACOS
igphyml_LDFLAGS += -L $(MACLLVM)
endif

if WANT_BLAS
OPT = BLAS
if MACOS
igphyml_LDFLAGS += -framework Accelerate
else
igphyml_LDADD += -llapack_atlas -llapack -lblas -latlas
endif	#MACOS
else
if WANT_OMP
OPT = OMP
if MACOS
igphyml_LDFLAGS += -I $(MACOMP)
else
igphyml_LDADD += -lgomp
endif	#MACOS
else
if WANT_BLAS_OMP
OPT = BLAS_OMP
if MACOS
igphyml_LDFLAGS += -framework Accelerate -I $(MACOMP)
else
igphyml_LDADD += -lgomp -llapack_atlas -llapack -lblas -latlas
endif	#MACOS
else
OPT=NBLAS
endif	#WANT_BLAS_OMP
endif 	#WANT_OMP
endif	#WANT_BLAS

#added by Ken 16/5/2016
#Includes directory for motifs tables.
.c.o:
	$(COMPILE) -D_GNU_SOURCE -DHTABLE='"$(CURDIR)/motifs/HTABLE_"' -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
	
#added by Ken 16/5/2016
igphyml_LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(igphyml_LDFLAGS) \
	$(LDFLAGS) -o igphyml
