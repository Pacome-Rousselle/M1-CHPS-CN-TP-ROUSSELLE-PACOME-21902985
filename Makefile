##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
# -- option to dedicated machine
#
# Rk: You should create a file such as your-machineName.mk
# Follow the example of the machine called "ambre" in the 
# file ambre.mk
#
HOSTNAME=$(shell hostname)
# todo later : make an option to rename current .mk file 
# to $(HOSTNAME).mk if not defined
include $(HOSTNAME).mk

# 
# -- Compiler Option
OPTC=${OPTCLOCAL}

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=${LIBSLOCAL}

# -- Include directories
INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)

#
#################################################################
# makefile
############
#
OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp_poisson1D_direct.o
#

all: testenv tpPoisson1D_iter tpPoisson1D_direct
	rm -f *.o
	./plotting.sh

testenv: bin/tp_testenv
	$^

tpPoisson1D_iter: bin/tpPoisson1D_iter
	$^

tpPoisson1D_direct: bin/tpPoisson1D_direct
	$^

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp_poisson1D_iter.o: $(TPDIRSRC)/tp_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_poisson1D_iter.c  

tp_poisson1D_direct.o: $(TPDIRSRC)/tp_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_poisson1D_direct.c  

bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tpPoisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tpPoisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tpPoisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tpPoisson1D_iter:
	bin/tpPoisson1D_iter

run_tpPoisson1D_direct:
	bin/tpPoisson1D_direct

clean:
	rm -f *.o *.dat 
	rm -f bin/*
	rm -f dat/direct/analytical/* dat/direct/experimental/*
	rm -f dat/iter/analytical/*.dat dat/iter/experimental/*.dat
	rm -f dat/iter/experimental/jacobi/* dat/iter/experimental/gauss-seidel/*
	rm -f dat/direct/*.dat dat/iter/*.dat
	rm -f graphs/direct/*.dat graphs/iter/*.dat
	rm -f graphs/direct/*.png graphs/iter/*.png
