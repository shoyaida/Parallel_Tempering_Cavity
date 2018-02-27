CC = gcc
CFLAGS  = -O3 -Wall -lm

default: CavityPT

CavityPT:  CavityPT_KA.o ran_uniform.o 
	$(CC) $(CFLAGS) -o CPT CavityPT_KA.o ran_uniform.o

CavityPT_KA.o:  CavityPT_KA.c ran_uniform.h
	$(CC) $(CFLAGS) -c CavityPT_KA.c

ran_uniform.o: ran_uniform.c ran_uniform.h 
	$(CC) $(CFLAGS) -c ran_uniform.c
clean: 
	$(RM) count *.o *~
