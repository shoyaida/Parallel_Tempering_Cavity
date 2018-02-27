# Parallel_Tempering_Cavity

Description
-----------
Perform parallel-tempering Monte Carlo (MC) simulation to efficiently sample configurations inside the cavity while keeping outside intact

Files
-----

- CavityPT_KA.c  
Main code: details of the algorithm are stipulated in L. Berthier, P. Charbonneau, and S. Yaida, "Efficient measurement of point-to-set correlations and overlap fluctuations in glass-forming liquids," J. Chem. Phys. 144, 024501 (2016).

- ran_uniform.c
Uniform random number generator code.
This is redistribution: Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura; GNU Library General Public License.

- ran_uniform.h
Header file for uniform random generator.

- makefile
Makefile that generates the executive file, CPT.

- test_SamplingParameters.dat
1st element=Ndump=dumping frequency of cavity snapshots;
2nd element=Nsnap=number of total snapshots;
3rd element=cavity_size=the size of cavity in which particles move.

- test_PTParameters.dat
1st element=n_replica=number of replicas;
the rest are columns consisting of temperatures and lambdas.

- test_CavitySample.dat
1st element=NA=number of A-species particles;
2nd element=NB=number of B-species particles;
3rd element=beta=inverse temperature;
4th element=IRcutoff=maximal linear distance from the core of the cavity for which particle configurations are recorded in test_CavitySample.dat;
the rest are columns consisting of (x,y,z,r,species)
