////////Perform parallel-tempering Monte Carlo (MC) simulation to efficiently sample configurations inside the cavity while keeping outside intact////////
////////Used to measure point-to-set correlations in glass-forming liquids////////
////////PARALLEL TEMPERING: number of replica must be greater than or equal to 2////////
////////Cavity is set to be an ``open ball"////////
////////Here illustrated for the Kob-Andersen binary Lennard-Jones liquid////////
////////For general potentials, one needs to change parts of the code inside the parameters inside Cavity_Equilibration_MC////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran_uniform.h" /*uniform random number generating code*/
#include <time.h>

////////Parameters for the cell list////////
#define NCELLMAX 50 /*max number of particles in one cell*/
#define NEIGHBORMAX 120 /*max number of neighbor particles per particle*/
#define INTERACTIONRANGE 2.5 /*chosen such that it is greater than or equal to the maximal interaction range of cutoff Lennard-Jones potential*/

////////Parameters for Monte Carlo////////
#define MAXMOVE 0.3 /*max_move_length*/
#define ADJUSTENERGY 100 /*Re-calculate energy profile every ADJUSTENERGY macro MC sweeps*/
#define TBALLISTIC 1000.0 /*Parallel-tempering swaps of replicas are attempted once per TBALLISTIC macro-MC sweeps*/

#define DIM 3 /*number of spatial dimensions*/
#define DIMPOWER(x) ((x)*(x)*(x)) /*x^DIM*/
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define EMPTY 0
#define LENGTH 255
#define TRUE 1
#define FALSE 0
void Cavity_Equilibration_MC(long long int macroMC_step,double lambda,double hard_wall_size,int N_mobile_tot,int NA_mobile,int NA_pinned,double max_move_length,int n,double hL, double l_cell, double beta,double etot_list[],double x_M[],double y_M[],double z_M[],int a_M[],double x_P[],double y_P[],double z_P[],int a_P[],int cell_has_which_Mparticle[],int cell_has_which_Pparticle[],int in_which_cell_M[],int in_which_cell_P[],int where_in_cell_list_M[],int where_in_cell_list_P[],int cell_has_how_much_Mparticle[],int cell_has_how_much_Pparticle[],int neighbor_search_facilitator_x[],int neighbor_search_facilitator_y[],int neighbor_search_facilitator_z[],double interaction_range);
void Energy_Profiler(double lambda,int N_mobile_tot,int N_pinned_tot,int NA_mobile,int NA_pinned,int n,double etot_list[],double x_M[],double y_M[],double z_M[],int a_M[],double x_P[],double y_P[],double z_P[],int a_P[],int cell_has_which_Mparticle[],int cell_has_which_Pparticle[],int in_which_cell_M[],int in_which_cell_P[],int cell_has_how_much_Mparticle[],int cell_has_how_much_Pparticle[],int neighbor_search_facilitator_x[],int neighbor_search_facilitator_y[],int neighbor_search_facilitator_z[],double interaction_range);
void Take_One_Replica_Out(int replica_index,int N_mobile_tot,int N_pinned_tot,int n_replica,int c_index_max,double x_M[],double y_M[],double z_M[],int in_which_cell_M[],int cell_has_which_Mparticle[],int where_in_cell_list_M[],int cell_has_how_much_Mparticle[], double etot_list[],double x_M_replicated[],double y_M_replicated[],double z_M_replicated[],int in_which_cell_M_replicated[],int cell_has_which_Mparticle_replicated[],int where_in_cell_list_M_replicated[],int cell_has_how_much_Mparticle_replicated[], double etot_list_replicated[]);
void Put_One_Replica_Back(int replica_index,int N_mobile_tot,int N_pinned_tot,int n_replica,int c_index_max,double x_M[],double y_M[],double z_M[],int in_which_cell_M[],int cell_has_which_Mparticle[],int where_in_cell_list_M[],int cell_has_how_much_Mparticle[], double etot_list[],double x_M_replicated[],double y_M_replicated[],double z_M_replicated[],int in_which_cell_M_replicated[],int cell_has_which_Mparticle_replicated[],int where_in_cell_list_M_replicated[],int cell_has_how_much_Mparticle_replicated[], double etot_list_replicated[]);	
int main(){
	
////////////////Initial assignment of variables////////////////
    
	long long int i; /*generic integer dummy variable*/
    char filename[LENGTH];
	FILE *infile;
	FILE *outfile;
    
    /*Initialize uniform random number generator*/
    double seed=(double)222.0; /*a seed number for uniform random number generator*/
    InitializeRandomNumberGenerator(seed);
    
    /*Each particle information*/
    double x_temp; /*temporal x-coordinate*/
    double y_temp; /*temporal y-coordinate*/
    double z_temp; /*temporal z-coordinate*/
    double radial_distance; /*radial distance from the core of the cavity*/
    int a_species;   /*particle species label: 0 for species A, 1 for species B*/

	/*Some counting of particles*/
	int NA;          /*number of A-species particles*/
	int NB;          /*number of B-species particles*/
    int N;           /*number of particles: N=NA+NB*/
    int NA_mobile;   /*number of A-species mobile particles inside the cavity*/
    int NB_mobile;   /*number of B-species mobile particles inside the cavity*/
    int N_mobile_tot;/*number of total mobile particles inside the cavity*/
    int NA_pinned;   /*number of A-species pinned particles outside the cavity, within interaction range*/
    int NB_pinned;   /*number of B-species pinned particles outside the cavity, within interaction range*/
    int N_pinned_tot;/*number of total pinned particles outside the cavity, within interaction range*/
    
    /*Cell list parameters*/
    int n; /*n counts the number of cells in one direction*/
    int n_cell_max=NCELLMAX; /*n_cell_max sets the maximum number of particles in one cell allowed*/
    double interaction_range=(double)INTERACTIONRANGE; /*maximum range of interaction; this code automatically EXITs if not chosen big enough*/
    double l_cell; /*linear size of the cell*/
    double IRcutoff; /*maximal linear distance from the core of the cavity for which particle configurations are recorded in test_CavitySample.dat*/
    double hL; /*half linear extent*/
	double cavity_size; /*size of the cavity*/
	double hard_wall_size; /*the place of the hard wall beyond which mobile particles cannot escape; hard_wall_size=cavity_size in this code*/

    int c_x; /*cell index in x direction*/
    int c_y; /*cell index in y direction*/
    int c_z; /*cell index in z direction*/
    int c_index; /*cell index=c_x+c_y*n+c_z*(n*n);*/
    int c_index_max; /*number of cells*/
    int index_within_cell; /*indexes particles within a given cell*/
    int n_particles_in_cell; /*counts the number of particles within a given cell*/
    
    int neighbor_cell_index;
    int neighbor_cell_index_max=DIMPOWER(3);
    int x_neighbor_search;
    int y_neighbor_search;
    int z_neighbor_search;
    
    /*Monte Carlo (MC) bookkeeping*/
    long long int Ndump; /*number of macro MC sweeps between dumping of snapshots*/
    long long int Nsnap; /*number of total snapshots to be taken over the course of simulation */
    long long int NEquiPro;/*NEquiPro=Ndump*Nsnap=number of total macro MC sweeps of the simulation*/
    long long int Nadjust=ADJUSTENERGY; /*re-calculate energy profile every Nadjust macro MC sweeps*/
    long long int counter; /*counts the number of macro MC sweeps*/
    long long int snap_index; /*indexes the snapshots*/
    long long int macroMC_step;/*macro MC step*/
    double energy_checker; /*check energy is not too far off (more than 10^(-6)) for every Nadjust macro MC sweeps*/
    
    /*Monte Carlo (MC) parallel-tempering (PT) parameters*/
    double max_move_length=(double)MAXMOVE;
    double temperature; /*temperature of the system*/
    double beta; /*inverse temperature*/
    double lambda; /*particle-shrinking factor; see J. Chem. Phys. 144, 024501 (2016) for more details*/
    double beta_interest; /*inverse temperature of interest whereat configurations need to be sampled*/
    
    int n_replica; /*number of replicas to be used*/
    int replica_index; /*index systems*/
    int thermal_index; /*index thermal properties*/
    
    double t_ballistic=TBALLISTIC;/*sets the average time one replica spends in the same thermal enviroment between swap trials*/
    double p_swap; /*(n_replica-1)/t_ballistic, setting the frequency of replcia swap attempts to be as stipulated above*/
    
    /*Temporal book-keeping variables used in parallel-tempering swap decisions*/
    int thermal_index_low;
    int thermal_index_high;
    int replica_index_low;
    int replica_index_high;
    double beta_low;
    double beta_high;
    double lambda_low;
    double lambda_high;
    double Etot_low_big;
    double Etot_high_small;
    double Etot_low_small;
    double Etot_high_big;
    double swap_decider;

    /*Variou FLAGs*/
    int SWAP_TRIAL_FLAG;
    int SWAP_FLAG;
    int ENERGY_ADJUSTMENT_FLAG;
    int SNAPSHOT_FLAG;
    int SAVE_FLAG;
    int count_down_till_save;
    
////////////////Read in cavity parameters to go over in this simulation sequence////////////////
    
    /*Set the sampling parameters*/
    sprintf(filename,"./test_SamplingParameters.dat");
    if((infile = fopen(filename, "r")) == NULL){printf("Error opening file %s\n",filename);exit(EXIT_FAILURE);}
    fscanf(infile,"%lld %lld %lf",&Ndump,&Nsnap,&cavity_size);
    NEquiPro=Nsnap*Ndump;
    hard_wall_size=cavity_size;
    fclose(infile);
    
    /*Set parallel-tempering parameters*/
    sprintf(filename,"./test_PTParameters.dat");
    if((infile = fopen(filename, "r")) == NULL){printf("Error opening file %s\n",filename);exit(EXIT_FAILURE);}
    fscanf(infile,"%d",&n_replica);
    /*Allocate memory for the string of replica parameters*/
    double *beta_string;
    double *lambda_string;
    beta_string = (double *)calloc(n_replica, sizeof(double));
    lambda_string = (double *)calloc(n_replica, sizeof(double));
    if (beta_string==NULL || lambda_string==NULL){printf("calloc of replica parameters failed. \n");exit(EXIT_FAILURE);}
    /*Assign string of replica parameters*/
    for (thermal_index=0;thermal_index<n_replica;thermal_index++){
        fscanf(infile, "%lf %lf", &temperature, &lambda);
        beta=(1.0/(temperature));
        beta_string[thermal_index]=beta;
        lambda_string[thermal_index]=lambda;
    }
    fclose(infile);
    
    p_swap=((double)n_replica-1)/t_ballistic; /*When p_swap is greater than 1, swap trial would occur less often than one per t_ballistic MC sweeps*/
    
////////////////Read in an original cavity(+its surrounding) configuration////////////////
	
	sprintf(filename,"./test_CavitySample.dat");
	if((infile = fopen(filename, "r")) == NULL){printf("Error opening file %s\n",filename);exit(EXIT_FAILURE);}
	
	/*Parameters of binary LJ liquid*/
	fscanf(infile,"%d %d %lf %lf",&NA,&NB,&beta_interest,&IRcutoff);
    N=NA+NB;
	
	/*Allocate memory for coordinate and radius*/
	double *x_ori;
	double *y_ori;
	double *z_ori;
	double *r_ori;
    int *a_ori;
	x_ori = (double *)calloc(NA+NB+1, sizeof(double));
	y_ori = (double *)calloc(NA+NB+1, sizeof(double));
	z_ori = (double *)calloc(NA+NB+1, sizeof(double));
	r_ori = (double *)calloc(NA+NB+1, sizeof(double));
    a_ori = (int *)calloc(NA+NB+1, sizeof(int));
	if (x_ori==NULL || y_ori==NULL || z_ori==NULL || r_ori==NULL || a_ori==NULL){printf("calloc of coordinates failed. \n");exit(EXIT_FAILURE);}

	/*Read in coordinate and radius*/
	for(i=1;i<NA+NB+1;i++){
		fscanf(infile,"%lf %lf %lf %lf %d",&x_temp,&y_temp,&z_temp,&radial_distance,&a_species);
		x_ori[i]=x_temp;
		y_ori[i]=y_temp;
		z_ori[i]=z_temp;
		r_ori[i]=radial_distance;
        a_ori[i]=a_species;
	}
	fclose(infile);
    if (hard_wall_size+interaction_range>IRcutoff){printf("We are missing relevant particles!\n");exit(EXIT_FAILURE);}
	
////////////////Identify unpinned particles////////////////
    
	NA_mobile=0;
    NB_mobile=0;
	for (i=1;i<(NA+NB+1);i++){
		if (r_ori[i]<cavity_size){
            a_species=a_ori[i];
            if(a_species==0){NA_mobile=NA_mobile+1;}
            if(a_species==1){NB_mobile=NB_mobile+1;}
		}
	}
	N_mobile_tot=NA_mobile+NB_mobile;
	if (N_mobile_tot==0){printf("No particle, no point. \n");exit(EXIT_FAILURE);}
	
	NA_pinned=NA-NA_mobile;
	NB_pinned=NB-NB_mobile;
	N_pinned_tot=NA_pinned+NB_pinned;
	
////////////////Separate interior and exterior information////////////////
	
	/*Allocate memory for coordinate and radius*/
	double *x_M;
	double *y_M;
	double *z_M;
    int *a_M;
	double *x_P;
	double *y_P;
	double *z_P;
    int *a_P;
	x_M = (double *)calloc(N_mobile_tot+1, sizeof(double));
	y_M = (double *)calloc(N_mobile_tot+1, sizeof(double));
	z_M = (double *)calloc(N_mobile_tot+1, sizeof(double));
    a_M = (int *)calloc(N_mobile_tot+1, sizeof(int));
	x_P = (double *)calloc(N_pinned_tot+1, sizeof(double));
	y_P = (double *)calloc(N_pinned_tot+1, sizeof(double));
    z_P = (double *)calloc(N_pinned_tot+1, sizeof(double));
    a_P = (int *)calloc(N_pinned_tot+1, sizeof(int));
	if (x_M==NULL || y_M==NULL || z_M==NULL || a_M==NULL || x_P==NULL || y_P==NULL || z_P==NULL || a_P==NULL){printf("calloc of coordinates failed. \n");exit(EXIT_FAILURE);}

	/*Read in coordinate*/
    N_mobile_tot=0;
    N_pinned_tot=0;
	for (i=1;i<(NA+NB+1);i++){
        if(r_ori[i]<cavity_size){
            x_M[i-N_pinned_tot]=x_ori[i];
            y_M[i-N_pinned_tot]=y_ori[i];
            z_M[i-N_pinned_tot]=z_ori[i];
            a_M[i-N_pinned_tot]=a_ori[i];
            N_mobile_tot++;
        }
        else{
            x_P[i-N_mobile_tot]=x_ori[i];
            y_P[i-N_mobile_tot]=y_ori[i];
            z_P[i-N_mobile_tot]=z_ori[i];
            a_P[i-N_mobile_tot]=a_ori[i];
            N_pinned_tot++;
        }
	}
	if (N_mobile_tot!=NA_mobile+NB_mobile || N_pinned_tot!=NA_pinned+NB_pinned){printf("Particle numbers screwed! \n");exit(EXIT_FAILURE);}
	free(x_ori);
	free(y_ori);
	free(z_ori);
	free(r_ori);
    free(a_ori);
	
////////////////Allocate memory for a cell structure and an energy profile ////////////////
	
	/*set cell structure parameters*/
	l_cell=interaction_range; /*linear size of one cell*/
	n=2*(ceil((IRcutoff/l_cell)-0.5))+1; /*n counts the number of cells in one direction*/
	hL=0.5*n*l_cell; /*half linear extent of the entire box spanned by cells*/
	c_index_max=DIMPOWER(n);
	if(max_move_length>l_cell){printf("MC move too large!\n");exit(EXIT_FAILURE);}
	
	/*Allocate memory for cell structure: mobile*/
	int *cell_has_which_Mparticle;/*list particles within a cell; 0 means empty*/
	int *in_which_cell_M; /*tells where a given particle resides*/
	int *cell_has_how_much_Mparticle; /*self-explanatory*/
	int *where_in_cell_list_M; /*facilitates cell reconstruction*/
	cell_has_which_Mparticle = (int *)calloc(c_index_max*n_cell_max, sizeof(int));
	in_which_cell_M = (int *)calloc(N_mobile_tot+1, sizeof(int));
	cell_has_how_much_Mparticle = (int *)calloc(c_index_max, sizeof(int));
	where_in_cell_list_M = (int *)calloc(N_mobile_tot+1, sizeof(int));
	if (cell_has_which_Mparticle==NULL || in_which_cell_M==NULL || cell_has_how_much_Mparticle==NULL || where_in_cell_list_M==NULL){printf("calloc of cell structure failed. \n");exit(EXIT_FAILURE);}

	/*Allocate memory for cell structure: pinned*/
	int *cell_has_which_Pparticle; /* 0 means empty*//*list particles within a cell*/
	int *in_which_cell_P; /*tells where a given particle resides*/
	int *cell_has_how_much_Pparticle; /*self-explanatory*/
	int *where_in_cell_list_P; /*facilitates cell reconstruction*/
	cell_has_which_Pparticle = (int *)calloc(c_index_max*n_cell_max, sizeof(int));
	in_which_cell_P = (int *)calloc(N_pinned_tot+1, sizeof(int));
	cell_has_how_much_Pparticle = (int *)calloc(c_index_max, sizeof(int));
	where_in_cell_list_P = (int *)calloc(N_pinned_tot+1, sizeof(int));
	if (cell_has_which_Pparticle==NULL || in_which_cell_P==NULL || cell_has_how_much_Pparticle==NULL || where_in_cell_list_P==NULL){printf("calloc of cell structure failed. \n");exit(EXIT_FAILURE);}
	
	/*Allocate memory for energy list*/
	double *etot_list; /*half-energy associated with the particular particle*/
	etot_list=(double *)calloc(N_mobile_tot+N_pinned_tot+1, sizeof(double));
	if (etot_list==NULL){printf("calloc of neighbor list structure failed. \n");exit(EXIT_FAILURE);}
	
    double *etot_list_checker; /*half-energy associated with the particular particle*/
    etot_list_checker=(double *)calloc(N_mobile_tot+N_pinned_tot+1, sizeof(double));
    if (etot_list_checker==NULL){printf("calloc of neighbor list structure failed. \n");exit(EXIT_FAILURE);}
	
////////////////Neighbor search facilitating device////////////////
	
	int *neighbor_search_facilitator_x; /*given neighbor_cell_index, get the associated index in x direction*/
	int *neighbor_search_facilitator_y; /*given neighbor_cell_index, get the associated index in y direction*/
	int *neighbor_search_facilitator_z; /*given neighbor_cell_index, get the associated index in z direction*/
	neighbor_search_facilitator_x=(int *)calloc(neighbor_cell_index_max, sizeof(int));
	neighbor_search_facilitator_y=(int *)calloc(neighbor_cell_index_max, sizeof(int));
	neighbor_search_facilitator_z=(int *)calloc(neighbor_cell_index_max, sizeof(int));
	if (neighbor_search_facilitator_x==NULL || neighbor_search_facilitator_y==NULL || neighbor_search_facilitator_z==NULL){printf("calloc of search facilitator failed. \n");exit(EXIT_FAILURE);}
	
	for (x_neighbor_search=-1;x_neighbor_search<2;x_neighbor_search++){
		for (y_neighbor_search=-1;y_neighbor_search<2;y_neighbor_search++){
			for (z_neighbor_search=-1;z_neighbor_search<2;z_neighbor_search++){
				neighbor_cell_index=(x_neighbor_search+1)+(y_neighbor_search+1)*3+(z_neighbor_search+1)*9;
				neighbor_search_facilitator_x[neighbor_cell_index]=x_neighbor_search;
				neighbor_search_facilitator_y[neighbor_cell_index]=y_neighbor_search;
				neighbor_search_facilitator_z[neighbor_cell_index]=z_neighbor_search;
			}
		}
	}
	
////////////////Assign a suitable cell structure////////////////

	/*Inside: evolve in time*/
	for (i=1;i<(N_mobile_tot+1);i++){
		x_temp=x_M[i];
		y_temp=y_M[i];
		z_temp=z_M[i];
		c_x=floor((x_temp+hL)/l_cell);
		c_y=floor((y_temp+hL)/l_cell);
		c_z=floor((z_temp+hL)/l_cell);
		c_index=c_x+c_y*n+c_z*(n*n);
		n_particles_in_cell=cell_has_how_much_Mparticle[c_index];
		index_within_cell=n_particles_in_cell;
		cell_has_which_Mparticle[c_index*n_cell_max+index_within_cell]=i;
		where_in_cell_list_M[i]=index_within_cell;
		n_particles_in_cell++;
		if (n_particles_in_cell>n_cell_max){printf("Play safer for n_cell_max, mobile! \n");exit(EXIT_FAILURE);}
		cell_has_how_much_Mparticle[c_index]=n_particles_in_cell;
		in_which_cell_M[i]=c_index;
	}
	
	/*Outside: never change*/
	for (i=1;i<(N_pinned_tot+1);i++){
		x_temp=x_P[i];
		y_temp=y_P[i];
		z_temp=z_P[i];
		c_x=floor((x_temp+hL)/l_cell);
		c_y=floor((y_temp+hL)/l_cell);
		c_z=floor((z_temp+hL)/l_cell);
		c_index=c_x+c_y*n+c_z*(n*n);
		n_particles_in_cell=cell_has_how_much_Pparticle[c_index];
		index_within_cell=n_particles_in_cell;
		cell_has_which_Pparticle[c_index*n_cell_max+index_within_cell]=i;
		where_in_cell_list_P[i]=index_within_cell;
		n_particles_in_cell++;
		if (n_particles_in_cell>n_cell_max){printf("Play safer for n_cell_max, pinned! \n");exit(EXIT_FAILURE);}
		cell_has_how_much_Pparticle[c_index]=n_particles_in_cell;
		in_which_cell_P[i]=c_index;
	}
	
////////////////Keep record of the reference configuration and parameters////////////////
	
	sprintf(filename,"./CavityRef.dat");
	if((outfile = fopen(filename, "w+")) == NULL){printf("Error opening file %s\n",filename);exit(EXIT_FAILURE);}
    fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\n",beta_string[0],IRcutoff,cavity_size);
    
	fprintf(outfile,"%d\t%d\n",NA_pinned,NB_pinned);
	for(i=1;i<N_pinned_tot+1;i++){
		fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\t%d\n",x_P[i],y_P[i],z_P[i],a_P[i]);
	}
	
	fprintf(outfile,"%d\t%d\n",NA_mobile,NB_mobile);
	for(i=1;i<N_mobile_tot+1;i++){
		fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\t%d\n",x_M[i],y_M[i],z_M[i],a_M[i]);
	}
	fclose(outfile);

////////////////Mapping of indices////////////////
	
	int *from_thermalI_to_repI;
	int *from_repI_to_thermalI;
	from_thermalI_to_repI = (int *)calloc(n_replica, sizeof(int));
	from_repI_to_thermalI = (int *)calloc(n_replica, sizeof(int));
    if (from_thermalI_to_repI==NULL || from_repI_to_thermalI==NULL){printf("calloc of index-mappers failed. \n");exit(EXIT_FAILURE);}
	
	
////////////////Allocate memory for REPLICAS////////////////
	
	double *x_M_replicated;
	double *y_M_replicated;
	double *z_M_replicated;
	x_M_replicated = (double *)calloc((N_mobile_tot+1)*n_replica, sizeof(double));
	y_M_replicated = (double *)calloc((N_mobile_tot+1)*n_replica, sizeof(double));
	z_M_replicated = (double *)calloc((N_mobile_tot+1)*n_replica, sizeof(double));
	if (x_M_replicated==NULL || y_M_replicated==NULL || z_M_replicated==NULL){printf("calloc of coordinates failed. \n");exit(EXIT_FAILURE);}
	
	int *cell_has_which_Mparticle_replicated; /* 0 means empty*/
	int *in_which_cell_M_replicated;
	int *cell_has_how_much_Mparticle_replicated;
	int *where_in_cell_list_M_replicated;
	cell_has_which_Mparticle_replicated = (int *)calloc((c_index_max*n_cell_max)*n_replica, sizeof(int));
	in_which_cell_M_replicated = (int *)calloc((N_mobile_tot+1)*n_replica, sizeof(int));
	cell_has_how_much_Mparticle_replicated = (int *)calloc((c_index_max)*n_replica, sizeof(int));
	where_in_cell_list_M_replicated = (int *)calloc((N_mobile_tot+1)*n_replica, sizeof(int));
	if (cell_has_which_Mparticle_replicated==NULL || in_which_cell_M_replicated==NULL || cell_has_how_much_Mparticle_replicated==NULL|| where_in_cell_list_M_replicated==NULL){printf("calloc of cell structure failed. \n");exit(EXIT_FAILURE);}
	
	double *etot_list_replicated; /*half-energy associated with the particular particle*/
	etot_list_replicated=(double *)calloc((N_mobile_tot+N_pinned_tot+1)*n_replica, sizeof(double));
	if (etot_list_replicated==NULL){printf("calloc of neighbor list structure failed. \n");exit(EXIT_FAILURE);}
	
	
	
	
	
	
////////////////Initialization: Xerox the original configuration////////////////
	
	for (replica_index=0;replica_index<n_replica;replica_index++){
		Put_One_Replica_Back(replica_index,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
	}
	
	for (replica_index=0;replica_index<n_replica;replica_index++){
		from_thermalI_to_repI[replica_index]=replica_index;
		from_repI_to_thermalI[replica_index]=replica_index;
	}


////////////////Re-calculate the energy profile with the shrinking factor before beggining////////////////
	
	for (thermal_index=0;thermal_index<n_replica;thermal_index++){
		replica_index=from_thermalI_to_repI[thermal_index];
		lambda=lambda_string[thermal_index];
		Take_One_Replica_Out(replica_index,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
	
		/*Update energy profile and record it*/
		Energy_Profiler(lambda,N_mobile_tot,N_pinned_tot,NA_mobile,NA_pinned,n,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);
		for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
			etot_list_replicated[(i)*(n_replica)+replica_index]=etot_list[i];
		}
	}	
	
	

	
	
	
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////Now we are ready to begin the sequence of cavity simulations!////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
	
	
	counter=0;
	snap_index=0;
	count_down_till_save=0;
	
	while (counter<NEquiPro){
		/*Determine how much macro-MC steps we take before snapshot, swap-trial, and/or energy-adjustment*/
		macroMC_step=0;
		SNAPSHOT_FLAG=FALSE;
		SWAP_TRIAL_FLAG=FALSE;
		ENERGY_ADJUSTMENT_FLAG=FALSE;
		SAVE_FLAG=FALSE;
		while(SNAPSHOT_FLAG==FALSE && SWAP_TRIAL_FLAG==FALSE && ENERGY_ADJUSTMENT_FLAG==FALSE){
			counter++;
			macroMC_step++;
			if ((counter%Ndump)==0){SNAPSHOT_FLAG=TRUE;}
			if (RandomNumber()<p_swap){SWAP_TRIAL_FLAG=TRUE;}
			if ((counter%Nadjust)==0){ENERGY_ADJUSTMENT_FLAG=TRUE;}
		}
		
		

////////////////Actual Equilibration Part////////////////

		for (thermal_index=0;thermal_index<n_replica;thermal_index++){
			replica_index=from_thermalI_to_repI[thermal_index];
			beta=beta_string[thermal_index];
			lambda=lambda_string[thermal_index];
			
			/*Replicate forward*/
			Take_One_Replica_Out(replica_index,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);

			/*Actual actual equilibration*/		
			
			Cavity_Equilibration_MC(macroMC_step,lambda,hard_wall_size,N_mobile_tot,NA_mobile,NA_pinned,max_move_length,n,hL,l_cell,beta,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,where_in_cell_list_M,where_in_cell_list_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);

			/*Replicated back*/
			Put_One_Replica_Back(replica_index,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
		
			/*SNAPSHOT*/
			if (SNAPSHOT_FLAG==TRUE && thermal_index==0){
				sprintf(filename,"./CavitySnap%lld.dat",snap_index);
				if((outfile = fopen(filename, "w+")) == NULL){printf("Error opening file CavitySnap.dat\n");exit(EXIT_FAILURE);}
				fprintf(outfile,"%lld\n",counter);
				for(i=1;i<N_mobile_tot+1;i++){
					fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\t%d\n",x_M[i],y_M[i],z_M[i],a_M[i]);
				}
				fclose(outfile);
				snap_index++;
				count_down_till_save++;
				if (count_down_till_save==1000){SAVE_FLAG=TRUE; count_down_till_save=0;}
			}
		}
		
////////////////POSSIBLE SWAPPING////////////////
		
		if (SWAP_TRIAL_FLAG==TRUE){
			SWAP_FLAG=FALSE;
			
			/*FLIP THE COIN 1: randomly choose which two adjacent temperatures to consider*/
			thermal_index_low=(int)(RandomNumber()*(n_replica-1));
			thermal_index_high=thermal_index_low+1;
			beta_low=beta_string[thermal_index_low]; /*``LOW": temperature is low, beta is high*/
			beta_high=beta_string[thermal_index_high]; /*``HIGH": temperature is high, beta is low*/
			lambda_low=lambda_string[thermal_index_low]; /*``LOW": lambda is big*/
			lambda_high=lambda_string[thermal_index_high]; /*``HIGH": lambda is small;shrunk more*/
			
			/*Get energies 1: I MUST include Frozen particle contirbution here (this can be made more efficient by just keeping track of things in MC)*/
			replica_index_low=from_thermalI_to_repI[thermal_index_low];
			Etot_low_big=0;
			for (i=1;i<(N_mobile_tot+N_pinned_tot+1);i++){
				Etot_low_big=Etot_low_big+etot_list_replicated[(i)*(n_replica)+replica_index_low];
			}
			
			replica_index_high=from_thermalI_to_repI[thermal_index_high];
			Etot_high_small=0;
			for (i=1;i<(N_mobile_tot+N_pinned_tot+1);i++){
				Etot_high_small=Etot_high_small+etot_list_replicated[(i)*(n_replica)+replica_index_high];
			}
			
			/*Get energies 2, with SWAPPED HAMILTONIAN: I MUST include Frozen particle contirbution here*/
			replica_index_low=from_thermalI_to_repI[thermal_index_low];
			Take_One_Replica_Out(replica_index_low,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
			Energy_Profiler(lambda_high,N_mobile_tot,N_pinned_tot,NA_mobile,NA_pinned,n,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);
			Etot_low_small=0;
			for (i=1;i<(N_mobile_tot+N_pinned_tot+1);i++){
				Etot_low_small=Etot_low_small+etot_list[i];
			}
			
			replica_index_high=from_thermalI_to_repI[thermal_index_high];
			Take_One_Replica_Out(replica_index_high,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
			Energy_Profiler(lambda_low,N_mobile_tot,N_pinned_tot,NA_mobile,NA_pinned,n,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);
			Etot_high_big=0;
			for (i=1;i<(N_mobile_tot+N_pinned_tot+1);i++){
				Etot_high_big=Etot_high_big+etot_list[i];
			}
			
			/*FLIP THE COIN 2*/
			swap_decider=beta_low*Etot_low_big-beta_high*Etot_low_small-beta_low*Etot_high_big+beta_high*Etot_high_small;
			if (swap_decider>=0){
				SWAP_FLAG=TRUE;
			}
			if(swap_decider<0){
				if(RandomNumber()<exp(swap_decider)){
					SWAP_FLAG=TRUE;
				}
			}
			
			/*If the swap is accepted, then swap*/
			if (SWAP_FLAG==TRUE){
				from_thermalI_to_repI[thermal_index_low]=replica_index_high;
				from_thermalI_to_repI[thermal_index_high]=replica_index_low;
				from_repI_to_thermalI[replica_index_low]=thermal_index_high;
				from_repI_to_thermalI[replica_index_high]=thermal_index_low;
				
				for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
					etot_list_replicated[(i)*(n_replica)+replica_index_high]=etot_list[i];
				}
				Take_One_Replica_Out(replica_index_low,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
				Energy_Profiler(lambda_high,N_mobile_tot,N_pinned_tot,NA_mobile,NA_pinned,n,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);
				for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
					etot_list_replicated[(i)*(n_replica)+replica_index_low]=etot_list[i];
				}
			}
		}
		
		
////////////////Energy Adjustment////////////////		
		if (ENERGY_ADJUSTMENT_FLAG==TRUE){
			for (thermal_index=0;thermal_index<n_replica;thermal_index++){
				replica_index=from_thermalI_to_repI[thermal_index];
				lambda=lambda_string[thermal_index];
				Take_One_Replica_Out(replica_index,N_mobile_tot,N_pinned_tot,n_replica,c_index_max,x_M,y_M,z_M,in_which_cell_M,cell_has_which_Mparticle,where_in_cell_list_M,cell_has_how_much_Mparticle,etot_list,x_M_replicated,y_M_replicated,z_M_replicated,in_which_cell_M_replicated,cell_has_which_Mparticle_replicated,where_in_cell_list_M_replicated,cell_has_how_much_Mparticle_replicated,etot_list_replicated);
				
				for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
					etot_list_checker[i]=etot_list[i];
				}
				/*Update energy profile and record it*/
				Energy_Profiler(lambda,N_mobile_tot,N_pinned_tot,NA_mobile,NA_pinned,n,etot_list,x_M,y_M,z_M,a_M,x_P,y_P,z_P,a_P,cell_has_which_Mparticle,cell_has_which_Pparticle,in_which_cell_M,in_which_cell_P,cell_has_how_much_Mparticle,cell_has_how_much_Pparticle,neighbor_search_facilitator_x,neighbor_search_facilitator_y,neighbor_search_facilitator_z,interaction_range);
				for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
					etot_list_replicated[(i)*(n_replica)+replica_index]=etot_list[i];
					energy_checker=etot_list_replicated[(i)*(n_replica)+replica_index]-etot_list_checker[i];
					if (energy_checker>0.000001 || energy_checker<-0.000001){printf("Adjust energy more often!\n");exit(EXIT_FAILURE);}					
				}
			}	
		}
		
////////////////Save the restart file once in a while////////////////
		if (SAVE_FLAG==TRUE){
			
			/*Restart file*/
			sprintf(filename,"./Restart%lld.dat",counter);
			if((outfile = fopen(filename, "w+")) == NULL){printf("Error opening file ReStart.dat\n");exit(EXIT_FAILURE);}
			fprintf(outfile,"%4.17f\t%4.17f\n",cavity_size,IRcutoff);
			fprintf(outfile,"%lld\t%lld\n",Ndump,Nsnap);
			fprintf(outfile,"%lld\t%lld\n",snap_index,counter);
			fprintf(outfile,"%d\t%d\n",NA_pinned,NB_pinned);
			for(i=1;i<N_pinned_tot+1;i++){
				fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\t%d\n",x_P[i],y_P[i],z_P[i],a_M[i]);
			}
			
			fprintf(outfile,"%d\t%d\t%d\n",NA_mobile,NB_mobile,n_replica);
			for (thermal_index=0;thermal_index<n_replica;thermal_index++){
				replica_index=from_thermalI_to_repI[thermal_index];
				fprintf(outfile,"%d\t%4.17f\t%4.17f\n",replica_index,beta_string[thermal_index],lambda_string[thermal_index]);
				for(i=1;i<N_mobile_tot+1;i++){
					fprintf(outfile,"%4.17f\t%4.17f\t%4.17f\t%d\n",x_M_replicated[(i)*(n_replica)+replica_index],y_M_replicated[(i)*(n_replica)+replica_index],z_M_replicated[(i)*(n_replica)+replica_index],a_M[i]);
				}
			}
			fclose(outfile);
		}
	}

	
////////Let there be freedom for memory!////////
	

	free(x_M);
	free(y_M);
	free(z_M);
    free(a_M);
	free(x_P);
	free(y_P);
	free(z_P);
    free(a_P);
	free(cell_has_which_Mparticle);
	free(cell_has_how_much_Mparticle);
	free(in_which_cell_M);
	free(where_in_cell_list_M);
	free(cell_has_which_Pparticle);
	free(cell_has_how_much_Pparticle);
	free(in_which_cell_P);
	free(where_in_cell_list_P);;
	free(cell_has_which_Mparticle_replicated);
	free(cell_has_how_much_Mparticle_replicated);
	free(in_which_cell_M_replicated);
	free(where_in_cell_list_M_replicated);
	free(x_M_replicated);
	free(y_M_replicated);
	free(z_M_replicated);
	free(etot_list_replicated);
	free(etot_list);
	free(neighbor_search_facilitator_x);
	free(neighbor_search_facilitator_y);
	free(neighbor_search_facilitator_z);
	
	free(beta_string);
	free(lambda_string);
    free(from_thermalI_to_repI);
	free(from_repI_to_thermalI);
	
	
	
	
////////////////////////////////////////
////////////////THE END!////////////////
////////////////////////////////////////
	return 1;
}






////////////////One trial move per particle////////////////
void Cavity_Equilibration_MC(long long int macroMC_step,double lambda,double hard_wall_size,int N_mobile_tot,int NA_mobile,int NA_pinned,double max_move_length,int n,double hL, double l_cell, double beta,double etot_list[],double x_M[],double y_M[],double z_M[],int a_M[],double x_P[],double y_P[],double z_P[],int a_P[],int cell_has_which_Mparticle[],int cell_has_which_Pparticle[],int in_which_cell_M[],int in_which_cell_P[],int where_in_cell_list_M[],int where_in_cell_list_P[],int cell_has_how_much_Mparticle[],int cell_has_how_much_Pparticle[],int neighbor_search_facilitator_x[],int neighbor_search_facilitator_y[],int neighbor_search_facilitator_z[],double interaction_range){

////////////////HARD-CODED Binary Lennard-Jones liquid parameters (see, e.g., J. Chem. Phys. 144, 024501 (2016))////////////////
    double sigmaAA=(double)1.00;
    double sigmaAB=(double)0.80;
    double sigmaBB=(double)0.88;
    double eAA=(double)1.0;
    double eAB=(double)1.5;
    double eBB=(double)0.5;
    double rcutAA=(double)2.5;
    double rcutAB=(double)2.0;
    double rcutBB=(double)2.2;
    /*Quickly cook Binary Lennard-Jones liquid parameters*/
    double sigmaAA_sixth=CUBE(SQR(sigmaAA));
    double sigmaAB_sixth=CUBE(SQR(sigmaAB));
    double sigmaBB_sixth=CUBE(SQR(sigmaBB));
    double sigmaAA_twelfth=SQR(sigmaAA_sixth);
    double sigmaAB_twelfth=SQR(sigmaAB_sixth);
    double sigmaBB_twelfth=SQR(sigmaBB_sixth);
    double rcutAA_square=SQR(rcutAA);
    double rcutAB_square=SQR(rcutAB);
    double rcutBB_square=SQR(rcutBB);
    double r_sixth;
    double r_twelfth;
    double VshiftAA;
    double VshiftAB;
    double VshiftBB;
    if (rcutAA>interaction_range || rcutAB>interaction_range || rcutBB>interaction_range){printf("Interaction range chosen too small\n");exit(EXIT_FAILURE);}
    r_sixth=CUBE(rcutAA_square);
    r_twelfth=SQR(r_sixth);
    VshiftAA=-((double)2.0)*eAA*(-(sigmaAA_sixth/r_sixth)+(sigmaAA_twelfth/r_twelfth));
    r_sixth=CUBE(rcutAB_square);
    r_twelfth=SQR(r_sixth);
    VshiftAB=-((double)2.0)*eAB*(-(sigmaAB_sixth/r_sixth)+(sigmaAB_twelfth/r_twelfth));
    r_sixth=CUBE(rcutBB_square);
    r_twelfth=SQR(r_sixth);
    VshiftBB=-((double)2.0)*eBB*(-(sigmaBB_sixth/r_sixth)+(sigmaBB_twelfth/r_twelfth));
	
////////////////Initial assignment of variables////////////////
	
	long long int t;
	int ACCEPT_FLAG;
	int WALL_FLAG;
	int i_move;
	int i_neighbor;
    int a_given;
    int a_neighbor;
	int i_shuffle;
	double l_move;
	double pi_double=2*M_PI;
	double costheta;
	double sintheta; /*always positive in the standard convention*/
	double phi;
	double sinphi;
	double cosphi;
	double x_new;
	double y_new;
	double z_new;
	double r_new;
	double x_ori;
	double y_ori;
	double z_ori;
	double x_temp;
	double y_temp;
	double z_temp;
	double radial_distance_square;
	double half_pair_energy;
	double energy_before;
	double energy_after;
	double energy_change;
	int c_index_ori; /*cell index*/
	int c_x_ori; /*cell index in x direction*/
	int c_y_ori; /*cell index in y direction*/
	int c_z_ori; /*cell index in z direction*/
	int c_index_new; /*cell index*/
	int c_x_new; /*cell index in x direction*/
	int c_y_new; /*cell index in y direction*/
	int c_z_new; /*cell index in z direction*/
	int c_index_neighbor; /*cell index*/
	int c_x_neighbor; /*cell index in x direction*/
	int c_y_neighbor; /*cell index in y direction*/
	int c_z_neighbor; /*cell index in z direction*/
	int Mneighbor_counter_new;
	int Mneighbor_counter_ori;
	int Pneighbor_counter_new;
	int Pneighbor_counter_ori;
	int relabelling_within_cell_delete;
	int relabelling_within_cell_shuffle;
	int n_cell_max=NCELLMAX; /*n_cell_max sets the maximum number of particles in one cell allowed*/
	
	int index_within_cell;
	int n_particles_in_cell;
	int neighbor_index;
	int neighbor_index_max=NEIGHBORMAX;
	int neighbor_cell_index;
	int neighbor_cell_index_max=DIMPOWER(3);
	

	double shrinking_modificationUU=1/(SQR(lambda)); /*When two particles are both unpinned and shrunk*/
	double shrinking_modificationUP=1/(SQR((1+lambda)/2)); /*When one particle is unpinned and shrunk*/
	
	double sigma_sixth_temp;
	double sigma_twelfth_temp;
	double Vshift_temp;
	double e_temp;
	double rcut_square_temp;

	
////////////////To reduce the repetitive energy calculations, keep the record (can be improved)////////////////

	double *temporal_energy_list_M_after;
	double *temporal_energy_list_M_before;
	int *temporal_neighbor_list_M_after;
	int *temporal_neighbor_list_M_before;
	temporal_energy_list_M_after=(double *)calloc(neighbor_index_max, sizeof(double));
	temporal_energy_list_M_before=(double *)calloc(neighbor_index_max, sizeof(double));
	temporal_neighbor_list_M_after=(int *)calloc(neighbor_index_max, sizeof(int));
	temporal_neighbor_list_M_before=(int *)calloc(neighbor_index_max, sizeof(int));
	if (temporal_energy_list_M_after==NULL || temporal_energy_list_M_before==NULL || temporal_neighbor_list_M_after==NULL || temporal_neighbor_list_M_before==NULL){printf("calloc of e_change_list failed. \n");exit(EXIT_FAILURE);}

	double *temporal_energy_list_P_after;
	double *temporal_energy_list_P_before;
	int *temporal_neighbor_list_P_after;
	int *temporal_neighbor_list_P_before;
	temporal_energy_list_P_after=(double *)calloc(neighbor_index_max, sizeof(double));
	temporal_energy_list_P_before=(double *)calloc(neighbor_index_max, sizeof(double));
	temporal_neighbor_list_P_after=(int *)calloc(neighbor_index_max, sizeof(int));
	temporal_neighbor_list_P_before=(int *)calloc(neighbor_index_max, sizeof(int));
	if (temporal_energy_list_P_after==NULL || temporal_energy_list_P_before==NULL || temporal_neighbor_list_P_after==NULL || temporal_neighbor_list_P_before==NULL){printf("calloc of e_change_list failed. \n");exit(EXIT_FAILURE);}
	
	
////////////////Now begin N_mobile_tot trial moves////////////////
	
	for (t=0;t<N_mobile_tot*macroMC_step;t++){
		
		ACCEPT_FLAG=FALSE;
		WALL_FLAG=FALSE;

////////////////Random Move Choice////////////////
		
		/*Pick particle at random*/
		i_move = (int)(RandomNumber()*N_mobile_tot)+1;
        a_given=a_M[i_move];
        
		/*Pick length at random*/
		l_move=RandomNumber()*max_move_length;
		
		/*Pick direction at random (uniformly on unit sphere)*/
		costheta=(2.0)*RandomNumber()-1;
		sintheta=sqrt(1-SQR(costheta));
		
		phi=RandomNumber()*pi_double;
		cosphi=cos(phi);
		sinphi=sin(phi);

		
////////////////Phase 1: a particle is (possibly) moving////////////////
		
		/*Prepare yourself*/
		energy_before=2*etot_list[i_move];
		energy_after=0;
		Mneighbor_counter_new=0;
		Pneighbor_counter_new=0;
		
		/*Where was it originally?*/
		x_ori=x_M[i_move];
		y_ori=y_M[i_move];
		z_ori=z_M[i_move];
		c_index_ori=in_which_cell_M[i_move];
		
		/*Try a move*/
		x_new=x_ori+l_move*sintheta*cosphi;
		y_new=y_ori+l_move*sintheta*sinphi;
		z_new=z_ori+l_move*costheta;
		r_new=sqrt(SQR(x_new)+SQR(y_new)+SQR(z_new));
		if (r_new>=hard_wall_size){
			WALL_FLAG=TRUE;
		}
		
		/*In which cell is it in now?*/
		c_x_new=floor((x_new+hL)/l_cell);
		c_y_new=floor((y_new+hL)/l_cell);
		c_z_new=floor((z_new+hL)/l_cell);
		c_index_new=c_x_new+c_y_new*n+c_z_new*(n*n);
			
		
////////////////Phase 2: calculate the new energy and FLIP THE COIN////////////////
	
		if (WALL_FLAG==FALSE){
			
			/*Look around neighbor cells of the NEW cell*/
			for (neighbor_cell_index=0;neighbor_cell_index<neighbor_cell_index_max;neighbor_cell_index++){
				c_x_neighbor=c_x_new+neighbor_search_facilitator_x[neighbor_cell_index];
				c_y_neighbor=c_y_new+neighbor_search_facilitator_y[neighbor_cell_index];
				c_z_neighbor=c_z_new+neighbor_search_facilitator_z[neighbor_cell_index];
				c_index_neighbor=c_x_neighbor+c_y_neighbor*n+c_z_neighbor*(n*n);
				
				/*Pick mobile-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Mparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Mparticle[c_index_neighbor*n_cell_max+index_within_cell];
                        a_neighbor=a_M[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                        else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                        else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                        }
						
						/*Calculate pair distance*/
						/*NO SELF ENERGY*/
						if (i_neighbor!=i_move){
							x_temp=x_new-x_M[i_neighbor];
							y_temp=y_new-y_M[i_neighbor];
							z_temp=z_new-z_M[i_neighbor];
							radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
							
							/*Possibly shrunk*/
							radial_distance_square=radial_distance_square*shrinking_modificationUU;
							
							if (radial_distance_square<rcut_square_temp){
								r_sixth=CUBE(radial_distance_square);
								r_twelfth=SQR(r_sixth);
								
								/*Calculate pair LJ energy, ``half of it"*/
								half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;

								/*Add it and count it*/
								energy_after=energy_after+2*half_pair_energy;
								temporal_energy_list_M_after[Mneighbor_counter_new]=half_pair_energy;
								temporal_neighbor_list_M_after[Mneighbor_counter_new]=i_neighbor;
								Mneighbor_counter_new++;
								if (Mneighbor_counter_new>neighbor_index_max){printf("Play safer for neighbor_index_max, mobile! \n");exit(EXIT_FAILURE);}	
							}
						}					
					}
				}
				
				/*Pick pinned-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Pparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Pparticle[c_index_neighbor*n_cell_max+index_within_cell];
						a_neighbor=a_P[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                        else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                        else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                        }
						
						/*Calculate pair distance*/
						x_temp=x_new-x_P[i_neighbor];
						y_temp=y_new-y_P[i_neighbor];
						z_temp=z_new-z_P[i_neighbor];
						radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
						
						/*Possibly shrunk*/
						radial_distance_square=radial_distance_square*shrinking_modificationUP;
						
						if (radial_distance_square<rcut_square_temp){
							r_sixth=CUBE(radial_distance_square);
							r_twelfth=SQR(r_sixth);
							
							/*Calculate pair LJ energy, ``half of it"*/
							half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;
							
							/*Add it and count it*/
							energy_after=energy_after+2*half_pair_energy;
							temporal_energy_list_P_after[Pneighbor_counter_new]=half_pair_energy;
							temporal_neighbor_list_P_after[Pneighbor_counter_new]=i_neighbor;
							Pneighbor_counter_new++;
							if (Pneighbor_counter_new>neighbor_index_max){printf("Play safer for neighbor_index_max, pinned! \n");exit(EXIT_FAILURE);}
						}
					}
				}
			}

			/*FLIP THE COIN*/
			energy_change=energy_after-energy_before;
			if (energy_change<=0){
				ACCEPT_FLAG=TRUE;
			}
			if(energy_change>0){
				if(RandomNumber()<exp(-beta*energy_change)){
					ACCEPT_FLAG=TRUE;
				}
			}

			
////////////////Phase 3: when move is accepted, calculate the original energy (THIS CAN BE IMPROVED by keeping the neighbor-list, at the cost of big data storage: I opt for the current scheme)////////////////
				
			if (ACCEPT_FLAG==TRUE){
				Mneighbor_counter_ori=0;
				Pneighbor_counter_ori=0;
				
				/*Look around neighbor cells of the ORIGINAL cell*/
				c_x_ori=c_index_ori%n;
				c_y_ori=((int)floor(c_index_ori/((double)n)))%n;	
				c_z_ori=floor(c_index_ori/(((double)n)*n));
				for (neighbor_cell_index=0;neighbor_cell_index<neighbor_cell_index_max;neighbor_cell_index++){
					c_x_neighbor=c_x_ori+neighbor_search_facilitator_x[neighbor_cell_index];
					c_y_neighbor=c_y_ori+neighbor_search_facilitator_y[neighbor_cell_index];
					c_z_neighbor=c_z_ori+neighbor_search_facilitator_z[neighbor_cell_index];
					c_index_neighbor=c_x_neighbor+c_y_neighbor*n+c_z_neighbor*(n*n);
					
					/*Pick mobile-neighbors, one by one*/
					n_particles_in_cell=cell_has_how_much_Mparticle[c_index_neighbor];
					if (n_particles_in_cell!=0){
						for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
							i_neighbor=cell_has_which_Mparticle[c_index_neighbor*n_cell_max+index_within_cell];
                            a_neighbor=a_M[i_neighbor];
                            
                            if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                            else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                            else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                            }
							
							/*Calculate pair distance*/
							/*NO SELF ENERGY*/
							if (i_neighbor!=i_move){
								x_temp=x_ori-x_M[i_neighbor];
								y_temp=y_ori-y_M[i_neighbor];
								z_temp=z_ori-z_M[i_neighbor];
								radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
								
								/*Possibly shrunk*/
								radial_distance_square=radial_distance_square*shrinking_modificationUU;
								
								if (radial_distance_square<rcut_square_temp){
									r_sixth=CUBE(radial_distance_square);
									r_twelfth=SQR(r_sixth);
									
									/*Calculate pair LJ energy, ``half of it"*/
									half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;

									/*Add it and count it*/
									temporal_energy_list_M_before[Mneighbor_counter_ori]=half_pair_energy;
									temporal_neighbor_list_M_before[Mneighbor_counter_ori]=i_neighbor;
									Mneighbor_counter_ori++;
									if (Mneighbor_counter_ori>neighbor_index_max){printf("Play safer for neighbor_index_max, mobile! \n");exit(EXIT_FAILURE);}
								}
							}					
						}
					}
					
					/*Pick pinned-neighbors, one by one*/
					n_particles_in_cell=cell_has_how_much_Pparticle[c_index_neighbor];
					if (n_particles_in_cell!=0){
						for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
							i_neighbor=cell_has_which_Pparticle[c_index_neighbor*n_cell_max+index_within_cell];
                            a_neighbor=a_P[i_neighbor];
                            
                            if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                            else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                            else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                            }
							
							/*Calculate pair distance*/
							x_temp=x_ori-x_P[i_neighbor];
							y_temp=y_ori-y_P[i_neighbor];
							z_temp=z_ori-z_P[i_neighbor];
							radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
							
							/*Possibly shrunk*/
							radial_distance_square=radial_distance_square*shrinking_modificationUP;
							
							if (radial_distance_square<rcut_square_temp){
								r_sixth=CUBE(radial_distance_square);
								r_twelfth=SQR(r_sixth);
								
								/*Calculate pair LJ energy, ``half of it"*/
								half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;

								/*Add it and count it*/
								temporal_energy_list_P_before[Pneighbor_counter_ori]=half_pair_energy;
								temporal_neighbor_list_P_before[Pneighbor_counter_ori]=i_neighbor;
								Pneighbor_counter_ori++;
								if (Pneighbor_counter_ori>neighbor_index_max){printf("Play safer for neighbor_index_max, pinned! \n");exit(EXIT_FAILURE);}
							}
						}
					}
				}
					
////////////////Phase 4: when move is accepted, update configuration, cell strucutre, and an energy profile////////////////
					
				/*Update configuration*/
				x_M[i_move]=x_new;
				y_M[i_move]=y_new;
				z_M[i_move]=z_new;
				
				/*Update cell structure, if necessary*/
				/*Order is very important in the following: do not mess with it*/
				if (c_index_new!=c_index_ori){
					
					/*delete from cell_has_which_Mparticle*/
					n_particles_in_cell=cell_has_how_much_Mparticle[c_index_ori];
					index_within_cell=where_in_cell_list_M[i_move];
					relabelling_within_cell_delete=c_index_ori*n_cell_max+index_within_cell;
					relabelling_within_cell_shuffle=c_index_ori*n_cell_max+n_particles_in_cell-1;
					i_shuffle=cell_has_which_Mparticle[relabelling_within_cell_shuffle];
					cell_has_which_Mparticle[relabelling_within_cell_delete]=i_shuffle;
					cell_has_which_Mparticle[relabelling_within_cell_shuffle]=EMPTY;
					
					/*update where_in_cell_list_M (where, without it, a bug once lived)*/
					where_in_cell_list_M[i_shuffle]=index_within_cell;
					
					/*delete from cell_has_how_much_Mparticle*/
					cell_has_how_much_Mparticle[c_index_ori]=n_particles_in_cell-1;
					
					/*update in_which_cell_M*/
					in_which_cell_M[i_move]=c_index_new;
					
					/*add to cell_has_which_Mparticle*/
					n_particles_in_cell=cell_has_how_much_Mparticle[c_index_new];
					index_within_cell=n_particles_in_cell;
					cell_has_which_Mparticle[c_index_new*n_cell_max+index_within_cell]=i_move;
					
					/*update where_in_cell_list_M*/
					where_in_cell_list_M[i_move]=index_within_cell;	
					
					/*add cell_has_how_much_Mparticle*/
					n_particles_in_cell++;
					if (n_particles_in_cell>n_cell_max){printf("Play safer for n_cell_max, mobile! \n");exit(EXIT_FAILURE);}
					cell_has_how_much_Mparticle[c_index_new]=n_particles_in_cell;
				}
				

				/*Subtraction MM*/
				if (Mneighbor_counter_ori!=0){
					for (neighbor_index=0;neighbor_index<Mneighbor_counter_ori;neighbor_index++){
						i_neighbor=temporal_neighbor_list_M_before[neighbor_index];
						half_pair_energy=temporal_energy_list_M_before[neighbor_index];
						etot_list[i_move]=etot_list[i_move]-half_pair_energy;
						etot_list[i_neighbor]=etot_list[i_neighbor]-half_pair_energy;
					}
				}
				
				/*Subtraction MP*/
				if (Pneighbor_counter_ori!=0){
					for (neighbor_index=0;neighbor_index<Pneighbor_counter_ori;neighbor_index++){
						i_neighbor=temporal_neighbor_list_P_before[neighbor_index];
						half_pair_energy=temporal_energy_list_P_before[neighbor_index];
						etot_list[i_move]=etot_list[i_move]-half_pair_energy;
						etot_list[N_mobile_tot+i_neighbor]=etot_list[N_mobile_tot+i_neighbor]-half_pair_energy;
					}
				}
				
				/*Addition MM*/
				if (Mneighbor_counter_new!=0){
					for (neighbor_index=0;neighbor_index<Mneighbor_counter_new;neighbor_index++){
						i_neighbor=temporal_neighbor_list_M_after[neighbor_index];
						half_pair_energy=temporal_energy_list_M_after[neighbor_index];
						etot_list[i_move]=etot_list[i_move]+half_pair_energy;
						etot_list[i_neighbor]=etot_list[i_neighbor]+half_pair_energy;
					}
				}
				
				/*Addition MP*/
				if (Pneighbor_counter_new!=0){
					for (neighbor_index=0;neighbor_index<Pneighbor_counter_new;neighbor_index++){
						i_neighbor=temporal_neighbor_list_P_after[neighbor_index];
						half_pair_energy=temporal_energy_list_P_after[neighbor_index];
						etot_list[i_move]=etot_list[i_move]+half_pair_energy;
						etot_list[N_mobile_tot+i_neighbor]=etot_list[N_mobile_tot+i_neighbor]+half_pair_energy;
					}
				}
					
					
					
////////////////Phase 5: zero temporal lists////////////////
				for (neighbor_index=0;neighbor_index<Mneighbor_counter_ori;neighbor_index++){
					temporal_neighbor_list_M_before[neighbor_index]=0;
					temporal_energy_list_M_before[neighbor_index]=0;
				}
				for (neighbor_index=0;neighbor_index<Pneighbor_counter_ori;neighbor_index++){
					temporal_neighbor_list_P_before[neighbor_index]=0;
					temporal_energy_list_P_before[neighbor_index]=0;
				}			
			}
			for (neighbor_index=0;neighbor_index<Mneighbor_counter_new;neighbor_index++){
				temporal_neighbor_list_M_after[neighbor_index]=0;
				temporal_energy_list_M_after[neighbor_index]=0;
			}
			for (neighbor_index=0;neighbor_index<Pneighbor_counter_new;neighbor_index++){
				temporal_neighbor_list_P_after[neighbor_index]=0;
				temporal_energy_list_P_after[neighbor_index]=0;
			}	
		}
	}
	free(temporal_energy_list_M_after);
	free(temporal_energy_list_P_after);
	free(temporal_energy_list_M_before);
	free(temporal_energy_list_P_before);
	free(temporal_neighbor_list_M_after);
	free(temporal_neighbor_list_P_after);
	free(temporal_neighbor_list_M_before);
	free(temporal_neighbor_list_P_before);			
}
		

void Energy_Profiler(double lambda,int N_mobile_tot,int N_pinned_tot,int NA_mobile,int NA_pinned,int n,double etot_list[],double x_M[],double y_M[],double z_M[],int a_M[],double x_P[],double y_P[],double z_P[],int a_P[],int cell_has_which_Mparticle[],int cell_has_which_Pparticle[],int in_which_cell_M[],int in_which_cell_P[],int cell_has_how_much_Mparticle[],int cell_has_how_much_Pparticle[],int neighbor_search_facilitator_x[],int neighbor_search_facilitator_y[],int neighbor_search_facilitator_z[],double interaction_range){
    
    ////////////////HARD-CODED Binary Lennard-Jones liquid parameters (see, e.g., J. Chem. Phys. 144, 024501 (2016))////////////////
    double sigmaAA=(double)1.00;
    double sigmaAB=(double)0.80;
    double sigmaBB=(double)0.88;
    double eAA=(double)1.0;
    double eAB=(double)1.5;
    double eBB=(double)0.5;
    double rcutAA=(double)2.5;
    double rcutAB=(double)2.0;
    double rcutBB=(double)2.2;
    /*Quickly cook Binary Lennard-Jones liquid parameters*/
    double sigmaAA_sixth=CUBE(SQR(sigmaAA));
    double sigmaAB_sixth=CUBE(SQR(sigmaAB));
    double sigmaBB_sixth=CUBE(SQR(sigmaBB));
    double sigmaAA_twelfth=SQR(sigmaAA_sixth);
    double sigmaAB_twelfth=SQR(sigmaAB_sixth);
    double sigmaBB_twelfth=SQR(sigmaBB_sixth);
    double rcutAA_square=SQR(rcutAA);
    double rcutAB_square=SQR(rcutAB);
    double rcutBB_square=SQR(rcutBB);
    double r_sixth;
    double r_twelfth;
    double VshiftAA;
    double VshiftAB;
    double VshiftBB;
    if (rcutAA>interaction_range || rcutAB>interaction_range || rcutBB>interaction_range){printf("Interaction range chosen too small\n");exit(EXIT_FAILURE);}
    r_sixth=CUBE(rcutAA_square);
    r_twelfth=SQR(r_sixth);
    VshiftAA=-((double)2.0)*eAA*(-(sigmaAA_sixth/r_sixth)+(sigmaAA_twelfth/r_twelfth));
    r_sixth=CUBE(rcutAB_square);
    r_twelfth=SQR(r_sixth);
    VshiftAB=-((double)2.0)*eAB*(-(sigmaAB_sixth/r_sixth)+(sigmaAB_twelfth/r_twelfth));
    r_sixth=CUBE(rcutBB_square);
    r_twelfth=SQR(r_sixth);
    VshiftBB=-((double)2.0)*eBB*(-(sigmaBB_sixth/r_sixth)+(sigmaBB_twelfth/r_twelfth));
    
    double sigma_sixth_temp;
    double sigma_twelfth_temp;
    double Vshift_temp;
    double e_temp;
    double rcut_square_temp;
    
////////////////Initial assignment of variables////////////////
	int i_given;
	int i_neighbor;
    int a_given;
    int a_neighbor;
	int c_index;
	int c_x;
	int c_y;
	int c_z;
	double x_temp;
	double y_temp;
	double z_temp;
	double radial_distance_square;
	double half_pair_energy;
	
	int c_index_neighbor; /*cell index*/
	int c_x_neighbor; /*cell index in x direction*/
	int c_y_neighbor; /*cell index in y direction*/
	int c_z_neighbor; /*cell index in z direction*/
	int n_cell_max=NCELLMAX; /*n_cell_max sets the maximum number of particles in one cell allowed*/
	int neighbor_cell_index;
	int neighbor_cell_index_max=DIMPOWER(3);
	int index_within_cell;
	int n_particles_in_cell;
	
	double shrinking_modificationUU=1/(SQR(lambda)); /*When two particles are both unpinned and shrunk*/
	double shrinking_modificationUP=1/(SQR((1+lambda)/2)); /*When one particle is unpinned and shrunk*/
	
	/*For a given mobile-particle,*/
	for (i_given=1;i_given<(N_mobile_tot+1);i_given++){
		c_index=in_which_cell_M[i_given];
		c_x=c_index%n;
		c_y=((int)(floor(c_index/(double)n)))%n;
		c_z=(int)floor(c_index/(((double)n)*n));
        a_given=a_M[i_given];

		/*refresh first*/
		etot_list[i_given]=0;
		
		/*Look around neighbor cells WITHIN the box of interest*/
		for (neighbor_cell_index=0;neighbor_cell_index<neighbor_cell_index_max;neighbor_cell_index++){
			c_x_neighbor=c_x+neighbor_search_facilitator_x[neighbor_cell_index];
			c_y_neighbor=c_y+neighbor_search_facilitator_y[neighbor_cell_index];
			c_z_neighbor=c_z+neighbor_search_facilitator_z[neighbor_cell_index];
			if (c_x_neighbor>=0 && c_x_neighbor<n && c_y_neighbor>=0 && c_y_neighbor<n && c_z_neighbor>=0 && c_z_neighbor<n){
				c_index_neighbor=c_x_neighbor+c_y_neighbor*n+c_z_neighbor*(n*n);
				
				/*Pick mobile-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Mparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Mparticle[c_index_neighbor*n_cell_max+index_within_cell];
                        a_neighbor=a_M[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
						else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
						else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
						}
						
						
						/*Calculate pair distance*/
						/*NO SELF ENERGY*/
						if (i_neighbor!=i_given){
							x_temp=x_M[i_given]-x_M[i_neighbor];
							y_temp=y_M[i_given]-y_M[i_neighbor];
							z_temp=z_M[i_given]-z_M[i_neighbor];
							radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
							
							/*Possibly shrunk*/
							radial_distance_square=radial_distance_square*shrinking_modificationUU;
							
							if (radial_distance_square<rcut_square_temp){
								r_sixth=CUBE(radial_distance_square);
								r_twelfth=SQR(r_sixth);
								
								/*Calculate pair LJ energy, ``half of it"*/
								half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;
								
								/*Add it to potential energy*/
								etot_list[i_given]=etot_list[i_given]+half_pair_energy;
							}
						}					
					}
				}
				
				/*Pick pinned-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Pparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Pparticle[c_index_neighbor*n_cell_max+index_within_cell];
						a_neighbor=a_P[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                        else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                        else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                        }
						
						/*Calculate pair distance*/
						
						x_temp=x_M[i_given]-x_P[i_neighbor];
						y_temp=y_M[i_given]-y_P[i_neighbor];
						z_temp=z_M[i_given]-z_P[i_neighbor];
						radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
						
						/*Possibly shrunk*/
						radial_distance_square=radial_distance_square*shrinking_modificationUP;
						
						if (radial_distance_square<rcut_square_temp){
							r_sixth=CUBE(radial_distance_square);
							r_twelfth=SQR(r_sixth);
							
							/*Calculate pair LJ energy, ``half of it"*/
							half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;
							
							/*Add it to potential energy*/
							etot_list[i_given]=etot_list[i_given]+half_pair_energy;
						}
						
					}
				}
			}			
		}
	}
	
	
	
	/*For a given pinned-particle,*/
	for (i_given=1;i_given<(N_pinned_tot+1);i_given++){
		c_index=in_which_cell_P[i_given];
		c_x=c_index%n;
		c_y=((int)(floor(c_index/(double)n)))%n;
		c_z=(int)floor(c_index/(((double)n)*n));
        a_given=a_P[i_given];
		
		/*refresh first*/
		etot_list[N_mobile_tot+i_given]=0;
		
		/*Look around neighbor cells WITHIN the box of interest*/
		for (neighbor_cell_index=0;neighbor_cell_index<neighbor_cell_index_max;neighbor_cell_index++){
			c_x_neighbor=c_x+neighbor_search_facilitator_x[neighbor_cell_index];
			c_y_neighbor=c_y+neighbor_search_facilitator_y[neighbor_cell_index];
			c_z_neighbor=c_z+neighbor_search_facilitator_z[neighbor_cell_index];
			if (c_x_neighbor>=0 && c_x_neighbor<n && c_y_neighbor>=0 && c_y_neighbor<n && c_z_neighbor>=0 && c_z_neighbor<n){
				c_index_neighbor=c_x_neighbor+c_y_neighbor*n+c_z_neighbor*(n*n);
				
				/*Pick mobile-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Mparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Mparticle[c_index_neighbor*n_cell_max+index_within_cell];
                        a_neighbor=a_M[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                        else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                        else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                        }
						
						/*Calculate pair distance*/
						x_temp=x_P[i_given]-x_M[i_neighbor];
						y_temp=y_P[i_given]-y_M[i_neighbor];
						z_temp=z_P[i_given]-z_M[i_neighbor];
						radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
						
						/*Possibly shrunk*/
						radial_distance_square=radial_distance_square*shrinking_modificationUP;
						
						if (radial_distance_square<rcut_square_temp){
							r_sixth=CUBE(radial_distance_square);
							r_twelfth=SQR(r_sixth);
							
							/*Calculate pair LJ energy, ``half of it"*/
							half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;
							
							/*Add it to potential energy*/
							etot_list[N_mobile_tot+i_given]=etot_list[N_mobile_tot+i_given]+half_pair_energy;
						}					
					}
				}
				
				/*Pick pinned-neighbors, one by one*/
				n_particles_in_cell=cell_has_how_much_Pparticle[c_index_neighbor];
				if (n_particles_in_cell!=0){
					for (index_within_cell=0;index_within_cell<n_particles_in_cell;index_within_cell++){
						i_neighbor=cell_has_which_Pparticle[c_index_neighbor*n_cell_max+index_within_cell];
                        a_neighbor=a_P[i_neighbor];
                        
                        if (a_given==0 && a_neighbor==0){sigma_sixth_temp=sigmaAA_sixth; sigma_twelfth_temp=sigmaAA_twelfth;Vshift_temp=VshiftAA; e_temp=eAA; rcut_square_temp=rcutAA_square;}
                        else {if (a_given==1 && a_neighbor==1){sigma_sixth_temp=sigmaBB_sixth; sigma_twelfth_temp=sigmaBB_twelfth;Vshift_temp=VshiftBB; e_temp=eBB; rcut_square_temp=rcutBB_square;}
                        else {sigma_sixth_temp=sigmaAB_sixth; sigma_twelfth_temp=sigmaAB_twelfth;Vshift_temp=VshiftAB; e_temp=eAB; rcut_square_temp=rcutAB_square;}
                        }
						
						/*Calculate pair distance*/
						/*NO SELF ENERGY*/
						if (i_neighbor!=i_given){
							x_temp=x_P[i_given]-x_P[i_neighbor];
							y_temp=y_P[i_given]-y_P[i_neighbor];
							z_temp=z_P[i_given]-z_P[i_neighbor];
							radial_distance_square=SQR(x_temp)+SQR(y_temp)+SQR(z_temp);
							
							/*DO NOT SHRINK: radial_distance_square=radial_distance_square;*/
							
							if (radial_distance_square<rcut_square_temp){
								r_sixth=CUBE(radial_distance_square);
								r_twelfth=SQR(r_sixth);
								
								/*Calculate pair LJ energy, ``half of it"*/
								half_pair_energy=2*e_temp*(-(sigma_sixth_temp/r_sixth)+(sigma_twelfth_temp/r_twelfth))+Vshift_temp;
								
								/*Add it to potential energy*/
								etot_list[N_mobile_tot+i_given]=etot_list[N_mobile_tot+i_given]+half_pair_energy;
							}
						}					
					}
				}
			}			
		}
	}
}
			
			
			
void Take_One_Replica_Out(int replica_index,int N_mobile_tot,int N_pinned_tot,int n_replica,int c_index_max,double x_M[],double y_M[],double z_M[],int in_which_cell_M[],int cell_has_which_Mparticle[],int where_in_cell_list_M[],int cell_has_how_much_Mparticle[], double etot_list[],double x_M_replicated[],double y_M_replicated[],double z_M_replicated[],int in_which_cell_M_replicated[],int cell_has_which_Mparticle_replicated[],int where_in_cell_list_M_replicated[],int cell_has_how_much_Mparticle_replicated[], double etot_list_replicated[]){
	
    int i;
    int c_index; /*cell index=c_x+c_y*n+c_z*(n*n);*/
    int index_within_cell; /*indexes particles within a given cell*/
    int n_cell_max=NCELLMAX; /*n_cell_max sets the maximum number of particles in one cell allowed*/

	for(i=1;i<N_mobile_tot+1;i++){
		x_M[i]=x_M_replicated[(i)*(n_replica)+replica_index];
		y_M[i]=y_M_replicated[(i)*(n_replica)+replica_index];
		z_M[i]=z_M_replicated[(i)*(n_replica)+replica_index];
		in_which_cell_M[i]=in_which_cell_M_replicated[(i)*(n_replica)+replica_index];
		where_in_cell_list_M[i]=where_in_cell_list_M_replicated[(i)*(n_replica)+replica_index];		
	}
	for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
		etot_list[i]=etot_list_replicated[(i)*(n_replica)+replica_index];			
	}
	
	for (c_index=0;c_index<c_index_max;c_index++){
		cell_has_how_much_Mparticle[c_index]=cell_has_how_much_Mparticle_replicated[(c_index)*(n_replica)+replica_index];
		for (index_within_cell=0;index_within_cell<n_cell_max;index_within_cell++){
			cell_has_which_Mparticle[c_index*n_cell_max+index_within_cell]=cell_has_which_Mparticle_replicated[(c_index*n_cell_max+index_within_cell)*(n_replica)+replica_index];
		}	
	}
}
			
			
void Put_One_Replica_Back(int replica_index,int N_mobile_tot,int N_pinned_tot,int n_replica,int c_index_max,double x_M[],double y_M[],double z_M[],int in_which_cell_M[],int cell_has_which_Mparticle[],int where_in_cell_list_M[],int cell_has_how_much_Mparticle[], double etot_list[],double x_M_replicated[],double y_M_replicated[],double z_M_replicated[],int in_which_cell_M_replicated[],int cell_has_which_Mparticle_replicated[],int where_in_cell_list_M_replicated[],int cell_has_how_much_Mparticle_replicated[], double etot_list_replicated[]){
    
	int i;
    int c_index; /*cell index=c_x+c_y*n+c_z*(n*n);*/
    int index_within_cell; /*indexes particles within a given cell*/
	int n_cell_max=NCELLMAX; /*n_cell_max sets the maximum number of particles in one cell allowed*/
	
	for(i=1;i<N_mobile_tot+1;i++){
		x_M_replicated[(i)*(n_replica)+replica_index]=x_M[i];
		y_M_replicated[(i)*(n_replica)+replica_index]=y_M[i];
		z_M_replicated[(i)*(n_replica)+replica_index]=z_M[i];
		in_which_cell_M_replicated[(i)*(n_replica)+replica_index]=in_which_cell_M[i];
		where_in_cell_list_M_replicated[(i)*(n_replica)+replica_index]=where_in_cell_list_M[i];
	}
	for(i=1;i<N_mobile_tot+N_pinned_tot+1;i++){
		etot_list_replicated[(i)*(n_replica)+replica_index]=etot_list[i];			
	}
	
	for (c_index=0;c_index<c_index_max;c_index++){
		cell_has_how_much_Mparticle_replicated[(c_index)*(n_replica)+replica_index]=cell_has_how_much_Mparticle[c_index];
		for (index_within_cell=0;index_within_cell<n_cell_max;index_within_cell++){
			cell_has_which_Mparticle_replicated[(c_index*n_cell_max+index_within_cell)*(n_replica)+replica_index]=cell_has_which_Mparticle[c_index*n_cell_max+index_within_cell];
		}	
	}
}
