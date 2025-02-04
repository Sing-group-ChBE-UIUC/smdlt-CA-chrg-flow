#include "Initialization.h"

void Initialization(int argc, char * argv[]){

    // clustero = malloc(sizeof(char)*100);
    // sprintf(clustero,"qsub/out.txt");
    // Parse command line inputs
    ParseInput(argc,argv);
    // Read text file inputs
    readInput();
    // Initialize periodic simulation cell size
    initBox();
    // Initialize KRBC lattice
    initLattice();
    // Initialize Verlet neighbor tables
    initVerlet();
    // Allocate arrays
    allocate();
    // Initialize grid space HI and CA TEA
    initEwald();
    // Print parameters to an output file
    printOutput();
}

void ParseInput(int argc, char * argv[]){
    // Default values if command line argument not given
    // double m = 0.58942199, b = 0.155783717; // Scaling parameters for log(Rg) = m*log(N) + b
    // double m = 0.56701369, b = 0.184561671; // Normal LJ, attractive Scaling parameters for log(Rg) = m*log(N) + b
    // double b = 0.7924, m = 0.6392; // Kremer-Grest, Rg = b*N_{b}^m
    c_norm = 0.1;
    Nb = 50;
    Nc = 20;
    // tmax = 31000001;
    tmax = (int)20e6;
    num_threads = 4;
    trace = 0;
    Ncharges = 1;
    lambda_d = 5.0;
    restart = 0;
    barrier = 3.0;
    barrier2 = 3.0;
    flowrate = 0.1;
    nt_avg = 8;
    fixSeed = 0;
    dt = 0.0005;
    dt_eq = 0.001;
    iterstart = 0;
    int option = 0;
    flowtype = 1;

    if(argc < 4) printf("Warning: Using defaults for unspecified command line arguments.\n");
    while((option = getopt(argc, argv, "c:b:a:t:d:p:n:r:f:s:i:m:g:h:l:x:y:o:")) != -1){
        switch(option){
            // Concentration normalized by c*
            case 'c':
                sscanf(optarg, "%lf", &c_norm);
                break;
            // Chain length
            case 'b':
                sscanf(optarg, "%d", &Nb);
                break;
            // Number of chains
            case 'a':
                sscanf(optarg, "%d", &Nc);
                break;
            // Total accumulated strain
            // case 'e':
            //     sscanf(optarg, "%lf", &total_strain);
            //     break;
            // Time step
            case 'd':
                sscanf(optarg, "%lf", &dt);
                break;
            // Number of time steps
            case 't':
                sscanf(optarg, "%lu", &tmax_temp);
                break;
            // Trajectory number
            case 'p':
                sscanf(optarg, "%d", &trace);
                break;
            // Number of parallel threads
            case 'n':
                sscanf(optarg, "%d", &num_threads);
                break;
            // Restart condition
            case 'r':
                sscanf(optarg, "%d", &restart);
                break;
            // Strain rate edot
            case 'f':
                sscanf(optarg, "%lf", &flowrate);
                break;
            // Number of bins per KRBC period
            case 's':
                sscanf(optarg, "%d", &nt_avg);
                break;
            // Grid resolution
            case 'g':
                sscanf(optarg, "%lf", &bin_size);
                break;
            // Code to fix seed or not
            // case 'z':
            //     sscanf(optarg, "%d", &fixSeed);
            //     break;
            // Starting iteration, 0 -> FD, >0 -> HI
            case 'i':
                sscanf(optarg, "%d", &iterstart);
                break;
            // Maximum iteration minus one
            case 'm':
                sscanf(optarg, "%d", &itermax);
                break;
            //Charge Stuff//  //also add variable to option list in line 60
            case 'h':
                sscanf(optarg, "%d", &Ncharges);
                break;
            case 'l':
                sscanf(optarg, "%lf", &lambda_d);
                break;
            case 'x':
                sscanf(optarg, "%lf", &barrier);
                break;
            case 'y':
                sscanf(optarg, "%lf", &barrier2);
                break;
            case 'o':
                sscanf(optarg, "%d", &flowtype);
                break;
            ////////////////
            case '?':
                printf("Unknown option -%c.\n Execution abort.", optopt);
                exit(EXIT_FAILURE);
        }
    }
    // double m = 0.56701369, b = 0.184561671; // Normal LJ, attractive Scaling parameters for log(Rg) = m*log(N) + b
    // double con = m*log10(Nb)+b ;
    // Rg = pow(10,con);

    double b = 0.7924, m = 0.6392; // Kremer-Grest, Rg = b*N_{b}^m
    Rg = b*pow(Nb,m);
    

    N = Nb*Nc; // Total number of beads
    // printf("%lf\n",Rg);

    if(Ncharges>N){
      printf("Error - Ncharges exceeds N, exiting\n");
      exit(1);
    }
    if(restart==1){
      printf("Error - restart=1 not yet implemented (See Original Charlie Code)\n");
      exit(1);
    }
}

void readInput(){
    FILE *inputfile;
    inputfile = fopen("Input.txt","r");
    fscanf(inputfile, "spring = %s\n", &spring);
    fscanf(inputfile, "epsilon = %lf\n", &epsilon);
    fscanf(inputfile, "kappab = %lf\n", &kappab);
    fscanf(inputfile, "kappas = %lf\n", &kappas);
    fscanf(inputfile, "qmax = %lf\n", &qmax);
    //Charge Stuff//
    fscanf(inputfile, "total_strain = %lf\n", &total_strain);
    fscanf(inputfile, "fixSeed = %d\n", &fixSeed);
    fscanf(inputfile, "lambda_b = %lf\n", &lambda_b);
    fscanf(inputfile, "chargehop_type = %d\n", &chargehop_type_temp);
    fscanf(inputfile, "MStep = %d\n", &MStep);
    fscanf(inputfile, "MSDstart = %d\n", &MSDstart);
    /////////////////
    fclose(inputfile);
}

void initBox(){
    int itertemp,tempthreads;
    unsigned long ttemp;
    // Set N, find box length for c = c_star*c_norm
    c_star = Nb/(M_PI*4/3*Rg*Rg*Rg); // Overlap concentration
    box_volume = N/(c_norm*c_star); // Box volume at the normalized concentration
    L = pow(box_volume,(double)1/3); // For cubic lattice, L = V**(1/3)
    box_side = L/2; // Box centered at origin, so sides at +- L/2

    outp = malloc(sizeof(char)*100);
    sprintf(outp, "txt/P%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace);
    // printf("%s\n",outp);

    // if(restart==0 && iterstart==0){
    //     outputfile = fopen(outp,"r");
    //     if(outputfile){
    //         printf("Error: Trying to write over existing simulation\n");
    //         printf("Try setting restart = 2 or check input\n");
    //         exit(1);
    //     }
    // }
    if(restart==1){
        outputfile = fopen(outp, "r");
        if(!outputfile){
            printf("Error: trying to restart a simulation with no existing files\n");
            printf("Try checking input or change to restart = 0\n");
            exit(1);
        }
        fclose(outputfile);
    }
    // if(restart==2){
    //     outputfile = fopen(outp, "r");
    //     if(!outputfile){
    //         printf("Error: trying to write over a simulation with no existing files\n");
    //         printf("Try checking input or change to restart = 0\n");
    //         exit(1);
    //     }
    //     fclose(outputfile);
    // }
}

void initLattice(){
    int i,j,treal;
    double theta;
    tp = (int)(EP/flowrate/dt); //0.2857 found to make sure angle>=45 deg
    // tp = (int)(0.2857/flowrate/dt); //0.2857 found to make sure angle>=45 deg
    // printf("%lf %lf %lf %d\n",EP,flowrate,dt,tp);
    double xmax,ymax,zmax;
    // Init un-rotated box
    theta = 90.0-(45.0+THETA0);
    theta *= M_PI/180;
    point[0][0] = -sqrt(2)/2*L*cos(theta); point[0][1] = sqrt(2)/2*L*sin(theta); //x,y component of point 0 
    point[1][0] = -sqrt(2)/2*L*sin(theta); point[1][1] = -sqrt(2)/2*L*cos(theta); //x,y component of point 1
    point[2][0] = sqrt(2)/2*L*cos(theta); point[2][1] = -sqrt(2)/2*L*sin(theta); //x,y component of point 2
    // printf("p[0][0]: %lf\n",point[0][0]);
    // printf("p[1][0]: %lf\n",point[1][0]);
    // printf("p[2][0]: %lf\n",point[2][0]);
    // printf("p[0][1]: %lf\n",point[0][1]);
    // printf("p[1][1]: %lf\n",point[1][1]);
    // printf("p[2][1]: %lf\n",point[2][1]);
    for(i = 0; i<3; ++i){
        for(j = 0; j<2; ++j){
            point0[i][j] = point[i][j];
        }
    }
    // ymax = -point0[1][1];
	// printf("ymax: %lf\n",ymax);
	//printf("ymax_change: %lf\n",L/2);
	ymax = L/2;
    // Find max width of un-rotated box
    treal = tp;
    // for(i = 0; i<3; ++i){
    //     point[i][0] = point0[i][0]*exp(flowrate*treal*dt); point[i][1] = point0[i][1]*exp(-flowrate*treal*dt);
    // }
    point[0][0] = point0[0][0] + point[0][1]*flowrate*treal*dt ; 
    point[1][0] = point0[1][0] + point[1][1]*flowrate*treal*dt ; 
    point[2][0] = point0[2][0] + point[2][1]*flowrate*treal*dt ; 

    //L1[0] = 0; L1[1] = L; L2[0] = L; L2[1] = 0;
    // xmax = point[2][0];
	// printf("xmax: %lf\n",xmax);
	//printf("xmax_change: %lf\n",L/2*flowrate*treal*dt);
	xmax = L/2*flowrate*treal*dt ; 
    zmax = L/2;
    // Return to undeformed box to find max height in un-rotated frame
    theta = 90.0-(45.0+THETA0);
    theta *= M_PI/180;
    point[0][0] = -sqrt(2)/2*L*cos(theta); point[0][1] = sqrt(2)/2*L*sin(theta);
    point[1][0] = -sqrt(2)/2*L*sin(theta); point[1][1] = -sqrt(2)/2*L*cos(theta);
    point[2][0] = sqrt(2)/2*L*cos(theta); point[2][1] = -sqrt(2)/2*L*sin(theta);
    for(i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
        //printf("L1[%d]=%lf L2[%d]=%lf\n",i,L1[i],i,L2[i]);
    }

    // printf("%lf %lf\n",xmax,ymax);
    // bin_size_x = 2.0;
    // bin_size_y = 2.0;
    // bin_size_z = 2.0;
    // bin_size_x = 1.0;
    // bin_size_y = 1.0;
    // bin_size_z = 1.0;
    // bin_size_x = 0.5;
    // bin_size_y = 0.5;
    // bin_size_z = 0.5;
    bin_size_x = bin_size;
    bin_size_y = bin_size;
    bin_size_z = bin_size;
    num_bins_x = 2*ceil(xmax/bin_size_x) + 1;
    bin_size_x = xmax/((num_bins_x - 1)/2);
    nhx = (num_bins_x - 1)/2;
    num_bins_y = 2*ceil(ymax/bin_size_y) + 1;
    bin_size_y = ymax/((num_bins_y - 1)/2);
    nhy = (num_bins_y - 1)/2;
    num_bins_z = 2*ceil(zmax/bin_size_z) + 1;
    bin_size_z = zmax/((num_bins_z - 1)/2);
    nhz = (num_bins_z - 1)/2;
    // printf("%lf %d %lf %d\n",zmax,num_bins_z,bin_size_z,nhz);
    // exit(1);
    num_bins = num_bins_x*num_bins_y*num_bins_z;
    //printf("box length %lf\n",L);
    printf("bin_size_x - %lf num_bins_x - %d bin_size_y - %lf num_bins - %d nhx - %d nhy - %d\n",bin_size_x,num_bins_x,bin_size_y,num_bins,nhx,nhy);
}

void initVerlet(){
    int i;
    sigma = 2.0; // 1 particle diameter = 2 radii, HS overlap, get LJ cutoff from this
    rc = 2.5*sigma; // LJ cutoff range (assumption: goes to zero outside this range)
    r2cut = rc*rc;
    rv =  3.0*sigma; // Verlet list radius. Particles inside this radius are added to the varlet list
    rv2 = rv*rv;
    rnew = 0.5*(rv-rc);
    nlist = calloc(N, sizeof(long int)); // Number of neighbors in neighbor list for each bead
    // List for neighbors in format list(bead,neighbor indices)
    // Size NxN because there are N beads requiring neighbor lists and a maximum of N possible neighbors for each bead
    list = calloc(N*N, sizeof(long int));
    // Displacement from the bead position the last time neighbor tables were updated
    dxv = calloc(N,sizeof(double));
    dyv = calloc(N,sizeof(double));
    dzv = calloc(N,sizeof(double));
}

void allocate(){
    int i,j,k,l,periods;
    printperiod = 1000;
    printprops = 1000;
    if(printprops<1000){
      printf("Warning - printprops is below 1000\n");
    }
    spp = printperiod*dt; // Samples per period
    xyz = malloc(sizeof(char)*100);
    map = malloc(sizeof(char)*100);
    rgt = malloc(sizeof(char)*100);
    mmatrix = malloc(sizeof(char)*100);
    decomp = malloc(sizeof(char)*100);
    beta = malloc(sizeof(char)*100);
    ext = malloc(sizeof(char)*100);
    grid = malloc(sizeof(char)*100);
    visc = malloc(sizeof(char)*100);
    tim = malloc(sizeof(char)*100);
    sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(map, "xyz/Mapped%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(mmatrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(grid, "mm/G%d_%d_%.3lf_%.5lf_%d.txt", Nb, Nc, c_norm,flowrate,nt_avg);
    sprintf(decomp, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(beta, "decomp/B%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(ext, "prop/E%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(visc, "prop/VF%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(tim, "prop/T%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    //Charge Stuff//
    sprintf(rgt, "prop/RG%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    str2 = malloc(sizeof(char)*300);
    outp1 = malloc(sizeof(char)*300);
    sprintf(str2,"prop/CH%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    sprintf(outp1,"prop/MSD%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,iterstart);
    Charge  = calloc(N, sizeof(int));
    Charge_track = calloc(Ncharges, sizeof(int));
    Charge_old = calloc(Ncharges, sizeof(int));
    Charge_indices = calloc(Ncharges, sizeof(int));
    check_ch = calloc(3*Ncharges, sizeof(int));
    i_values = calloc(N, sizeof(int));
    dxt = calloc(Ncharges*printprops, sizeof(long double));
    dyt = calloc(Ncharges*printprops, sizeof(long double));
    dzt = calloc(Ncharges*printprops, sizeof(long double));
    rbi = calloc(3, sizeof(double));
    rbic = calloc(3*Ncharges, sizeof(double));
    ////////////////

    rx = calloc(N, sizeof(double));
    ry = calloc(N, sizeof(double));
    rz = calloc(N, sizeof(double));
    sx = calloc(N, sizeof(double));
    sy = calloc(N, sizeof(double));
    sz = calloc(N, sizeof(double));
    rb = calloc(3*N, sizeof(double));
    dxu = calloc(N, sizeof(double));
    dyu = calloc(N, sizeof(double));
    dzu = calloc(N, sizeof(double));
    drb = calloc(3*N, sizeof(double));
    px = calloc(N, sizeof(double));
    py = calloc(N, sizeof(double));
    pz = calloc(N, sizeof(double));
    rg = calloc(Nc,sizeof(double));
    comx = calloc(Nc,sizeof(double));
    comy = calloc(Nc,sizeof(double));
    comz = calloc(Nc,sizeof(double));
    comx_r = calloc(Nc,sizeof(double));
    comy_r = calloc(Nc,sizeof(double));
    comz_r = calloc(Nc,sizeof(double));
    reex = calloc(Nc,sizeof(double));
    reey = calloc(Nc,sizeof(double));
    reez = calloc(Nc,sizeof(double));
    rext = calloc(Nc,sizeof(double));
    fx = calloc(N, sizeof(double));
    fy = calloc(N, sizeof(double));
    fz = calloc(N, sizeof(double));
    f = calloc(3*N, sizeof(double));
    R = calloc(3*N, sizeof(double));
    D = calloc(nt_avg,sizeof(double));
    for(i=0;i<nt_avg;i++){
        D[i] = calloc(9*num_bins,sizeof(double));
    }
    M = calloc(9*N*N,sizeof(double));
    ebins = 20; // 20 bins for extension, increments of 1.0/20 = 0.05
    ext_interval = 1.0/(double)ebins;
    sbins = 20;
    strain_ss = 10.0;
    if(strain_ss > total_strain){
        strain_ss = total_strain;
    }
    strain_interval = strain_ss/sbins;
    samp_per_inter = 20;
    strain_sample = (int)floor(strain_interval/(samp_per_inter*flowrate*dt));
    //printf("strain_ss - %lf strain_interval - %lf samp_per_inter - %d samp_rate - %d\n",strain_ss,strain_interval,samp_per_inter,strain_sample);
    //tmax = total_strain/(dt*flowrate);
    //printf("flow simulation time steps - %lu\n",tmax);
    // exit(1);
    Ct = calloc(3*N,sizeof(double));
    C = calloc(sbins+1,sizeof(double));
    Crun = calloc(sbins+1,sizeof(double));
    for(i=0;i<sbins+1;++i){
        C[i] = calloc(3*Nb,sizeof(double));
        Crun[i] = calloc(3*Nb,sizeof(double));
    }
    tauxx = 0.0;
    tauyy = 0.0;
    tauxy = 0.0;
    rowsum = calloc(3*N,sizeof(double));
    count = calloc(num_bins,sizeof(int));
    start_ij = calloc(N*N,sizeof(int));
    bp = calloc(sbins+1,sizeof(double));
    brun = calloc(sbins+1,sizeof(double));
    gwcount = calloc(sbins+1,sizeof(int));
    p = sqrt(2.0/dt);
    pc = sqrt(2.0*dt);
    pc_eq = sqrt(2.0*dt_eq);
    srpy = 100;
    idum = malloc(sizeof(long)); // = initRan(); // for testing, set seed
    if(fixSeed==0){
        *idum = initRan();
    }
    else if(fixSeed==1){
        *idum = 830062920;
    }
    else{
        printf("Exit code: fixSeed must by specified\n");
        exit(1);
    }
    if(*idum==2 || fixSeed==1){
        printf("Warning! Seed is fixed!\n");
    }
    //qmax = 2.0*1.5;
    qmax2 = qmax*qmax;
    H = 30.0*epsilon/4.0;
    tau_FD = 2.0*554.108;
    equil_start = (int)(10.0*tau_FD/dt_eq);
    periods = ceil(equil_start/tp);
    // equil_start = periods*tp;
    equil_start = 0;
    // equil_start = tp;
    tmax += equil_start;
}

void initEwald(){
    int i,j,k,l,m,n,o,kk,kx,ky,kz,a,b,start,atemp,btemp,ctemp,dtemp,nt;
    double rkx,rky,rkz,rkk,m2,kcoeff;
    unsigned long ttemp;
    eps = 0.0;
    // If FD preaveraging run was not finished, load in running average TEA parameters
    if(restart==1 &&iterstart==0){
        DCfile = fopen(decomp,"r");
        if(DCfile){
            DCfile = fopen(decomp,"r");
            fscanf(DCfile,"%lu\n",&ttemp);
            for(i=0;i<sbins+1;++i){
                fscanf(DCfile,"%d %d %le\n",&atemp,&gwcount[i],&brun[i]);
                for(j=0;j<Nb;j++){
                    fscanf(DCfile,"%d %le %le %le\n",&atemp,&Crun[i][3*j],&Crun[i][3*j+1],&Crun[i][3*j+2]);
                }
                brun[i] *= gwcount[i];
                for(j=0;j<3*N;j++){
                    Crun[i][j] *= (Nc*gwcount[i]);
                }
            }
            fclose(DCfile);
        }
        else{
            fclose(DCfile);
            printf("Error - couldn't find decomp file %s to restart with, exiting\n",decomp);
            exit(1);
        }
    }
    // Check the values if something is going wrong during time integration
    // for(i=0;i<sbins+1;++i){
    //  printf("%d %d %lf\n",i,gwcount[i],brun[i]);
    //  for(j=0;j<Nb;++j){
    //      printf("%d %lf %lf %lf\n",j,Crun[i][3*j],Crun[i][3*j+1],Crun[i][3*j+2]);
    //  }
    // }
    // exit(1);
    // If starting an HI iteration, load in look table (grid space) HI
    if(iterstart>0){
        printf("%s\n",grid);
        gridfile = fopen(grid,"r");
        for(nt=0;nt<nt_avg;nt++){
            // printf("nt - %d\n",nt);
            // fscanf(gridfile,"nt %d\n",&atemp);
            for(i=0;i<num_bins_x;i++){
                for(j=0;j<num_bins_y;j++){
                    for(k=0;k<num_bins_z;k++){
                        start = i*num_bins_y*num_bins_z + j*num_bins_z + k;
                        fscanf(gridfile,"%d %d %d %d\n",&atemp,&btemp,&ctemp,&dtemp);
                        for(l=0;l<3;l++){
                            for(m=0;m<3;m++){
                                // fscanf(MMfile,"%lf ",&D[nt][9*start + 3*l + m]);
                                fscanf(gridfile,"%le ",&D[nt][9*start + 3*l + m]);
                            }
                            fscanf(gridfile,"\n");
                        }
                        // printf("bin_x - %d bin_y - %d bin_z - %d D - %lf\n",i,j,k,D[nt][9*start]);
                    }
                }
            }
        }
        fclose(gridfile);
        // printf("Ds[0] - %lf D[origin] - %lf D[0] - %lf\n",Ds[0],D[1][9*(num_bins_y*num_bins_z*nhx + num_bins_z*nhy + nhz)],D[1][0]);
        sprintf(decomp, "decomp/DC%d_%d_%.3lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, flowrate, trace, iterstart-1);
        DCfile = fopen(decomp,"r");
        if(!DCfile){
          printf("Error - couldn't find decomp file %s to restart with, exiting\n",decomp);
          exit(1);
        }
        fscanf(DCfile,"%lu\n",&ttemp);
        for(i=0;i<sbins+1;++i){
            fscanf(DCfile,"%d %d %le\n",&atemp,&btemp,&bp[i]);
            // printf("%d %d %lf\n",i,atemp,bp[i]);
            for(j=0;j<Nb;j++){
                fscanf(DCfile,"%d %le %le %le\n",&atemp,&C[i][3*j],&C[i][3*j+1],&C[i][3*j+2]);
                // printf("%d %lf %lf %lf\n",atemp,C[i][3*j],C[i][3*j+1],C[i][3*j+2]);
            }
        }
        fclose(DCfile);
        // printf("%lf %lf %lf\n",D[0][9*(nhx*num_bins_y*num_bins_z + nhy*num_bins_z + nhz)],bp[0],C[0][0]);
        // exit(1);
        sprintf(decomp, "decomp/DC%d_%d_%.3lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, flowrate, trace, iterstart);
    }
}

void printOutput(){
    omp_set_num_threads(num_threads);

    outputfile = fopen(outp, "w");
    fprintf(outputfile, "epsilon = %lf\n", epsilon);
    fprintf(outputfile,"spring = %s\n",spring);
    fprintf(outputfile, "kappab = %lf\n", kappab);
    fprintf(outputfile,"kappas = %lf\n",kappas);
    fprintf(outputfile, "H = %lf\n", H);
    fprintf(outputfile, "qmax = %lf\n", qmax);
    fprintf(outputfile, "c_norm = %lf\n", c_norm);
    
    fprintf(outputfile, "N = %d\n", N);
    fprintf(outputfile, "Nb = %d\n", Nb);
    fprintf(outputfile, "Nc = %d\n", Nc);
    ///charge stuff
    fprintf(outputfile,"conc = %lf\n", conc);
    fprintf(outputfile,"Ncharges = %d\n",Ncharges);
    fprintf(outputfile,"lambda_d = %lf\n",lambda_d);
    fprintf(outputfile,"lambda_b = %lf\n",lambda_b);
    fprintf(outputfile,"MStep = %d\n",MStep);
    fprintf(outputfile,"chargehop_type = %d\n",chargehop_type_temp);
    fprintf(outputfile,"barrier = %lf\n",barrier);
    fprintf(outputfile,"barrier2 = %lf\n",barrier2);
    fprintf(outputfile,"flowtype = %d\n",flowtype);
    fprintf(outputfile,"MSDstart = %d\n",MSDstart);
    //
    fprintf(outputfile, "edot = %lf\n", flowrate);
    if(flowrate==0){fprintf(outputfile, "tp = %d\n", 0);}
    else{fprintf(outputfile, "tp = %d\n", tp);}
    fprintf(outputfile, "Rg = %lf\n", Rg);
    fprintf(outputfile, "dt = %lf\n", dt);
    fprintf(outputfile, "tmax = %lu\n", tmax_temp);
    fprintf(outputfile, "equil_start = %lu\n", equil_start);
    fprintf(outputfile, "L = %lf\n", L);
    fprintf(outputfile, "rv = %lf\n", rv);
    fprintf(outputfile, "num_threads = %d\n",num_threads);
    fprintf(outputfile, "trace = %d\n", trace);
    fprintf(outputfile, "restart = %d\n",restart);
    fprintf(outputfile, "bin_size_x = %lf\n",bin_size_x);
    fprintf(outputfile, "bin_size_y = %lf\n",bin_size_y);
    fprintf(outputfile, "bin_size_z = %lf\n",bin_size_z);
    fprintf(outputfile, "num_bins_x = %d\n",num_bins_x);
    fprintf(outputfile, "num_bins_y = %d\n",num_bins_y);
    fprintf(outputfile, "num_bins_z = %d\n",num_bins_z);
    fprintf(outputfile, "num_bins = %d\n",num_bins);
    fprintf(outputfile, "strain_ss = %lf\n",strain_ss);
    fprintf(outputfile, "sbins = %d\n",sbins);
    fprintf(outputfile, "total_strain = %lf\n",total_strain);
    fprintf(outputfile, "iterstart = %d\n",iterstart);
    fprintf(outputfile, "itermax = %d\n",itermax);
    fprintf(outputfile, "\n");
    fprintf(outputfile, "\n");
    fprintf(outputfile, "SEED %ld\n", *idum);
    fclose(outputfile);
}

float ran1(long *idum){
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum <= 0)
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for(j=NTAB+7;j>=0;--j)
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if(*idum<0) *idum+=IM1;
			if(j<NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if(*idum<0) *idum += IM1;
	k=idum2/IQ2;
	if(*idum<0) idum2+= IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if(iy<1) iy += IMM1;
	if((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

long initRan(){
    //time_t seconds;
    //time(&seconds);
    //return -1*(unsigned long)(seconds/12345); This is bad.  :(

    //This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
    //http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000; //careful here.  Another 0 might break the ran1 (long long instead of just long)
}
