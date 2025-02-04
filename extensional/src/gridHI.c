#include "gridHI.h"

int main(int argc, char * argv[]){
	int a,b,i,j,k,nt,kx,ky,kz,count,nx,ny,nz;
	int ylower,yupper,zlower,zupper,start;
	double rkx,rky,rkz,rkxy,rkk,m2,dx,dy,dz,term1,term2,C1,C2,C3,C4,Dt[9],M1[9],coskr,rr,r;
	Vector3D_t dr;
	ParseInput(argc,argv);
	// readInput();
	initBox();
	initLattice();
	//printf("tp: %d\n",tp);
	allocate();
	for(nt=0;nt<nt_avg;nt++){
		t = nt*tp/nt_avg;
		// t1 = nt*tp/nt_avg;
		// t2 = (nt+1)*tp/nt_avg;
		// t = (t1+t2)/2;
		updateLattice(t);
		// printf("t - %lu nt - %d tp - %d treal - %d\n",t,nt,tp,treal);
		if(treal <= tp/2){
			// rce = 1.0*1.5*L2yp;
			rce = L2yp;
			nmax = 1;
		}
		else{
			// rce = 2.0*1.5*L2yp;
			rce = 2.0*L2yp;
			nmax = 2;
		}
		// M_HI = 10.0;
		M_HI = 6.0;
		// nmax = ceil(rce/(1.5*L2yp));
		alpha = M_HI/rce;
		rkcut = 2.0*M_HI*M_HI/rce;
		rkkcut = rkcut*rkcut;
		alpha2 = alpha*alpha; alpha3 = alpha2*alpha; alpha4 = alpha2*alpha2; alpha5 = alpha4*alpha; alpha7 = alpha4*alpha2*alpha;
		selfconst = 1-6/rtpi*alpha + 40/(3*rtpi)*alpha3;
		kxmax = ceil(rkcut*detL/(M_PI*2.0*L*L2[1]));
		kymax = ceil(rkcut*detL/(M_PI*2.0*L*L1[0]));
		kzmax = ceil(rkcut*L/(M_PI*2.0));
		int *kyupper = calloc(kxmax+1,sizeof(int));
		int *kylower = calloc(kxmax+1,sizeof(int));
		int *kzupper = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
		int *kzlower = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
		printf("treal - %d rce - %lf L1[0] = %lf L2[0] = %lf L1[1] = %lf L2[1] = %lf\n",treal,rce,L1[0],L2[0],L1[1],L2[1]);
		kcoeff = M_PI*2.0/detL;
		for(kx=0;kx<kxmax+1;kx++){
			ylower = 0; yupper = 0;
			for(ky=-kymax;ky<kymax+1;ky++){
				rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
				rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
				rkxy = rkx*rkx + rky*rky;
				zlower = 0; zupper = 0;
				for(kz=-kzmax;kz<kzmax+1;kz++){
					rkz = kz*M_PI*2.0/L;
					rkk = rkxy + rkz*rkz;
					// printf("kx - %d ky - %d kz - %d rkx - %lf rky - %lf rkz - %lf rkxy - %lf rkk - %lf\n",kx,ky,kz,rkx,rky,rkz,rkxy,rkk);
					if(rkk < rkkcut && zlower==0){
						kzlower[(2*kymax+1)*kx + ky + kymax] = kz;
						zlower = 1;
					// break;
					}
					if(rkk > rkkcut && zlower==1 && zupper==0){
						kzupper[(2*kymax+1)*kx + ky + kymax] = kz - 1;
						zupper = 1;
					}
				}
				if(rkxy < rkkcut && ylower==0){
					kylower[kx] = ky;
					ylower = 1;
					// printf("ylower %d %d %lf %d\n",kx,ky,rkxy,kylower[kx]);
				}
				if(rkxy > rkkcut && ylower==1 && yupper==0){
					// printf("yupper %d %d %lf\n",kx,ky,rkxy);
					kyupper[kx] = ky - 1;
					yupper = 1;
				}
			}
		}
		kx = 0;
		while(kyupper[kx]!=-kymax && kx <= (int)(rkcut*detL/(M_PI*2.0*L*L2yp))){
			// printf("%d %d\n",kx,kyupper[kx]);
			kxmax = kx;
			kx++;
		}
		count = 0;
		kxmax++;
		for(kx=0;kx<kxmax+1;kx++){
			for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
				for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
					count++;
				}
			}
		}
		// Compute M2 coefficients for each wave vectors
		double *M2 = calloc(9*count,sizeof(double));
		count = 0;
		for(kx = 0; kx < kxmax+1; kx++){
			for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
				rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
				rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
				rkxy = rkx*rkx + rky*rky;
				for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
					rkz = kz*M_PI*2.0/L;
					rkk = rkxy + rkz*rkz;
					m2 = 2.0*(1.0-rkk/3.0)*(1.0+rkk/(4.0*alpha2)+rkk*rkk/(8.0*alpha4))*6.0*M_PI/rkk*exp(-rkk/(4.0*alpha2));
					start = 9*count;
					M2[start] = m2*(1-rkx*rkx/rkk)/box_volume;
					M2[start+1] = m2*(-rkx*rky/rkk)/box_volume;
					M2[start+2] = m2*(-rkx*rkz/rkk)/box_volume;
					M2[start+3] = m2*(-rky*rkx/rkk)/box_volume;
					M2[start+4] = m2*(1-rky*rky/rkk)/box_volume;
					M2[start+5] = m2*(-rky*rkz/rkk)/box_volume;
					M2[start+6] = m2*(-rkz*rkx/rkk)/box_volume;
					M2[start+7] = m2*(-rkz*rky/rkk)/box_volume;
					M2[start+8] = m2*(1-rkz*rkz/rkk)/box_volume;
					count++;
				}
			}
		}
		// Divide by 2 for x = 0 to account for symmetry in x. Set 0,0,0 equal to 0
		count = 0;
		for(ky = kylower[0]; ky < kyupper[0] + 1; ky++){
			for(kz = kzlower[(2*kymax + 1)*0 + ky + kymax]; kz < kzupper[(2*kymax + 1)*0 + ky + kymax] + 1; kz++){
				start = 9*count;
				for(i=0;i<9;i++){
					M2[start + i] /= 2.0;
				}
				if(ky==0 && kz==0){
					for(i=0;i<9;i++){
						M2[start+i] = 0.0;
					}
				}
				count++;
			}
		}
		for(i=0;i<num_bins;i++){
			if(i!=origin){
				dr = invbin(i);
			}
			else{
				dr.x = 0.6;
				dr.y = 0.6;
				dr.z = 0.6;
			}
			// printf("bin - %d dx - %lf dy - %lf dz - %lf\n",i,dr.x,dr.y,dr.z);

			//Compute M1 coefficients for each displacement
			memset(M1,0.0,9*sizeof(double));
			memset(Dt,0.0,9*sizeof(double));
			for(nx=-nmax;nx<nmax+1;nx++){
				for(ny=-nmax;ny<nmax+1;ny++){
					for(nz=-nmax;nz<nmax+1;nz++){
						dx = dr.x + nx*L1[0] + ny*L2[0];
						dy = dr.y + nx*L1[1] + ny*L2[1];
						dz = dr.z + nz*L;
						rr = dx*dx + dy*dy + dz*dz;
						r = sqrt(rr);
						if(r<rce){
							if(r>2.0){
								C1 = 0.75/r+0.5/(r*rr);
								C2 = 4.0*alpha7*rr*rr + 3.0*alpha3*rr - 20.0*alpha5*rr - 4.5*alpha + 14.0*alpha3 + alpha/rr;
								C3 = 0.75/r-1.5/(r*rr);
								C4 = -4.0*alpha7*rr*rr - 3.0*alpha3*rr + 16.0*alpha5*rr + 1.5*alpha - 2.0*alpha3 - 3.0*alpha/rr;
							}
							else{
								C1 = 1.0 - 0.28125*r;
								C2 = 64.0*alpha7 + 12.0*alpha3 - 80.0*alpha5 - 4.5*alpha + 14.0*alpha3 + 0.25*alpha;
								C3 = 0.09375*r;
								C4 = -64.0*alpha7 - 12.0*alpha3 + 64.0*alpha5 + 1.5*alpha - 2.0*alpha3 - 0.75*alpha;
							}
							term1 = C1*erfc(alpha*r) + C2*exp(-alpha2*rr)/rtpi;
							term2 = C3*erfc(alpha*r) + C4*exp(-alpha2*rr)/rtpi;
							// if(i==27449){
							// 	if(nx==0 && ny==0 && nz==0){
							// 		printf("%lf %lf %lf %lf %lf %lf\n",C1,C2,C3,C4,term1,term2);
							// 	}
							// }
							M1[0] = term1 + term2*dx*dx/rr;
							M1[1] = term2*dx*dy/rr;
							M1[2] = term2*dx*dz/rr;
							M1[3] = M1[1];
							M1[4] = term1 + term2*dy*dy/rr;
							M1[5] = term2*dy*dz/rr;
							M1[6] = M1[2];
							M1[7] = M1[5];
							M1[8] = term1 + term2*dz*dz/rr;
							Dt[0] += M1[0];
							Dt[1] += M1[1];
							Dt[2] += M1[2];
							Dt[3] += M1[3];
							Dt[4] += M1[4];
							Dt[5] += M1[5];
							Dt[6] += M1[6];
							Dt[7] += M1[7];
							Dt[8] += M1[8];
						}
					}
				}
			}
			count = 0;
			for(kx = 0; kx < kxmax + 1; kx++){
				for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
					rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
					rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
					rkxy = rkx*rkx + rky*rky;
					for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
						rkz = kz*M_PI*2.0/L;
						rkk = rkxy + rkz*rkz;
						start = 9*count;
						// coskr = cos(rkx*dx + rky*dy + rkz*dz);
						coskr = cos(rkx*dr.x + rky*dr.y + rkz*dr.z);
						Dt[0] += M2[start]*coskr;
						Dt[1] += M2[start+1]*coskr;
						Dt[2] += M2[start+2]*coskr;
						Dt[3] += M2[start+3]*coskr;
						Dt[4] += M2[start+4]*coskr;
						Dt[5] += M2[start+5]*coskr;
						Dt[6] += M2[start+6]*coskr;
						Dt[7] += M2[start+7]*coskr;
						Dt[8] += M2[start+8]*coskr;
						count++;
					}
				}
			}
			// if(i==27449){
			// 	printf("bin - %d dx - %lf dy - %lf dz - %lf Dt[8] - %lf\n",i,dr.x,dr.y,dr.z,Dt[8]);
			// 	exit(1);
			// }
			// printf("Dt[0] - %.12lf\n",Dt[0]);
			for(a=0;a<9;a++){
				D[nt][9*i + a] = Dt[a];
			}
		}
		free(kyupper);
		free(kylower);
		free(kzupper);
		free(kzlower);
		free(M2);
	}
	mmfile = fopen(mm,"w");
	for(nt=0;nt<nt_avg;nt++){
		for(i=0;i<num_bins_x;i++){
			for(j=0;j<num_bins_y;j++){
				for(k=0;k<num_bins_z;k++){
					start = num_bins_y*num_bins_z*i + num_bins_z*j + k;
					fprintf(mmfile,"%d %d %d 1\n",i,j,k);
					for(a=0;a<3;a++){
						for(b=0;b<3;b++){
							fprintf(mmfile,"%.6e ",D[nt][9*start + 3*a + b]);
						}
						fprintf(mmfile,"\n");
					}
				}
			}
		}
	}
	// for(nt=0;nt<nt_avg;nt++){
	// 	for(i=0;i<num_bins;i++){
	// 		// start = num_bins_y*num_bins_z*i + num_bins_z*j + k;
	// 		fprintf(mmfile,"%d\n",i);
	// 		for(a=0;a<3;a++){
	// 			for(b=0;b<3;b++){
	// 				fprintf(mmfile,"%.6e ",D[nt][9*i + 3*a + b]);
	// 			}
	// 			fprintf(mmfile,"\n");
	// 		}
	// 	}
	// }
	fclose(mmfile);
	return 0;
}
void ParseInput(int argc, char * argv[]){
    // Default values if command line argument not given
    // double m = 0.58942199, b = 0.155783717; // Scaling parameters for log(Rg) = m*log(N) + b
    // double m = 0.56701369, b = 0.184561671; // Normal LJ, attractive Scaling parameters for log(Rg) = m*log(N) + b
    double b = 0.7924, m = 0.6392; // Kremer-Grest, Rg = b*N_{b}^m
    c_norm = 1.0;
    Nb = 50;
    Nc = 20;
    num_threads = 1;
    trace = 0;
    restart = 0;
    flowrate = 0.001;
	dt = 0.001;
	nt_avg = 8;

    int option = 0;

    if(argc < 4) printf("Warning: Using defaults for unspecified command line arguments.\n");
    while((option = getopt(argc, argv, "c:b:a:t:p:n:f:s:g:")) != -1){
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
            // Trajectory number
            case 'p':
                sscanf(optarg, "%d", &trace);
                break;
            // Number of parallel threads
            case 'n':
                sscanf(optarg, "%d", &num_threads);
                break;
            // Restart condition
            // case 'r':
            //     sscanf(optarg, "%d", &restart);
            //     break;
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
            case '?':
                printf("Unknown option -%c.\n Execution abort.", optopt);
                exit(EXIT_FAILURE);
        }
    }

    printf("bin_size %lf\n",bin_size);

    N = Nb*Nc; // Total number of beads
    // Rg = pow(10,m*log10(Nb)+b);
    Rg = b*pow(Nb,m);
    // Rg = 12.705374;
    // printf("%lf\n",Rg);
}
// void readInput(){
//     FILE *inputfile;
//     inputfile = fopen("Input.txt","r");
//     fscanf(inputfile, "epsilon = %lf\n", &epsilon);
//     fscanf(inputfile, "kappas = %lf\n", &kappas);
//     fscanf(inputfile, "kappab = %lf\n", &kappab);
//     // fscanf(inputfile, "Rg = %lf\n", &Rg); // Mean Squared Radius of Gyration
//     fscanf(inputfile, "dt = %lf\n", &dt);
//     fscanf(inputfile, "itermax = %d\n", &itermax);
//     fscanf(inputfile, "iterstart = %d\n", &iterstart);
//     fclose(inputfile);
// }
void initBox(){
    // Set N, find box length for c = c_star*c_norm
    c_star = Nb/(M_PI*4/3*Rg*Rg*Rg); // Overlap concentration
    // c_star = Nb/(4/3*M_PI*Rg*Rg*Rg); // *** THIS IS WRONG!!! *** For some reason, 4/3* ignored in this form, use M_PI*4/3
    box_volume = N/(c_norm*c_star); // Box volume at the normalized concentration
    L = pow(box_volume,(double)1/3); // For cubic lattice, L = V**(1/3)
    box_side = L/2; // Box centered at origin, so sides at +- L/2
}
void initLattice(){
    int i,j,treal;
    tp = (int)(EP/flowrate/dt);
    printf("tp: %d\n",tp);
    double xmax,ymax,zmax;
    // Init un-rotated box
    theta = 90.0-(45.0+THETA0);
    theta *= M_PI/180;
    point[0][0] = -sqrt(2)/2*L*cos(theta); point[0][1] = sqrt(2)/2*L*sin(theta);
    point[1][0] = -sqrt(2)/2*L*sin(theta); point[1][1] = -sqrt(2)/2*L*cos(theta);
    point[2][0] = sqrt(2)/2*L*cos(theta); point[2][1] = -sqrt(2)/2*L*sin(theta);
    for(i = 0; i<3; ++i){
        for(j = 0; j<2; ++j){
            point0[i][j] = point[i][j];
        }
    }
    ymax = -point0[1][1];
    // Find max width of un-rotated box
    treal = tp;
    for(i = 0; i<3; ++i){
        point[i][0] = point0[i][0]*exp(flowrate*treal*dt); point[i][1] = point0[i][1]*exp(-flowrate*treal*dt);
    }
    xmax = point[2][0];
    zmax = L/2;
    // Return to undeformed box to find max height in un-rotated frame
    theta = 90.0-(45.0+THETA0);
    theta *= M_PI/180;
    point[0][0] = -sqrt(2)/2*L*cos(theta); point[0][1] = sqrt(2)/2*L*sin(theta);
    point[1][0] = -sqrt(2)/2*L*sin(theta); point[1][1] = -sqrt(2)/2*L*cos(theta);
    point[2][0] = sqrt(2)/2*L*cos(theta); point[2][1] = -sqrt(2)/2*L*sin(theta);
    for(i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
        // printf("%d %lf %lf\n",i,L1[i],L2[i]);
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
	origin = nhx*num_bins_y*num_bins_z + nhy*num_bins_z + nhz;
	printf("nx - %d ny - %d nz - %d\n",num_bins_x,num_bins_y,num_bins_z);
	printf("bx - %lf by - %lf bz - %lf\n",bin_size_x,bin_size_y,bin_size_z);
	// exit(1);
}
void allocate(){
	int i;
	rtpi = sqrt(M_PI);
	// nt_avg = 16;
	mm = malloc(sizeof(char)*100);
	// sprintf(mm, "mm/G%d_%d_%.3lf_%d_%.3lf.txt", Nb, Nc, c_norm, nt_avg,bin_size_x);
	sprintf(mm, "mm/G%d_%d_%.3lf_%.5lf_%d.txt", Nb, Nc, c_norm, flowrate,nt_avg);
	// printf("%s\n",mm);
	// mmfile = fopen(mm,"w");
	// fprintf(mmfile,"");
	// fclose(mmfile);
	D = calloc(nt_avg,sizeof(double));
	for(i=0;i<nt_avg;i++){
		D[i] = calloc(9*num_bins,sizeof(double));
	}
}
void updateLattice(unsigned long t){
	int i;
    treal = t%tp;
    for(i = 0; i<3; ++i){
        point[i][0] = point0[i][0]*exp(flowrate*treal*dt); point[i][1] = point0[i][1]*exp(-flowrate*treal*dt);
    }
    for(i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
    }
	theta = atan(L1[1]/L1[0]);
	L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
	L2xp = L2[0]*cos(theta) + L2[1]*sin(theta);
	L1yp = -L1[0]*sin(theta) + L1[1]*cos(theta);
	L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
	detL = L*(-L1yp*L2xp + L1xp*L2yp);
}
Vector3D_t invbin(int bin){
	Vector3D_t dr;
	int bin_x,bin_y,bin_z,test;
	test = bin;
	bin_x = (int)(round(bin/(num_bins_y*num_bins_z)));
	bin -= bin_x*num_bins_y*num_bins_z;
	bin_x -= nhx;
	bin_y = (int)(round(bin/(num_bins_z)));
	bin -= bin_y*num_bins_z;
	bin_y -= nhy;
	bin_z = bin;
	bin_z -= nhz;
	// printf("%d %d %d\n",bin_x+nhx,bin_y+nhyz,bin_z+nhyz);
	dr.x = bin_x*bin_size_x;
	dr.y = bin_y*bin_size_y;
	dr.z = bin_z*bin_size_z;
	// if(test==27449){
	// 	printf("bin - %d bin_x - %d bin_y - %d bin_z - %d dx - %lf dy - %lf dz - %lf\n",test,bin_x,bin_y,bin_z,dr.x,dr.y,dr.z);
	// }
	return dr;
}
