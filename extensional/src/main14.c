#include "main.h"

/////////////////////////
// This has the Ewald sum range "corrected" to be 1.0*L2yp rather than 1.5*L2yp
// Also M_HI is reduced to 3.5
////////////////////////

int main (int argc,char * argv[]){

	clock_t tic,toc;
	double start_time,time_step,ksteps;
	// int i,j,a,offset,treal;
	int i,j,startavg,startij,a,offset,offsum,offiter,treal;
	ljtimetest = 0.0;
	//printf("before init\n");
	Initialization(argc,argv); // Some initialization steps. See Initialization.c
	//printf("after init\n");
	mult_time = 0.0; bin_time = 0.0;
	start_time = omp_get_wtime();
	ksteps = omp_get_wtime();
	for(itercount=iterstart;itercount<itermax;itercount++){
		sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(mmatrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(decomp, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(beta, "decomp/B%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(ext, "prop/E%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(visc, "prop/VF%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(tim, "prop/T%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		sprintf(str2,"prop/CH%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
    sprintf(outp1,"prop/MSD%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.txt", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
		if(restart==0 || restart==2){
			xyzfile = fopen(xyz,"w");
			fprintf(xyzfile,"");
			fclose(xyzfile);
			//Charge Stuff/////
			chfile = fopen(str2,"w");
			fprintf(chfile, "");
			fclose(chfile);
			// outputfile=fopen(outp1,"w");
			// fprintf(outputfile, "");
			// fclose(outputfile);
			//////////////////
		}
		//printf("restart: %d itercount: %d\n",restart,itercount);
		initChains();
		updateLattice(tstart);
		transformCoords(); // Rotate to cartesian frame
		applyPBC();
		verletlist(); // Create the initial neighbor lists
		inverseCoords();  // Invert rotation back to the KRBC's frame
		total_time = omp_get_wtime();

		tcount=tstart%printprops;
		//printf("tstart: %d equil_start: %d\n",tstart,equil_start);

		tmax = tmax_temp;
		if(itercount==0){
            printf("No Charge Hopping in iteration 0\n");
			chargehop_type = 0;
		}
		else{
			chargehop_type = chargehop_type_temp;
		}
		

        tstart=0 ; 
		printf("Starting simulation from time step %lu -> %d\n",tstart,tmax);
	    for(t = tstart; t<(tmax+1); ++t){
				transformCoords();
				initPos();
				inverseCoords(); // Particles rotated to KRBC frame
				resetForce(); // Reset forces after each time step
	    	if(flowtype!=0){
					//if(t==equil_start){printf("Equilibrium Reached");}
	    		// printf("HI %lu\n",t);
		    	strain = (t-equil_start)*dt*flowrate;
				updateLattice(t); // Update box due to flow deformation according to KRBCs
				transformCoords(); // Particles rotated to cartesian frame
				applyPBC(); //PBC applied w/ KRBC frame
				// printprops=1;
				// if(t==0){printf("Warning - printprops = 1\n");}

				//Calculate Unwrapped Positions//
				for(i=0;i<Nc;i++){
					px[i*Nb] = 0.0; py[i*Nb] = 0.0; pz[i*Nb] = 0.0;
					for(j=0;j<Nb-1;j++){
						Vector3D_t NID = getNID(i*Nb + j,i*Nb + j + 1);
						px[i*Nb + j + 1] = px[i*Nb + j] + NID.x;
						py[i*Nb + j + 1] = py[i*Nb + j] + NID.y;
						pz[i*Nb + j + 1] = pz[i*Nb + j] + NID.z;
					}
				}
				//////////////
				if(t%printprops==0){
					calcExt(); // calculate extension
					printExt();
					//printTiming();
				}
				checkVerlet();
				bondforce(); // Calculate bonded forces from  connectivity (eg stretching, bending)
				LJforce(); // Calculate EV forces using Verlet list
				coulombForce();
				treal = t%tp;
				getNoise(); // Get new random velocities

                calcCOM(); //Calculate C.O.M of Unwrapped Chain to Apply Flow Field//

				// Update without HI if first iteration
				if(itercount==0){
					inverseCoords();  // Invert rotation back to the KRBC's frame
					calcDisFD(); //Calculate displacements
					updateFD(); //Update positions
				}
				else{
					updateHI(); // Invert rotation back to the KRBC's frame + update positions w/ mobility + store displacements
				}
				if(chargehop_type==1){
					transformCoords(); // Particles rotated to cartesian frame
					chargeHop();
					inverseCoords();  // Invert rotation back to the KRBC's frame
				}
				else if(chargehop_type==2){
					printf("Warning - Displacements not yet configured for this case, exiting\n");
                    exit(1); 
					transformCoords(); // Particles rotated to cartesian frame
					chargeHop_nonint();
					inverseCoords();  // Invert rotation back to the KRBC's frame
				}

				if(t%printprops==0){
					calcCH();
					for(i=0;i<Ncharges*printprops;i++){
						dxt[i] = 0;
						dyt[i] = 0;
						dzt[i] = 0;
					}
					tcount=0;
				}
				// Sample CA TEA parameters
				// Note CA of TEA parameters only performed for iter = 0
				if((t-equil_start)%strain_sample==0 && itercount==0 && strain<=strain_ss){
					ewaldBin();
					CATEA();
					printTEA();
				}
				else if(itercount==0 && strain>strain_ss && t%10000==0){
					ewaldBin();
					CATEA();
					printTEA();
				}
				//printperiod=1 ;
                //if(itercount==1){printperiod=1 ;}
				//if(t==0){printf("Warning - printperiod = 1\n");}
				if(t%printperiod==0){
					if(restart==0 || restart==2 || (restart==1 && t>tstart)){
						printTrajectory(); // Save trajectory to xyzfile every tau = 1/dt
					}
				}

				if(t==equil_start){
					xyzstart = 0; // start code
					printRestart();
				}
				// if(t%1000==0 || t==tstart){
				// 	// printf("ts - %lu strain - %lf 1k ts in %lf seconds, ext_avg - %lf, rx[100] = %.12lf\n",t,strain,omp_get_wtime()-ksteps,ext_avg,rx[100]);
				// 	printf("ts - %lu strain - %lf 1k ts in %lf s, ext_avg - %lf, bin_timepts - %lf s, mult_timepts - %lf s\n",t,strain,omp_get_wtime()-ksteps,ext_avg,bin_time/((double)t-(double)tstart),mult_time/((double)t-(double)tstart));
				// 	ksteps = omp_get_wtime();
				// 	// mult_time = 0.0;
				// 	// bin_time = 0.0;
				// }
	    	}
	    	else{
	    		// printf("FD %lu\n",t);
				if(t==0){printf("Equilibrium Simulation\n");}
	    		strain = 0.0;
                resetForce(); // Reset forces after each time step
	    		updateLattice(0);
	    		transformCoords();
	    		applyPBC();
	    		checkVerlet();
	    		bondforce();
	    		LJforce();
				coulombForce();
	    		getNoise();
                calcCOM(); //Calculate C.O.M of Unwrapped Chain to Apply Flow Field//
                // Update without HI if first iteration
				if(itercount==0){
					inverseCoords();  // Invert rotation back to the KRBC's frame
					calcDisFD(); //Calculate displacements
					updateFD(); //Update positions
				}
				else{
					updateHI(); // Invert rotation back to the KRBC's frame + update positions w/ mobility + store displacements
				}
				if(chargehop_type==1){
					transformCoords(); // Particles rotated to cartesian frame
					chargeHop();
					inverseCoords();  // Invert rotation back to the KRBC's frame
				}
				else if(chargehop_type==2){
					printf("Warning - Displacements not yet configured for this case, exiting\n");
                    exit(1); 
					transformCoords(); // Particles rotated to cartesian frame
					chargeHop_nonint();
					inverseCoords();  // Invert rotation back to the KRBC's frame
				}

				if(t%printprops==0){
					calcCH();
					for(i=0;i<Ncharges*printprops;i++){
						dxt[i] = 0;
						dyt[i] = 0;
						dzt[i] = 0;
					}
					tcount=0;
				}
                // Sample CA TEA parameters
				// Note CA of TEA parameters only performed for iter = 0
				if((t-equil_start)%strain_sample==0 && itercount==0 && strain<=strain_ss){
					ewaldBin();
					CATEA();
					printTEA();
				}
				else if(itercount==0 && strain>strain_ss && t%10000==0){
					ewaldBin();
					CATEA();
					printTEA();
				}

	    		// if(t%1000==0){
	    		// 	printf("ts - %lu strain - %lf 1k ts in %lf seconds, fx[0] - %lf, rx[100] = %.12lf\n",t,strain,omp_get_wtime()-ksteps,fx[0],rx[100]);
	    		// 	ksteps = omp_get_wtime();
	    		// 	mult_time = 0.0;
	    		// }
	    		if(t==equil_start){
	    			xyzstart = 0; // start code
	    			printRestart();
	    		}
				if(t%printperiod==0){
					if(restart==0 || restart==2 || (restart==1 && t>tstart)){
						printTrajectory(); // Save trajectory to xyzfile every tau = 1/dt
					}
				}
	    	}
			printVisc();
			// printTiming();
	    }
	    xyzstart = 1; // end code
	    printRestart(); // print the final frame for restarting simulations eg from steady to relax
	    tstart = 0;
			//printf("Warning - calcMSD() turned off\n");
			if(itercount>0) {
				calcMSD();
			}
			//calcMSDcom(); //Calculate MSD of C.O.M
	    if(itercount==0){
            printTEA();
            resetAverage();
		}
			//restart = 0;
			printf("%lu %.3e seconds \n",t,omp_get_wtime()-start_time);
		}
	    return 0;
}

void initChains(){
    int i,j,k,index,monindex,test,tempN,initcount,inittest;
    double phi,theta,dx,dy,dz;
    char tempname[1024];
    // If starting a new run, clear the xyz file and randomly distribute non-overlapping chains in the box
    // if(itercount==0 && (restart==1 || restart==2)){
    if(itercount==0 && restart==1){
    	xyzfile = fopen(xyz,"r");
    	if(!xyzfile){
    	    printf("Error: couldn't find xyz file %s to initialize from for code restart = %d or itercount = %d\n",xyz,restart,itercount);
    	    printf("Try checking input or change to restart = 0\n");
    	    exit(1);
    	}
    	while(!feof(xyzfile)){
    	    fscanf(xyzfile, "%d\n", &tempN);
    	    fscanf(xyzfile, "%lu\n", &tstart);
    	    for(j=0;j<tempN;j++){
    	        fscanf(xyzfile, "A %lf %lf %lf\n",&rx[j],&ry[j],&rz[j]);
    	    }
    	    // if(ttemp==tstart){
    	    //     break;
    	    // }
    	}
    	tstart++;
    }
    else if(itercount>0 && (restart==0 || restart==2)){
			if(restart==0){
				sprintf(xyz, "res/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg, trace);
			}
			else if (restart==2){
				sprintf(xyz, "res/SS%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg,trace);
			}
			//sprintf(xyz, "res/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg, trace);


	    xyzfile = fopen(xyz,"r");
	    if(!xyzfile){
	        printf("Error: couldn't find xyz file %s to initialize from for code restart = %d or itercount = %d\n",xyz,restart,itercount);
	        // printf("Try checking input or change to restart = 0\n");
	        exit(1);
	    }
	    while(!feof(xyzfile)){
	        fscanf(xyzfile, "%d\n", &tempN);
	        fscanf(xyzfile, "%lu\n", &tstart);
	        for(j=0;j<tempN;j++){
	            fscanf(xyzfile, "A %lf %lf %lf\n",&rx[j],&ry[j],&rz[j]);
	        }
	        // if(ttemp==tstart){
	        //     break;
	        // }
	    }
	    //tstart--;
		sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);


    	xyzfile = fopen(xyz,"w");
    	fprintf(xyzfile,"");
    	fclose(xyzfile);
	    extfile = fopen(ext,"w");
	    fprintf(extfile,"");
	    fclose(extfile);
	    viscfile = fopen(visc,"w");
	    fprintf(viscfile,"");
	    fclose(viscfile);
	    betafile = fopen(beta,"w");
	    fprintf(betafile,"");
	    fclose(betafile);
    }
    else if(itercount>0 && restart==1){
			sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
			xyzfile = fopen(xyz,"r");
        if(!xyzfile){
        	printf("Error: couldn't find xyz file %s to initialize from for code restart = %d or itercount = %d\n",xyz,restart,itercount);
            printf("Try checking input or change to restart = 0\n");
            exit(1);
        }
        while(!feof(xyzfile)){
            fscanf(xyzfile, "%d\n", &tempN);
            fscanf(xyzfile, "%lu\n", &tstart);
            for(j=0;j<tempN;j++){
                fscanf(xyzfile, "A %lf %lf %lf\n",&rx[j],&ry[j],&rz[j]);
            }
            // if(ttemp==tstart){
            //     break;
            // }
        }
        fclose(xyzfile);
        // sprintf(xyz, "xyz/R%d_%d_%.3lf_%.5lf_%d_%d_%d.xyz", Nb, Nc, c_norm, flowrate, nt_avg, trace, itercount);
        tstart++;
    }
    // else if(restart==0 || restart==2 || (restart==1 && itercount>iterstart)){
    else if(itercount==0 && (restart==0 || restart==2)){
    // else if(itercount==0 && restart==0){
		//sprintf(xyz, "res/R%d_%d_%.3lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, flowrate, nt_avg, trace);
	    //xyzfile = fopen(xyz,"r");
	    // if(xyzfile){
	    //     // printf("Error: couldn't find xyz file %s to initialize from for code restart = %d or itercount = %d\n",xyz,restart,itercount);
	    //     // printf("Try checking input or change to restart = 0\n");
	    //     // exit(1);
	    //     // while(!feof(xyzfile)){
	    //     //     fscanf(xyzfile, "%d\n", &tempN);
	    //     //     fscanf(xyzfile, "%lu\n", &tstart);
	    //     //     for(j=0;j<tempN;j++){
	    //     //         fscanf(xyzfile, "A %lf %lf %lf\n",&rx[j],&ry[j],&rz[j]);
	    //     //     }
	    //     //     // if(ttemp==tstart){
	    //     //     //     break;
	    //     //     // }
	    //     // }
	    //     //tstart++;
	    //     sprintf(xyz, "xyz/R%d_%d_%.3lf_%.5lf_%d_%d_%d.xyz", Nb, Nc, c_norm, flowrate, nt_avg, trace, itercount);
	    //     //if(restart==2){
	    //   	xyzfile = fopen(xyz,"w");
	    //   	fprintf(xyzfile,"");
	    //   	fclose(xyzfile);
	    //   	extfile = fopen(ext,"w");
	    //   	fprintf(extfile,"");
	    //   	fclose(extfile);
	    //   	viscfile = fopen(visc,"w");
	    //   	fprintf(viscfile,"");
	    //   	fclose(viscfile);
	    //   	betafile = fopen(beta,"w");
	    //   	fprintf(betafile,"");
	    //   	fclose(betafile);
	    //     //}
	    // }
	   // else{
		 sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
 		 xyzfile = fopen(xyz,"w");
	        fprintf(xyzfile,"");
	        fclose(xyzfile);
	        extfile = fopen(ext,"w");
	        fprintf(extfile,"");
	        fclose(extfile);
	        viscfile = fopen(visc,"w");
	        fprintf(viscfile,"");
	        fclose(viscfile);
	        betafile = fopen(beta,"w");
	        fprintf(betafile,"");
	        fclose(betafile);
	        tstart = 0;
	        inittest = 0;
	        while(inittest==0){
	            inittest = 1;
	            for(i = 0; i<Nc; ++i){
	                initcount = 0;
	                test = 0;
	                index = Nb*i; // Index of first monomer in  i
	                while(test==0){
	                    test = 1;
											if(Nc==1){
		                    rx[index] = 0;//ran1(idum)*box_length - box_side;
		                    ry[index] = 0;//ran1(idum)*box_length - box_side;
		                    rz[index] = 0;//ran1(idum)*box_length - box_side;
		                  }
		                  else{
		                    rx[index] = ran1(idum)*L - box_side;
		                    ry[index] = ran1(idum)*L - box_side;
		                    rz[index] = ran1(idum)*L - box_side;
		                  }
	                    for(k = 0; k<index; k++){
	                        dx = rx[index]-rx[k]; dy = ry[index]-ry[k]; dz = rz[index]-rz[k];
	                        dx -= L*round(dx/L); dy -= L*round(dy/L); dz -= L*round(dz/L);
	                        if(dx*dx+dy*dy+dz*dz<10.0)
	                        {
	                            test = 0;
	                        }
	                    }
	                    initcount++;
	                    if(initcount>1e3){
	                        // printf("max attempts reached\n");
	                        // exit(1);
	                        inittest = 0;
	                        break;
	                    }
	                }
	                if(inittest==0){
	                    break;
	                }
	                for(j=1;j<Nb;j++){
	                    monindex = index+j; // Index of monomer j in  i
	                    test = 0;
	                    while(test==0){
	                        test = 1;
	                        theta = ran1(idum)*2.0*M_PI;
	                        phi = acos(2.0*ran1(idum)-1.0);
	                        rx[monindex] = rx[monindex-1] + 2.05*sin(phi)*cos(theta);
	                        ry[monindex] = ry[monindex-1] + 2.05*sin(phi)*sin(theta);
	                        rz[monindex] = rz[monindex-1] + 2.05*cos(phi);
	                        rx[monindex] -= L*round(rx[monindex]/L);
	                        ry[monindex] -= L*round(ry[monindex]/L);
	                        rz[monindex] -= L*round(rz[monindex]/L);
	                        for(k = 0; k<monindex; k++){
	                            dx = rx[monindex]-rx[k]; dy = ry[monindex]-ry[k]; dz = rz[monindex]-rz[k];
	                            dx -= L*round(dx/L); dy -= L*round(dy/L); dz -= L*round(dz/L);
	                            if(dx*dx+dy*dy+dz*dz<4.0)
	                            {
	                                test = 0;
	                            }
	                        }
	                        initcount++;
	                        if(initcount>1e3){
	                            // printf("max attempts reached\n");
	                            // exit(1);
	                            inittest = 0;
	                            break;
	                        }
	                    }
	                    if(inittest==0){
	                        break;
	                    }
	                }
	                if(inittest==0){
	                    break;
	                }
	            }
	        }
	        //Inverse Transform
	        double rxi,ryi;
	        double theta = atan(L1[1]/L1[0]);
	        for(i=0; i<N; ++i){
	            rxi = rx[i]; ryi = ry[i];
	            rx[i] = rxi*cos(theta) - ryi*sin(theta);
	            ry[i] = rxi*sin(theta) + ryi*cos(theta);
	        }
	    //}
    }
    else{
        printf("Error - Initialize chains restart code %d %d not recognized\n",restart,itercount);
        exit(1);
    }

		//Charge Stuff//
		for(i=0;i<N;i++){
			Charge[i] = 0 ;
		}
		for(i=0;i<Ncharges;i++){
			cnumber=ran1(idum)*N;
			if(Charge[cnumber]==0){
				Charge[cnumber]=1;
				Charge_indices[i]=cnumber;
				Charge_track[i]=cnumber;
				Charge_old[i]=cnumber;
			}
			else
			{
				i=i-1;
			}
		}
		////
}
void resetForce(){
	int i;
	for(i = 0; i<N; ++i){
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
	}
}
float gasdev(long *idum){
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset = 0;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0*ran1(idum)-1.0;
			v2 = 2.0*ran1(idum)-1.0;
			rsq = v1*v1+v2*v2;
		}
		while(rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}
void bondforce(){
	double dx1, dy1, dz1, rr, r, Fs, dx2, dy2, dz2, amp1, amp2, var, thetab, coeff,dx,dy,dz;
	int i, j;
	// unsigned long tmod = t%printperiod;
	for(i = 0; i<Nc; ++i){
		for(j = 1; j < Nb; j++){
			Vector3D_t NID = getNID(Nb*i+j,Nb*i+j-1);
			// r = sqrt(NID.x*NID.x + NID.y*NID.y + NID.z*NID.z);
			rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
			r = sqrt(rr);
			if(strcmp(spring,"FENE")==0){
				if(rr>qmax2){
					// inverseCoords();
					// printTrajectory();
					printf("FENE overstretch %lu %d %d %lf %lf %lf %lf %lf %ld\n",t,i,j,sqrt(rr),qmax2,NID.x,NID.y,NID.z,*idum);
					exit(1);
				}
				Fs = -kappas/(1.0-rr/qmax2);
			}
			else if(strcmp(spring,"Hookean")==0){
				Fs = -kappas*(r-qmax)/r;
			}
			else{
				printf("Error - spring type %s not recognized, exiting\n",spring);
				exit(1);
			}
			// Fs = -kappas*(r-2.0)/r;

			fx[Nb*i+j] += Fs*NID.x;
			fy[Nb*i+j] += Fs*NID.y;
			fz[Nb*i+j] += Fs*NID.z;
			fx[Nb*i+j-1] += -Fs*NID.x;
			fy[Nb*i+j-1] += -Fs*NID.y;
			fz[Nb*i+j-1] += -Fs*NID.z;
			tauxx -= 2.0*NID.x*Fs*NID.x;
			tauyy -= 2.0*NID.y*Fs*NID.y;
		}
	}
	if(kappab>0.0){
		// #pragma omp parallel
		{
			// int i,j;
			// double dx1, dy1, dz1, rr, r, Fs, dx2, dy2, dz2, amp1, amp2, var, thetab, coeff;
			// #pragma omp for schedule(static)
				for(i = 0; i<Nc; ++i){
					for(j = 1; j < Nb-1; j++){
						Vector3D_t NID1 = getNID(Nb*i+j,Nb*i+j-1);
						Vector3D_t NID2 = getNID(Nb*i+j,Nb*i+j+1);
						amp1 = sqrt(NID1.x*NID1.x + NID1.y*NID1.y + NID1.z*NID1.z);
						amp2 = sqrt(NID2.x*NID2.x + NID2.y*NID2.y + NID2.z*NID2.z);
						var = (NID1.x*NID2.x + NID1.y*NID2.y + NID1.z*NID2.z)/(amp1*amp2);
						thetab = M_PI - acos(var); // Bond angle, Flory coordNb*i+jnates conventNb*i+jon
						coeff = -kappab*thetab/sqrt(1-var*var); // theta_0 = 0 Nb*i+js the relaxed state
						// #pragma omp atomic
							fx[Nb*i+j] += coeff*((NID1.x+NID2.x)/(amp1*amp2)-var*(NID1.x/(amp1*amp1)+NID2.x/(amp2*amp2)));
							fx[Nb*i+j-1] += coeff*(var*NID1.x/(amp1*amp1)-NID2.x/(amp1*amp2));
							fx[Nb*i+j+1] += coeff*(var*NID2.x/(amp2*amp2)-NID1.x/(amp1*amp2));
							fy[Nb*i+j] += coeff*((NID1.y+NID2.y)/(amp1*amp2)-var*(NID1.y/(amp1*amp1)+NID2.y/(amp2*amp2)));
							fy[Nb*i+j-1] += coeff*(var*NID1.y/(amp1*amp1)-NID2.y/(amp1*amp2));
							fy[Nb*i+j+1] += coeff*(var*NID2.y/(amp2*amp2)-NID1.y/(amp1*amp2));
							fz[Nb*i+j] += coeff*((NID1.z+NID2.z)/(amp1*amp2)-var*(NID1.z/(amp1*amp1)+NID2.z/(amp2*amp2)));
							fz[Nb*i+j-1] += coeff*(var*NID1.z/(amp1*amp1)-NID2.z/(amp1*amp2));
							fz[Nb*i+j+1] += coeff*(var*NID2.z/(amp2*amp2)-NID1.z/(amp1*amp2));
							tauxx += NID1.x*coeff*(var*NID1.x/(amp1*amp1)-NID2.x/(amp1*amp2));
							tauxx -= NID1.x*coeff*((NID1.x)/(amp1*amp2)-var*(NID1.x/(amp1*amp1)));
							tauxx += NID2.x*coeff*((NID2.x)/(amp1*amp2)-var*(NID2.x/(amp2*amp2)));
							tauxx -= NID2.x*coeff*(var*NID2.x/(amp2*amp2)-NID1.x/(amp1*amp2));
							tauyy += NID1.y*coeff*(var*NID1.y/(amp1*amp1)-NID2.y/(amp1*amp2));
							tauyy -= NID1.y*coeff*((NID1.z)/(amp1*amp2)-var*(NID1.z/(amp1*amp1)));
							tauyy += NID2.y*coeff*((NID2.z)/(amp1*amp2)-var*(NID2.z/(amp2*amp2)));
							tauyy -= NID2.y*coeff*(var*NID2.y/(amp2*amp2)-NID1.y/(amp1*amp2));
					}
				}
		}
	}
}
void LJforce(){
	double dx, dy, dz, rr, coeff, r6, ratio;
	int i, j, k;
	// #pragma omp parallel for schedule(static,1) private(j,k,dx,dy,dz,rr,coeff,r6,ratio) // not thread safe, need critical section
	unsigned long tmod = t%printperiod;
	#pragma omp parallel
	{
		int j,k;
		double dx,dy,dz,rr,coeff,r6,ratio,tauxxt,tauyyt;
		double *fxt = calloc(N,sizeof(double));
		double *fyt = calloc(N,sizeof(double));
		double *fzt = calloc(N,sizeof(double));
		tauxxt = 0.0; tauyyt = 0.0;
		#pragma omp for schedule(static,1)
			for(i = 0; i<N; ++i){
				for(j = 0; j<nlist[i]; ++j){
					k = list[i*N+j];
					Vector3D_t NID = getNID(i,k);
					rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
					if(rr<r2cut){
						// if(sqrt(rr)<1.5){
						// 	printf("%lu %d %d %lf\n",t,i,k,sqrt(rr));
						// 	exit(1);
						// }
						ratio = 4.00/rr;
						r6 = ratio*ratio*ratio;
						// if(r6>3){
						// 	if(t<1e6){
						// 		r6 = 3;
						// 	}
						// 	else if(t<1e7){
						// 		if(r6>5) r6=5;
						// 	}
						// }
						coeff = (48.0*epsilon/rr)*(r6*r6-0.5*r6); //Lennard Jones
						 //coeff = (12*epsilon/rr)*(r6*r6-r6);
						// #pragma omp atomic
							fxt[i] += coeff*NID.x;
							fyt[i] += coeff*NID.y;
							fzt[i] += coeff*NID.z;
							fxt[k] -= coeff*NID.x;
							fyt[k] -= coeff*NID.y;
							fzt[k] -= coeff*NID.z;
							tauxxt -= 2.0*NID.x*coeff*NID.x;
							tauyyt -= 2.0*NID.y*coeff*NID.y;
					}
				}
			}
		#pragma omp critical
		{
			for(i=0;i<N;i++){
				fx[i] += fxt[i];
				fy[i] += fyt[i];
				fz[i] += fzt[i];
			}
			tauxx += tauxxt;
			tauyy += tauyyt;
			free(fxt); free(fyt); free(fzt);
		}
	}
}
void verletlist(){
	int i, j;
	double dx, dy, dz, r,start_time;
	for(i=0;i<N;i++){
		nlist[i] = 0; // Reset lists on update
		dxv[i] = 0.0; dyv[i] = 0.0; dzv[i] = 0.0;
	}
	#pragma omp parallel
	{
		// int id = omp_get_thread_num();
		int i,j;
		double dx,dy,dz,r;
		#pragma omp for schedule(static,1)
			for(i=0;i<N-1;i++){
				for(j=i+1;j<N;j++){
					Vector3D_t NID = getNID(i,j);
					r = sqrt(NID.x*NID.x + NID.y*NID.y + NID.z*NID.z);
					if(r < rv){
						list[i*N+nlist[i]] = j;
						// #pragma omp atomic
							nlist[i]++;
					}
				}
			}
	}
}
void getNoise(){
	int i;
	for(i = 0; i<N*3; ++i){
		R[i] = gasdev(idum);
	}
}
void updateEquil(){
	int i;
	double dx,dy,dz;
	// printf("%lf\n",rx[100]);
	for(i=0;i<N;i++){
		dx = 0.0; dy = 0.0; dz = 0.0;
		dx = dt_eq*fx[i] + pc_eq*R[3*i];
		dy = dt_eq*fy[i] + pc_eq*R[3*i+1];
		dz = dt_eq*fz[i] + pc_eq*R[3*i+2];
		rx[i] += dx;
		ry[i] += dy;
		rz[i] += dz;
		dxv[i] += dx; dyv[i] += dy; dzv[i] += dz;
	}
	// printf("%lf\n",rx[100]);
	// exit(1);
}
void updateFD(){
	int i,j;
	double dx,dy,dz,COMx,COMy,COMz;
    double drelxi,drelyi,drelx,drely;
    double theta = atan(L1[1]/L1[0]);

	//tcount++ ;
	//if(t==0){printf("Warning - tcount++ is in UpdateFD\n");}
	for(i=0;i<N;i++){
		dx = 0.0; dy = 0.0; dz = 0.0;
        drelx = 0; drely = 0; drelxi = 0; drelyi = 0;
		
        dx = dt*(fx[i] + flowrate*rx[i]) + pc*R[3*i];
		dy = dt*(fy[i] - flowrate*ry[i]) + pc*R[3*i+1];
		dz = dt*fz[i] + pc*R[3*i+2];

		// j = Charge_track[0];
		// if(i==j){
		// 	//printf("force_u: %.9lf %.9lf %.9lf\n",dt*fx[i],dt*fy[i],dt*fz[i]);
		// 	//printf("flow_u: %.9lf %.9lf %.9lf\n",dt*flowrate*(rx[i]-comx[0]),dt*flowrate*(ry[i]-comy[0]),0.0);

		// 	//printf("total_u: %.9lf %.9lf %.9lf\n",dx,dy,dz);
		// 	//printf("force_u: %.9lf %.9lf %.9lf\n",dt*fx[i],dt*fy[i],dt*fz[i]);
		// 	//printf("flow_u: %.9lf %.9lf %.9lf\n",dt*flowrate*(rx[i]-comx[0]),dt*flowrate*(ry[i]-comy[0]),0.0);

		// 	// printf("pos_u: %lu %.9lf %.9lf %.9lf\n",t,rx[j],ry[j],rz[j]);
		// 	// printf("map_u: %lu %.9lf %.9lf %.9lf\n",t,sx[j],sy[j],sz[j]);
		// }

		rx[i] += dx;
		ry[i] += dy;
		rz[i] += dz;
		dxv[i] += dx; dyv[i] += dy; dzv[i] += dz;
	}

	//calcCOM(); //Calculates center of mass of each chain


	// Subtract off the center off mass//
	// for(i=0;i<Nb;++i){
	// 	rx[i] -= comx[0];
	// 	ry[i] -= comy[0];
	// 	rz[i] -= comz[0];
	// }

	//if(t==0){printf("Warning - C.O.M subtracted from chain\n");}

}
void updateHI(){
	int i,j,k,l,m,n,a,startavg,startij,startji,bin,ii,jj,ext_bin,strain_bin;
	double dx,dy,dz,start_time,betaii;
	int offset,offsum,offiter;
	// ------------------- 2) Construct square MM then multiply ---------------------
	double *dxu = calloc(N, sizeof(double));
	double *dyu = calloc(N, sizeof(double));
	double *dzu = calloc(N, sizeof(double));
	if(t%srpy==0 || t==tstart){
		start_time = omp_get_wtime();
		calcHI();
		bin_time += omp_get_wtime() - start_time;
		strain_bin = (int)floor(strain/strain_interval);
		if(strain_bin>sbins){
			strain_bin = sbins; // steady state bin
		}
		bt = bp[strain_bin];
		for(i=0;i<3*N;++i){
			Ct[i] = C[strain_bin][i%(3*Nb)];
		}
		// printf("test binning\n");
		// if(t%1000==0){
			// printf("bin - %d ts - %lu bt - %lf\n",strain_bin,t,bt);
		// 	for(i=0;i<N;++i){
		// 		printf("%d %lf %lf %lf\n",i,Ct[3*i],Ct[3*i+1],Ct[3*i+2]);
		// 	}
		// 	for(i=0;i<3*N;++i){
		// 		printf("%d %lf\n",i,M[i]);
		// 	}
		// 	exit(1);
		// }
		// printf("%lu bin time %lf\n",t,omp_get_wtime()-start_time);
	}
	inverseCoords();  // Invert rotation back to the KRBC's frame
	start_time = omp_get_wtime();
	#pragma omp parallel
	{
		double *dxn = calloc(N, sizeof(double));
		double *dyn = calloc(N, sizeof(double));
		double *dzn = calloc(N, sizeof(double));
		int m,n,j,k,startij;
		double betaii;
		#pragma omp for schedule(dynamic)
			for(i=0;i<N;i++){
				m = i*3;
				for(j=i+1;j<N;j++){
					n = j*3;
					startij = 9*N*i+9*j;
					dxn[i] += dt*(M[startij]*(fx[j]+p*bt*Ct[m]*R[n])+M[startij+1]*(fy[j]+p*bt*Ct[m]*R[n+1])+M[startij+2]*(fz[j]+p*bt*Ct[m]*R[n+2]));
					dyn[i] += dt*(M[startij+3]*(fx[j]+p*bt*Ct[m+1]*R[n])+M[startij+4]*(fy[j]+p*bt*Ct[m+1]*R[n+1])+M[startij+5]*(fz[j]+p*bt*Ct[m+1]*R[n+2]));
					dzn[i] += dt*(M[startij+6]*(fx[j]+p*bt*Ct[m+2]*R[n])+M[startij+7]*(fy[j]+p*bt*Ct[m+2]*R[n+1])+M[startij+8]*(fz[j]+p*bt*Ct[m+2]*R[n+2]));
					dxn[j] += dt*(M[startij]*(fx[i]+p*bt*Ct[n]*R[m])+M[startij+1]*(fy[i]+p*bt*Ct[n]*R[m+1])+M[startij+2]*(fz[i]+p*bt*Ct[n]*R[m+2]));
					dyn[j] += dt*(M[startij+3]*(fx[i]+p*bt*Ct[n+1]*R[m])+M[startij+4]*(fy[i]+p*bt*Ct[n+1]*R[m+1])+M[startij+5]*(fz[i]+p*bt*Ct[n+1]*R[m+2]));
					dzn[j] += dt*(M[startij+6]*(fx[i]+p*bt*Ct[n+2]*R[m])+M[startij+7]*(fy[i]+p*bt*Ct[n+2]*R[m+1])+M[startij+8]*(fz[i]+p*bt*Ct[n+2]*R[m+2]));
				}
				startij = 9*N*i+9*i;
				dxn[i] += dt*(M[startij]*(fx[i]+p*Ct[m]*R[m])+M[startij+1]*(fy[i]+p*Ct[m]*bt*R[m+1])+M[startij+2]*(fz[i]+p*Ct[m]*bt*R[m+2]));
				dyn[i] += dt*(M[startij+3]*(fx[i]+p*Ct[m+1]*bt*R[m])+M[startij+4]*(fy[i]+p*Ct[m+1]*R[m+1])+M[startij+5]*(fz[i]+p*Ct[m+1]*bt*R[m+2]));
				dzn[i] += dt*(M[startij+6]*(fx[i]+p*Ct[m+2]*bt*R[m])+M[startij+7]*(fy[i]+p*Ct[m+2]*bt*R[m+1])+M[startij+8]*(fz[i]+p*Ct[m+2]*R[m+2]));
			}
		#pragma omp critical
		{
			for(k=0;k<N;k++){
				dxu[k] += dxn[k];
				dyu[k] += dyn[k];
				dzu[k] += dzn[k];
			}
			free(dxn); free(dyn); free(dzn);
		}
	}
    //////Store Displacement from Forces for Charges////
    double drelxi=0;
    double drelyi=0;
    int chain;

    tcount++;
    for(i=0;i<Ncharges;i++){
        j = Charge_track[i];
		dx = 0.0; dy = 0.0; dz = 0.0;

        //Calculate relative displacement between unwrapped chain + C.O.M
		chain = (int) floor(j/Nb) ; 
		//printf("j: %d chain: %d\n",j,chain);
        drelxi = comx[chain]-px[j] ; 
        drelyi = comy[chain]-py[j] ; 

        dx = dxu[j] + dt*flowrate*drelxi ; 
        dy = dyu[j] - dt*flowrate*drelyi ; 
        dz = dzu[j] ; 

        dxt[i*printprops+tcount] = dx;
		dyt[i*printprops+tcount] = dy;
		dzt[i*printprops+tcount] = dz;
    
    }
    ////////////////////////////////////////////////
	mult_time += omp_get_wtime() - start_time;
	// for(i=0;i<N;++i){
	// 	printf("%d %lf %lf %lf\n",i,dxu[i],dyu[i],dzu[i]);
	// }
	// exit(1);
	for(i=0;i<N;i++){
		dxu[i] += dt*flowrate*rx[i];
		dyu[i] -= dt*flowrate*ry[i];
		dxv[i] += dxu[i];
		rx[i] += dxu[i];
		dyv[i] += dyu[i];
		ry[i] += dyu[i];
		dzv[i] += dzu[i];
		rz[i] += dzu[i];
	}
	free(dxu); free(dyu); free(dzu);
}
void checkVerlet(){
	// Update when dr > 1/2(rv-rc)
	int i,j;
	double drcheck,dx,dy,dz;
	for(i=0;i<N;i++){
		drcheck = sqrt(dxv[i]*dxv[i] + dyv[i]*dyv[i] + dzv[i]*dzv[i]);
		if(drcheck > rnew){ // Condition met, need to update lists
			verletlist();
			break;
		}
	}
}
void printTrajectory(){
    //if(t==0){printf("Warning - xyz set up for unwrapped positions, px,py,pz\n");}
	int i,j,ind;
	xyzfile = fopen(xyz, "a");
	fprintf(xyzfile, "%d\n", N);
	fprintf(xyzfile,"Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\" Origin=\"%lf %lf %lf\" Properties=\"\" Time=%lu\n",L1[0],L1[1],0.0,L2[0],L2[1],0.0,0.0,0.0,L,-(L1[0]+L2[0])/2,-(L1[1]+L2[1])/2,-L/2,t);
	for(i = 0; i<Nc; ++i){
		for(j=0;j<Nb;j++){
				ind = Nb*i +j;
				if(Charge[ind] == 0){
					fprintf(xyzfile, "A %lf %lf %lf\n", rx[Nb*i+j], ry[Nb*i+j], rz[Nb*i+j]);
				}
				else if(Charge[ind] == 1){
					fprintf(xyzfile, "B %lf %lf %lf\n", rx[Nb*i+j], ry[Nb*i+j], rz[Nb*i+j]);
				}
		}
	}
    // for(i = 0; i<Nc; ++i){
	// 	for(j=0;j<Nb;j++){
	// 			ind = Nb*i +j;
	// 			if(Charge[ind] == 0){
	// 				fprintf(xyzfile, "A %lf %lf %lf\n", px[Nb*i+j], py[Nb*i+j], pz[Nb*i+j]);
	// 			}
	// 			else if(Charge[ind] == 1){
	// 				fprintf(xyzfile, "B %lf %lf %lf\n", px[Nb*i+j], py[Nb*i+j], pz[Nb*i+j]);
	// 			}
	// 	}
	// }
	// fprintf(xyzfile, "M %lf %lf %lf\n", comx[0], comy[0], comz[0]) ; //Unwrapped Center of Mass
	//fprintf(xyzfile, "M %lf %lf %lf\n", comx_r[0], comy_r[0], comz_r[0]) ; //Real Center of Mass
	fclose(xyzfile);
}
void printRestart(){
	int i,j;
	if(xyzstart==0){
		sprintf(xyz, "res/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg, trace);

	}
	else if(xyzstart==1){
		sprintf(xyz, "res/SS%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg, trace);

	}
	else{
		printf("Warning, unrecognized xyzstart value %d\n",xyzstart);
		sprintf(xyz, "res/TMP%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, nt_avg, trace);
		printf("Assigning temporary filename %s - you might want to fix this\n",xyz);
	}
	xyzfile = fopen(xyz,"w");
	fprintf(xyzfile, "%d\n%lu\n", N, t);
	for(i = 0; i<Nc; ++i){
		for(j=0;j<Nb;j++){
				fprintf(xyzfile, "A %lf %lf %lf\n", rx[Nb*i+j], ry[Nb*i+j], rz[Nb*i+j]);
		}
	}
	fclose(xyzfile);
	sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz", Nb, Nc, c_norm, L,Ncharges,barrier,barrier2,lambda_d,flowrate, trace,itercount);
}
void calcExt(){
	int i,j;
	double maxx,minx,rxi,ryi;
    // double theta = atan(L1[1]/L1[0]);
	// for(i=0;i<Nc;i++){
	// 	px[i*Nb] = 0.0; py[i*Nb] = 0.0; pz[i*Nb] = 0.0;
	// 	for(j=0;j<Nb-1;j++){
	// 		Vector3D_t NID = getNID(i*Nb + j,i*Nb + j + 1);
	// 		px[i*Nb + j + 1] = px[i*Nb + j] + NID.x;
	// 		py[i*Nb + j + 1] = py[i*Nb + j] + NID.y;
	// 		pz[i*Nb + j + 1] = pz[i*Nb + j] + NID.z;
	// 	}
	// }
    // for(i = 0; i<N; ++i){
    //     rxi = px[i]; ryi = py[i];
    //     px[i] = rxi*cos(theta) + ryi*sin(theta);
    //     py[i] = -rxi*sin(theta) + ryi*cos(theta);
    // }
	for(i=0;i<Nc;i++){
		maxx = px[i*Nb]; minx = px[i*Nb];
		for(j=1;j<Nb;j++){
			maxx = max(maxx,px[i*Nb + j]);
			minx = min(minx,px[i*Nb + j]);
		}
		rext[i] = maxx - minx;  //Extension in x-direction only
	}
	ext_avg = 0.0;
	for(i=0;i<Nc;++i){
		ext_avg += rext[i];
	}
	ext_avg /= ((Nb-1)*qmax*Nc);
}
void printExt(){
	int i;
	double extt;
	extfile = fopen(ext,"a");
	//fprintf(extfile,"%lu %lf \n",t,strain);
	for(i=0;i<Nc;i++){
		extt = rext[i]/(Nb*qmax) ;
		fprintf(extfile,"%lu %lf %lf %lf\n",t,strain,extt,pow(extt,2)); //time, strain, x/L, (x/L)^2
	}
	//fprintf(extfile,"\n");
	fclose(extfile);
}
void ewaldBin(){
	int a,i,j,ii,kx,ky,kz,offset,start,kxmax,kymax,kzmax,nmax,nx,ny,nz,nn,kcount,ylower,zlower,yupper,zupper,treal,startii,startij,startji,bij,bji;
	double alpha,rkx,rky,rkz,rkk,m2,kcoeff,rkxx,rkxy,rtpi,coskr,selfconst,alpha2,alpha3,alpha4,alpha5,alpha7,rce,rkcut,rkkcut,M_HI;
	double dx,dy,dz,dxo,dyo,dzo,rr,r,C1,C2,C3,C4,term1,term2,theta,incr,rxi,ryi,rxj,ryj,esps;
	double Dtr[3],M1[9],Dt[9];
	theta = atan(L1[1]/L1[0]);
	treal = t%tp;
	if(treal <= tp/2){
	    // rce = 1.0*1.5*L2yp;
	    rce = 1.0*L2yp;
	    M_HI = 3.5;
	    nmax = 1;
	}
	else{
	    // rce = 2.0*1.5*L2yp;
	    rce = 2.0*L2yp;
	    M_HI = 3.5;
	    nmax = 2;
	}
	// nmax = ceil(rce/(1.5*L2yp));
    alpha = M_HI/rce;
    rkcut = 2.0*M_HI*M_HI/rce;
    rkkcut = rkcut*rkcut;
    rtpi = sqrt(M_PI);
    alpha2 = alpha*alpha; alpha3 = alpha2*alpha; alpha4 = alpha2*alpha2; alpha5 = alpha4*alpha; alpha7 = alpha4*alpha2*alpha;
    selfconst = 1-6/rtpi*alpha + 40/(3*rtpi)*alpha3;
    kxmax = ceil(rkcut*detL/(M_PI*2.0*L*L2[1]));
    kymax = ceil(rkcut*detL/(M_PI*2.0*L*L1[0]));
    kzmax = ceil(rkcut*L/(M_PI*2.0));
    int *kyupper = calloc(kxmax+1,sizeof(int));
    int *kylower = calloc(kxmax+1,sizeof(int));
    int *kzupper = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
    int *kzlower = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
    kcoeff = M_PI*2.0/detL;
    // Determine range of wave vectors that need to be included
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
                if(rkk < rkkcut && zlower==0){
                    kzlower[(2*kymax+1)*kx + ky + kymax] = kz;
                    zlower = 1;
                }
                if(rkk > rkkcut && zlower==1 && zupper==0){
                    kzupper[(2*kymax+1)*kx + ky + kymax] = kz - 1;
                    zupper = 1;
                }
            }

            if(rkxy < rkkcut && ylower==0){
                kylower[kx] = ky;
                ylower = 1;
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
        kxmax = kx;
        kx++;
    }
    kcount = 0;
    kxmax++;
    for(kx=0;kx<kxmax+1;kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
                kcount++;
            }
        }
    }
    // Compute M2 coefficients for each wave vectors
    double *M2 = calloc(9*kcount,sizeof(double));
    kcount = 0;
    for(kx = 0; kx < kxmax+1; kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
            rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
            rkxy = rkx*rkx + rky*rky;
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
                rkz = kz*M_PI*2.0/L;
                rkk = rkxy + rkz*rkz;
                m2 = 2.0*(1.0-rkk/3.0)*(1.0+rkk/(4.0*alpha2)+rkk*rkk/(8.0*alpha4))*6.0*M_PI/rkk*exp(-rkk/(4.0*alpha2));
                start = 9*kcount;
                M2[start] = m2*(1-rkx*rkx/rkk)/box_volume;
                M2[start+1] = m2*(-rkx*rky/rkk)/box_volume;
                M2[start+2] = m2*(-rkx*rkz/rkk)/box_volume;
                M2[start+3] = m2*(-rky*rkx/rkk)/box_volume;
                M2[start+4] = m2*(1-rky*rky/rkk)/box_volume;
                M2[start+5] = m2*(-rky*rkz/rkk)/box_volume;
                M2[start+6] = m2*(-rkz*rkx/rkk)/box_volume;
                M2[start+7] = m2*(-rkz*rky/rkk)/box_volume;
                M2[start+8] = m2*(1-rkz*rkz/rkk)/box_volume;
                kcount++;
            }
        }
    }
    // Divide by 2 for x = 0 to account for symmetry in x. Set 0,0,0 equal to 0
    kcount = 0;
    for(ky = kylower[0]; ky < kyupper[0] + 1; ky++){
        for(kz = kzlower[(2*kymax + 1)*0 + ky + kymax]; kz < kzupper[(2*kymax + 1)*0 + ky + kymax] + 1; kz++){
            start = 9*kcount;
            for(i=0;i<9;i++){
                M2[start + i] /= 2.0;
            }
            if(ky==0 && kz==0){
                for(i=0;i<9;i++){
                    M2[start+i] = 0.0;
                }
            }
            kcount++;
        }
    }
    // Compute self-mobility components
    Dtr[0] = selfconst; Dtr[1] = selfconst; Dtr[2] = selfconst;
    kcount = 0;
    for(kx = 0; kx < kxmax+1; kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
            	// printf("%d %d %d %d\n",kx,ky,kz,kcount);
                start = 9*kcount;
                Dtr[0] += M2[start];
                Dtr[1] += M2[start+4];
                Dtr[2] += M2[start+8];
                kcount++;
            }
        }
    }
    for(nx = -nmax; nx < nmax + 1; nx++){
        for(ny = -nmax; ny < nmax + 1; ny++){
            for(nz = -nmax; nz < nmax + 1; nz++){
                if(nx!= 0 && ny != 0 && nz != 0){
                    // dx = nx*L1xp + ny*L2xp; dy = ny*L2yp; dz = nz*L;
                    dx = nx*L1[0] + ny*L2[0]; dy = nx*L1[1] + ny*L2[1]; dz = nz*L;
                    rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
                    C1 = 0.75/r+0.5/(r*rr);
                    C2 = 4.0*alpha7*rr*rr + 3.0*alpha3*rr - 20.0*alpha5*rr - 4.5*alpha + 14.0*alpha3 + alpha/rr;
                    C3 = 0.75/r-1.5/(r*rr);
                    C4 = -4.0*alpha7*rr*rr - 3.0*alpha3*rr + 16.0*alpha5*rr + 1.5*alpha - 2.0*alpha3 - 3.0*alpha/rr;
                    term1 = C1*erfc(alpha*r) + C2*exp(-alpha2*rr)/rtpi;
                    term2 = C3*erfc(alpha*r) + C4*exp(-alpha2*rr)/rtpi;
                    Dtr[0] += term1 + term2*dx*dx/rr;
                    Dtr[1] += term1 + term2*dy*dy/rr;
                    Dtr[2] += term1 + term2*dy*dy/rr;
                }
            }
        }
    }
	// double start_time = omp_get_wtime();
	// omp_set_num_threads(num_threads);

	//Compute off-chain mobility
	#pragma omp parallel
	{
    	int i,j,kcount,start,a,kx,ky,kz,offset,nx,ny,nz,bij,bji;
    	double dx,dy,dz,dxo,dyo,dzo,rr,r,C1,C2,C3,C4,term1,term2,rkx,rky,rkz,rkxy,rkk,coskr,epsp,incr,rxi,ryi,rxj,ryj;
    	double M1[9],Dt[9];
		kcoeff = M_PI*2.0/detL;
		epsp = 0.0;
		#pragma omp for schedule(dynamic) // schedule(static,1)
			for(i=0;i<N;i++){
				// rxi = rx[i]*cos(theta) - ry[i]*sin(theta);
				// ryi = rx[i]*sin(theta) + ry[i]*cos(theta);
				for(j=((int)(i/Nb)+1)*Nb;j<N;j++){
					// rxj = rx[j]*cos(theta) - ry[j]*sin(theta);
					// ryj = rx[j]*sin(theta) + ry[j]*cos(theta);
				    memset(M1,0.0,9*sizeof(double));
				    memset(Dt,0.0,9*sizeof(double));
					// bij = binij(i,j);
					// bji = binij(j,i);
					// #pragma omp atomic
					// 	count[bij]++;
					// 	count[bji]++;
		    		// dxo = rxj - rxi;
		    		// dyo = ryj - ryi;
					dxo = rx[j] - rx[i];
					dyo = ry[j] - ry[i];
		    		dzo = rz[j] - rz[i];
				    for(nx = -nmax; nx < nmax + 1; nx++){
				        for(ny = -nmax; ny < nmax + 1; ny++){
				            for(nz = -nmax; nz < nmax + 1; nz++){
				                dx = dxo + nx*L1[0] + ny*L2[0]; dy = dyo + nx*L1[1] + ny*L2[1]; dz = dzo + nz*L;
				                rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
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
					// Loop over wavevectors for Drecip
					kcount = 0;
				    for(kx = 0; kx < kxmax + 1; kx++){
				        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
				            rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
				            rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
				            rkxy = rkx*rkx + rky*rky;
				            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
				                rkz = kz*M_PI*2.0/L;
				                rkk = rkxy + rkz*rkz;
				                start = 9*kcount;
			                    // coskr = cos(rkx*dx + rky*dy + rkz*dz);
			                    coskr = cos(rkx*dxo + rky*dyo + rkz*dzo);
			                    Dt[0] += M2[start]*coskr;
			                    Dt[1] += M2[start+1]*coskr;
			                    Dt[2] += M2[start+2]*coskr;
			                    Dt[3] += M2[start+3]*coskr;
			                    Dt[4] += M2[start+4]*coskr;
			                    Dt[5] += M2[start+5]*coskr;
			                    Dt[6] += M2[start+6]*coskr;
			                    Dt[7] += M2[start+7]*coskr;
			                    Dt[8] += M2[start+8]*coskr;
			                    kcount++;
				            }
				        }
				    }
					// printf("%lf\n",Dt[0]);
					for(a=0;a<9;a++){
						incr = Dt[a]/Dtr[(int)(floor(a/3))];
						// printf("%d %d %lf %lf %lf %lf %lf %lf\n",i,j,NID.x,NID.y,NID.z,Dt[a],Dii[(int)(floor(a/3))],incr);
						epsp += 2.0*incr;
						#pragma omp atomic
							// Drun[9*bij + a] += Dt[a];
							// Drun[9*bji + a] += Dt[a];
							rowsum[3*i + (int)(floor(a/3))] += incr*incr;
							rowsum[3*j + (int)(floor(a/3))] += incr*incr;
					}
				}
			}
		#pragma omp critical
		{
			eps += epsp;
		}
	}
	for(ii=0;ii<Nc;ii++){
		selfcount++;
		for(i=0;i<Nb;i++){
			// rxi = rx[ii*Nb+i]*cos(theta) - ry[ii*Nb+i]*sin(theta);
			// ryi = rx[ii*Nb+i]*sin(theta) + ry[ii*Nb+i]*cos(theta);
			for(j=i+1;j<Nb;j++){
				// rxj = rx[ii*Nb+j]*cos(theta) - ry[ii*Nb+j]*sin(theta);
				// ryj = rx[ii*Nb+j]*sin(theta) + ry[ii*Nb+j]*cos(theta);
				memset(M1,0.0,9*sizeof(double));
				memset(Dt,0.0,9*sizeof(double));
				// dxo = rxj - rxi;
				// dyo = ryj - ryi;
				dxo = rx[ii*Nb+j] - rx[ii*Nb+i];
				dyo = ry[ii*Nb+j] - ry[ii*Nb+i];
				dzo = rz[ii*Nb+j] - rz[ii*Nb+i];
			    for(nx = -nmax; nx < nmax + 1; nx++){
			        for(ny = -nmax; ny < nmax + 1; ny++){
			            for(nz = -nmax; nz < nmax + 1; nz++){
			                dx = dxo + nx*L1[0] + ny*L2[0]; dy = dyo + nx*L1[1] + ny*L2[1]; dz = dzo + nz*L;
			                rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
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
				// printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",ii*Nb+i,ii*Nb+j,Dt[0],dxo,dyo,dzo,rx[ii*Nb+i],ry[ii*Nb+i],rz[ii*Nb+i],rx[ii*Nb+j],ry[ii*Nb+j],rz[ii*Nb+j]);
				kcount = 0;
			    for(kx = 0; kx < kxmax + 1; kx++){
			        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
			            rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
			            rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
			            rkxy = rkx*rkx + rky*rky;
			            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
			                rkz = kz*M_PI*2.0/L;
			                rkk = rkxy + rkz*rkz;
			                start = 9*kcount;
		                    // coskr = cos(rkx*dx + rky*dy + rkz*dz);
		                    coskr = cos(rkx*dxo + rky*dyo + rkz*dzo);
		                    Dt[0] += M2[start]*coskr;
		                    Dt[1] += M2[start+1]*coskr;
		                    Dt[2] += M2[start+2]*coskr;
		                    Dt[3] += M2[start+3]*coskr;
		                    Dt[4] += M2[start+4]*coskr;
		                    Dt[5] += M2[start+5]*coskr;
		                    Dt[6] += M2[start+6]*coskr;
		                    Dt[7] += M2[start+7]*coskr;
		                    Dt[8] += M2[start+8]*coskr;
		                    kcount++;
			            }
			        }
			    }
				startij = 9*Nb*i+9*j;
				startji = 9*Nb*j+9*i;
				// printf("%d %d %lf\n",ii*Nb+i,ii*Nb+j,Dt[0]);
				// exit(1);
				for(a=0;a<9;a++){
					incr = Dt[a]/Dtr[(int)(floor(a/3))];
					eps += 2.0*incr;
					rowsum[3*Nb*ii + 3*i + (int)(floor(a/3))] += incr*incr;
					rowsum[3*Nb*ii + 3*j + (int)(floor(a/3))] += incr*incr;
				}
			}
		}
	}
    free(kyupper);
    free(kylower);
    free(kzupper);
    free(kzlower);
    free(M2);
}
void CATEA(){
	int i,j,k,ext_bin,strain_bin;
	int ROW = 3*N;
	double btemp;
	// calcExt();
	// ext_bin = (int)floor(ext_avg/ext_interval);
	strain_bin = (int)floor(strain/strain_interval);
	if(strain_bin>sbins){
		strain_bin = sbins; // steady state bin
	}
	eps /= (9*N*N-3*N);
	eps2 = eps*eps;
	btemp = (1.0-sqrt(1.0-(ROW*eps2-ROW*eps)))/(ROW*eps2-ROW*eps);
	// printf("ts - %lu beta - %lf ext_avg - %lf  int - %d ext_bin - %d\n",t,btemp,ext_avg,(int)floor(ext_avg/ext_interval),ext_bin);
	//printf("ts - %lu beta - %lf strain - %lf bin - %d brun[strain_bin] - %lf gwcount[bin] - %d ext_avg - %lf\n",t,btemp,strain,strain_bin,brun[strain_bin],wcount[strain_bin],ext_avg);
	brun[strain_bin] += btemp;
	// printf("Beta calc done\n");
	// printf("%lu %d %.12e %lf %.12e %.12e %.12e\n",t,t%tp,eps,bp,Dii[0],Dii[1],Dii[2]);
	for(i=0;i<3*N;i++){
		Crun[strain_bin][i%(3*Nb)] += sqrt(1.0/(1.0+rowsum[i]*btemp*btemp));
		// printf("%d %d %lf %lf\n",i,i%(3*Nb),sqrt(1.0/(1.0+rowsum[i]*btemp*btemp)),Crun[ext_bin][i%(3*Nb)]);
		// Ctemp[i] = sqrt(1.0/(1.0+rowsum[i]*bp*bp));
	}
	// for(i=0;i<Nb;++i){
	// 	printf("%d %lf %lf %lf\n",i,Crun[ext_bin][3*i],Crun[ext_bin][3*i+1],Crun[ext_bin][3*i+2]);
	// }
	// exit(1);
	// printf("C calc done\n");
	gwcount[strain_bin]++;
	//printf("ts - %lu beta - %lf strain - %lf bin - %d brun[strain_bin] - %lf gwcount[bin] - %d ext_avg - %lf\n",t,btemp,strain,strain_bin,brun[strain_bin],gwcount[strain_bin],ext_avg);
	//printf("Ts - %lu eps - %lf brun - %lf Crun[0] - %lf\n",t,eps,brun/gwcount,Crun[0]/gwcount);
	betafile = fopen(beta,"a");
	fprintf(betafile,"%lu %lf\n",t,btemp);
	for(i=0;i<N;++i){
		fprintf(betafile,"%d %lf %lf %lf\n",i,sqrt(1.0/(1.0+rowsum[3*i]*btemp*btemp)),sqrt(1.0/(1.0+rowsum[3*i+1]*btemp*btemp)),sqrt(1.0/(1.0+rowsum[3*i+2]*btemp*btemp)));
	}
	fclose(betafile);
	eps = 0.0; eps2 = 0.0;
	for(i=0;i<3*N;i++){
		rowsum[i] = 0.0;
	}
}
void updateLattice(unsigned long t){
	int i;
    int treal = t%tp;
	// printf("%d %lu\n",treal,t);
	// exit(1);
    //if(treal==0) printf("t=%d\n",t);
	// printf("%lf %d\n",flowrate,treal);
    for(i = 0; i<3; ++i){
        point[i][0] = point0[i][0]*exp(flowrate*treal*dt); point[i][1] = point0[i][1]*exp(-flowrate*treal*dt);
    }
    for(i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
    }
    //printf("x = %f y = %f\n", point[1][0], point[1][1]);
	double theta = atan(L1[1]/L1[0]);
	L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
	L2xp = L2[0]*cos(theta) + L2[1]*sin(theta);
	L1yp = -L1[0]*sin(theta) + L1[1]*cos(theta);
	L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
	detL = L*(-L1yp*L2xp + L1xp*L2yp);
	thetabox = atan(L2yp/L2xp);
	L2pmag = sqrt(L2xp*L2xp + L2yp*L2yp);
	// printf("%lu %lf %lf %lf %lf %lf\n",t,L1xp,L1yp,L2xp,L2yp,detL);
	// exit(1);
}
void applyPBC(){
	int i;
    double theta = atan(L1[1]/L1[0]);
    for(i = 0; i<N; ++i){
        rx[i] -= L2xp*round(ry[i]/L2yp);
        ry[i] -= L2yp*round(ry[i]/L2yp);
        rx[i] -= L1xp*round((rx[i]-ry[i]*L2xp/L2yp)/L1xp);
		// rx[i] -= L1xp*round(rx[i]/L1xp);
        rz[i] -= L*round(rz[i]/L);
    }
}
Vector3D_t getNID(int i,int j){
	Vector3D_t NID;
	double dxp,dyp;
	// printf("%lf %lf\n",L1[1],L1[0]);
    double theta = atan(L1[1]/L1[0]);
	// printf("%lf\n",theta);
    double rx1 = rx[i]; double rx2 = rx[j];
    double ry1 = ry[i]; double ry2 = ry[j];
    double rz1 = rz[i]; double rz2 = rz[j];
    double dx = rx1 - rx2; double dy = ry1 - ry2; double dz = rz1 - rz2;
    dx -= L2xp*round(dy/L2yp);
    dy -= L2yp*round(dy/L2yp);
    dx -= L1xp*round(dx/L1xp);
    dz -= L*round(dz/L);
	NID.x = dx*cos(theta) - dy*sin(theta);
	NID.y = dx*sin(theta) + dy*cos(theta);
	NID.z = dz;
	return NID;
}
Vector3D_t getNIDbin(int i,int j){
	Vector3D_t NID;
	int nx,ny,nmax;
	double dx,dy,dz,dxp,dyp,theta,dxL2,dxmin,dymin,rp,rmin;
	theta = atan(L1[1]/L1[0]);
	dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j];
	int treal = t%tp;
	if(treal <= tp/2){
	    nmax = 1;
	}
	else{
	    nmax = 2;
	}
	dz -= L*round(dz/L);
	dxmin = dx; dymin = dy;
	rmin = sqrt(dx*dx + dy*dy);
	for(nx=-nmax;nx<nmax+1;nx++){
		for(ny=-nmax;ny<nmax+1;ny++){
			dxp = dx + nx*L1xp + ny*L2xp;
			dyp = dy + ny*L2yp;
			rp = sqrt(dxp*dxp + dyp*dyp);
			if(rp < rmin){
				dxmin = dxp;
				dymin = dyp;
				rmin = rp;
			}
		}
	}
	NID.x = dxmin*cos(theta) - dymin*sin(theta);
	NID.y = dxmin*sin(theta) + dymin*cos(theta);
	NID.z = dz;
	return NID;
}
void transformCoords(){
	///Rotate coordinates based on rotation matrix//
	int i;
	double rxi,ryi;
    double theta = atan(L1[1]/L1[0]);
    for(i = 0; i<N; ++i){
        rxi = rx[i]; ryi = ry[i];
        rx[i] = rxi*cos(theta) + ryi*sin(theta);
        ry[i] = -rxi*sin(theta) + ryi*cos(theta);
    }
}
void inverseCoords(){
	///Rotate coordinates based on inverse rotation matrix//
	int i;
	double rxi,ryi;
    double theta = atan(L1[1]/L1[0]);
	for(i = 0; i<N; ++i){
        rxi = rx[i]; ryi = ry[i];
        rx[i] = rxi*cos(theta) - ryi*sin(theta);
        ry[i] = rxi*sin(theta) + ryi*cos(theta);
    }
}
int binij(int i,int j){
	int bin,bin_x,bin_y,bin_z;
	double dx,dy,dz,thetaij,phipbc,dxp,dyp,r,drproj,theta;
	// Vector3D_t NID = getNID(j,i);
	Vector3D_t NID = getNIDbin(j,i);
	bin_x = round(NID.x/bin_size_x) + nhx;
	bin_y = round(NID.y/bin_size_y) + nhy;
	bin_z = round(NID.z/bin_size_z) + nhz;
	bin = num_bins_y*num_bins_z*bin_x + num_bins_z*bin_y + bin_z;
	// printf("i - %d j - %d NID.x - %lf NID.y - %lf NID.z - %lf\n",i,j,NID.x,NID.y,NID.z);
	return bin;
}
void calcHI(){
	int a,i,j,k,ii,kx,ky,kz,offset,start,kxmax,kymax,kzmax,nmax,nx,ny,nz,nn,kcount,ylower,zlower,yupper,zupper,treal,startii,startij,startji,bij,bji,bin;
	int bin_x,bin_y,bin_z,startavg,nt;
	double alpha,rkx,rky,rkz,rkk,m2,kcoeff,rkxx,rkxy,rtpi,coskr,selfconst,alpha2,alpha3,alpha4,alpha5,alpha7,rce,rkcut,rkkcut,M_HI;
	double dx,dy,dz,dxo,dyo,dzo,rr,r,C1,C2,C3,C4,term1,term2,theta,incr,rxi,ryi,rxj,ryj,esps,xx,tf;
	double Dtr[3],M1[9],Dt[9];
	theta = atan(L1[1]/L1[0]);
	treal = t%tp;
	xx = treal/tp*nt_avg;
	if(treal <= tp/2){
	    rce = 1.0*L2yp;
	    M_HI = 4.0;
	    nmax = 1;
	}
	else{
	    rce = 2.0*L2yp;
	    M_HI = 4.0;
	    nmax = 2;
	}
	tf = (double)treal/tp;
	nt = (int)round(tf*nt_avg)%nt_avg;
    alpha = M_HI/rce;
    rkcut = 2.0*M_HI*M_HI/rce;
    rkkcut = rkcut*rkcut;
    rtpi = sqrt(M_PI);
    alpha2 = alpha*alpha; alpha3 = alpha2*alpha; alpha4 = alpha2*alpha2; alpha5 = alpha4*alpha; alpha7 = alpha4*alpha2*alpha;
    selfconst = 1-6/rtpi*alpha + 40/(3*rtpi)*alpha3;
    kxmax = ceil(rkcut*detL/(M_PI*2.0*L*L2[1]));
    kymax = ceil(rkcut*detL/(M_PI*2.0*L*L1[0]));
    kzmax = ceil(rkcut*L/(M_PI*2.0));
    int *kyupper = calloc(kxmax+1,sizeof(int));
    int *kylower = calloc(kxmax+1,sizeof(int));
    int *kzupper = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
    int *kzlower = calloc((kxmax+1)*(2*kymax+1),sizeof(int));
    kcoeff = M_PI*2.0/detL;
    // Determine range of wave vectors that need to be included
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
                if(rkk < rkkcut && zlower==0){
                    kzlower[(2*kymax+1)*kx + ky + kymax] = kz;
                    zlower = 1;
                }
                if(rkk > rkkcut && zlower==1 && zupper==0){
                    kzupper[(2*kymax+1)*kx + ky + kymax] = kz - 1;
                    zupper = 1;
                }
            }

            if(rkxy < rkkcut && ylower==0){
                kylower[kx] = ky;
                ylower = 1;
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
        kxmax = kx;
        kx++;
    }
    kcount = 0;
    kxmax++;
    for(kx=0;kx<kxmax+1;kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
                kcount++;
            }
        }
    }
    // Compute M2 coefficients for each wave vectors
    double *M2 = calloc(9*kcount,sizeof(double));
    kcount = 0;
    for(kx = 0; kx < kxmax+1; kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
            rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
            rkxy = rkx*rkx + rky*rky;
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
                rkz = kz*M_PI*2.0/L;
                rkk = rkxy + rkz*rkz;
                m2 = 2.0*(1.0-rkk/3.0)*(1.0+rkk/(4.0*alpha2)+rkk*rkk/(8.0*alpha4))*6.0*M_PI/rkk*exp(-rkk/(4.0*alpha2));
                start = 9*kcount;
                M2[start] = m2*(1-rkx*rkx/rkk)/box_volume;
                M2[start+1] = m2*(-rkx*rky/rkk)/box_volume;
                M2[start+2] = m2*(-rkx*rkz/rkk)/box_volume;
                M2[start+3] = m2*(-rky*rkx/rkk)/box_volume;
                M2[start+4] = m2*(1-rky*rky/rkk)/box_volume;
                M2[start+5] = m2*(-rky*rkz/rkk)/box_volume;
                M2[start+6] = m2*(-rkz*rkx/rkk)/box_volume;
                M2[start+7] = m2*(-rkz*rky/rkk)/box_volume;
                M2[start+8] = m2*(1-rkz*rkz/rkk)/box_volume;
                kcount++;
            }
        }
    }
    // Divide by 2 for x = 0 to account for symmetry in x. Set 0,0,0 equal to 0
    kcount = 0;
    for(ky = kylower[0]; ky < kyupper[0] + 1; ky++){
        for(kz = kzlower[(2*kymax + 1)*0 + ky + kymax]; kz < kzupper[(2*kymax + 1)*0 + ky + kymax] + 1; kz++){
            start = 9*kcount;
            for(i=0;i<9;i++){
                M2[start + i] /= 2.0;
            }
            if(ky==0 && kz==0){
                for(i=0;i<9;i++){
                    M2[start+i] = 0.0;
                }
            }
            kcount++;
        }
    }
    Dtr[0] = selfconst; Dtr[1] = selfconst; Dtr[2] = selfconst;
    kcount = 0;
    for(kx = 0; kx < kxmax+1; kx++){
        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
            	// printf("%d %d %d %d\n",kx,ky,kz,count);
                start = 9*kcount;
                Dtr[0] += M2[start];
                Dtr[1] += M2[start+4];
                Dtr[2] += M2[start+8];
                kcount++;
            }
        }
    }
    for(nx = -nmax; nx < nmax + 1; nx++){
        for(ny = -nmax; ny < nmax + 1; ny++){
            for(nz = -nmax; nz < nmax + 1; nz++){
                if(nx!= 0 && ny != 0 && nz != 0){
                    // dx = nx*L1xp + ny*L2xp; dy = ny*L2yp; dz = nz*L;
                    dx = nx*L1[0] + ny*L2[0]; dy = nx*L1[1] + ny*L2[1]; dz = nz*L;
                    rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
                    C1 = 0.75/r+0.5/(r*rr);
                    C2 = 4.0*alpha7*rr*rr + 3.0*alpha3*rr - 20.0*alpha5*rr - 4.5*alpha + 14.0*alpha3 + alpha/rr;
                    C3 = 0.75/r-1.5/(r*rr);
                    C4 = -4.0*alpha7*rr*rr - 3.0*alpha3*rr + 16.0*alpha5*rr + 1.5*alpha - 2.0*alpha3 - 3.0*alpha/rr;
                    term1 = C1*erfc(alpha*r) + C2*exp(-alpha2*rr)/rtpi;
                    term2 = C3*erfc(alpha*r) + C4*exp(-alpha2*rr)/rtpi;
                    Dtr[0] += term1 + term2*dx*dx/rr;
                    Dtr[1] += term1 + term2*dy*dy/rr;
                    Dtr[2] += term1 + term2*dy*dy/rr;
                }
            }
        }
    }
    for(i=0;i<N;i++){
    	startij = 9*N*i + 9*i;
    	M[startij] = Dtr[0];
    	M[startij+4] = Dtr[1];
    	M[startij+8] = Dtr[2];
    }
	#pragma omp parallel
	{
		int a,i,j,k,kx,ky,kz,start,nx,ny,nz,kcount,startij,bin;
		int bin_x,bin_y,bin_z,startavg,nxmin,nymin;
		double rkx,rky,rkz,rkk,rkxx,rkxy,coskr,rxi,ryi,rxj,ryj,rtest,dxmin,dymin,dzmin,rmin,drrbmin;
		double dx,dy,dz,dxo,dyo,dzo,rr,r,C1,C2,C3,C4,term1,term2;
		double Dtr[3],M1[9],Dt[9];
		#pragma omp for schedule(dynamic)
			for(i=0;i<N;i++){
				rxi = rx[i]*cos(theta) - ry[i]*sin(theta);
				ryi = rx[i]*sin(theta) + ry[i]*cos(theta);
				// for(j=(((int)(i/Nb)+1)*Nb);j<N;j++){
				for(j=i+1;j<N;j++){
					startij = 9*N*i+9*j;
					rxj = rx[j]*cos(theta) - ry[j]*sin(theta);
					ryj = rx[j]*sin(theta) + ry[j]*cos(theta);
					Vector3D_t drbmin = getNIDbin(j,i);
					// Vector3D_t drbmin = getNID(j,i);
					drrbmin = sqrt(drbmin.x*drbmin.x + drbmin.y*drbmin.y + drbmin.z*drbmin.z);
					// if(drrbmin < 10.0){
					// if(drbmin.y > (L1[1]+L2[1])/2.0){
					if(drrbmin < 12.0){
						memset(M1,0.0,9*sizeof(double));
						memset(Dt,0.0,9*sizeof(double));
						// dxo = rx[j] - rx[i];
						// dyo = ry[j] - ry[i];
						dxo = rxj - rxi;
						dyo = ryj - ryi;
						dzo = rz[j] - rz[i];
					    for(nx = -nmax; nx < nmax + 1; nx++){
					        for(ny = -nmax; ny < nmax + 1; ny++){
					            for(nz = -nmax; nz < nmax + 1; nz++){
					                dx = dxo + nx*L1[0] + ny*L2[0]; dy = dyo + nx*L1[1] + ny*L2[1]; dz = dzo + nz*L;
					                rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
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
						kcount = 0;
					    for(kx = 0; kx < kxmax + 1; kx++){
					        for(ky = kylower[kx]; ky < kyupper[kx] + 1; ky++){
					            rkx = kcoeff*L*(kx*L2[1] - ky*L1[1]);
					            rky = kcoeff*L*(-kx*L2[0] + ky*L1[0]);
					            rkxy = rkx*rkx + rky*rky;
					            for(kz = kzlower[(2*kymax + 1)*kx + ky + kymax]; kz < kzupper[(2*kymax + 1)*kx + ky + kymax] + 1; kz++){
					                rkz = kz*M_PI*2.0/L;
					                rkk = rkxy + rkz*rkz;
					                start = 9*kcount;
				                    // coskr = cos(rkx*dx + rky*dy + rkz*dz);
				                    coskr = cos(rkx*dxo + rky*dyo + rkz*dzo);
				                    Dt[0] += M2[start]*coskr;
				                    Dt[1] += M2[start+1]*coskr;
				                    Dt[2] += M2[start+2]*coskr;
				                    Dt[3] += M2[start+3]*coskr;
				                    Dt[4] += M2[start+4]*coskr;
				                    Dt[5] += M2[start+5]*coskr;
				                    Dt[6] += M2[start+6]*coskr;
				                    Dt[7] += M2[start+7]*coskr;
				                    Dt[8] += M2[start+8]*coskr;
				                    kcount++;
					            }
					        }
					    }
						memcpy(&M[startij],&Dt[0],9*sizeof(double));
					}
					else{
						bin = binij(i,j);
						memcpy(&M[startij],&D[nt][9*bin],9*sizeof(double));
					}
				}
			}
	}
	// printf("Ts - %lu Perc error - %lf Diff - %lf Denom - %lf\n",t,error/denomin*100.0,error,denomin);
	// printf("bins - %d conts - %d\n",bins,conts);
	// exit(1);
    free(kyupper);
    free(kylower);
    free(kzupper);
    free(kzlower);
    free(M2);
}
void printTEA(){
	int i,j;
	DCfile = fopen(decomp,"w");
	fprintf(DCfile,"%lu\n",t);
	for(i=0;i<sbins+1;++i){
		if(gwcount[i]>0){
			fprintf(DCfile,"%d %d %.12e\n",i,gwcount[i],brun[i]/gwcount[i]);
			for(j=0;j<Nb;j++){
				fprintf(DCfile,"%d %.12e %.12e %.12e\n",j,Crun[i][3*j]/(Nc*gwcount[i]),Crun[i][3*j+1]/(Nc*gwcount[i]),Crun[i][3*j+2]/(Nc*gwcount[i]));
			}
		}
		else{
			fprintf(DCfile,"%d %d %.12e\n",i,0,0.0);
			for(j=0;j<Nb;++j){
				fprintf(DCfile,"%d %.12e %.12e %.12e\n",j,0.0,0.0,0.0);
			}
		}
	}
	fclose(DCfile);
}
void resetAverage(){
	int i,j,nt,k,atemp,btemp,ctemp,dtemp,l,m,start;
	// bp = brun/gwcount;
	for(i=0;i<sbins+1;++i){
		bp[i] = brun[i]/gwcount[i];
		brun[i] = 0.0;
		for(j=0;j<3*Nb;j++){
			C[i][j] = Crun[i][j]/(Nc*gwcount[i]);
			Crun[i][j] = 0.0;
		}
		gwcount[i] = 0;
	}
	if(itercount==0){
		gridfile = fopen(grid,"r");
		if(!gridfile)
		{
				printf("Error - grid file not found, make sure input parameters match %s\n",grid);
				exit(1);
		}
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
	}
}
void printVisc(){
	// viscfile = fopen(visc,"a");
	// fprintf(viscfile,"%lu %lf\n",t,tau12);
	// fclose(viscfile);
	// tau12 = 0.0;
	int i;
	if(t>tstart && t%printperiod==0){
		viscfile = fopen(visc,"a");
		fprintf(viscfile,"%lu %lf %lf\n",t,tauxx/((double)printperiod),tauyy/((double)printperiod));
		fclose(viscfile);
		tauxx = 0.0;
		tauyy = 0.0;
	}
}
void printTiming(){
	if(t%printperiod==0){
		double elapsed,tpts;
		elapsed = omp_get_wtime() - total_time;
		tpts = elapsed/(t+1-tstart);
		timfile = fopen(tim,"w");
		fprintf(timfile,"Time steps - %lu\nWall time - %lf\nTime per time step - %lf\nBin time per time step - %lf\nMult time per time step - %lf\n",t,elapsed,tpts,bin_time/((double)t-(double)tstart),mult_time/((double)t-(double)tstart));
		fclose(timfile);
	}
}
void coulombForce(){
	double dx,dy,dz,rr,r,coeff,r6,ratio,Fc,box_length;
	int i,j,k,l;

	// for(i=0;i<N;i++){
	// 	printf("Charge[i]: %d \n", Charge[i]) ;
	// }
	// printf("lambda_d: %lf \n", lambda_d);
	kappa_d = 1/lambda_d ;

	for (i=0; i<N; ++i){
		for (j=i+1; j<N; ++j){
			if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
			{
				Vector3D_t NID = getNID(i,j);
				rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
				r = sqrt(rr) ;

				Fc = lambda_b*Charge[i]*Charge[j]*exp(-kappa_d*r)*(1+r*kappa_d)/rr;
				//printf("%lf \n",Fc);

				//E_initial += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Initial Energy

				fx[i] += Fc*NID.x/r;
				fy[i] += Fc*NID.y/r;
				fz[i] += Fc*NID.z/r;
				fx[j] -= Fc*NID.x/r;
				fy[j] -= Fc*NID.y/r;
				fz[j] -= Fc*NID.z/r;


			}
		}
	}
}
void chargeHop(){

	
	int i, j, k,c,l, test, test2, cnumber, cnumber2, count, NCharges,num_of_particles,cindex;
	double E_initial, E_final, dx, dy, dz, r, rr , mm, dice_roll, rbix,rbiy,rbiz;
	double rxtempi, rxtempf,rytempi,rytempf,rztempi,rztempf;
	rxtempi = 0, rxtempf = 0, rytempi = 0, rytempf = 0, rztempi = 0, rztempf = 0;

	NCharges = Ncharges;
	num_of_particles = Nb*Nc ;
	int num_of_p_particles = Nb*Nc ;
	double theta = atan(L1[1]/L1[0]);
    //if(t==0){printf("Warning - nonadjacent hopping turned off completely\n");}

	if (t%MStep==0) // Do Monte Carlo Move every MSteps
	{
		for (j=0;j<NCharges;++j)
		{
			for (i=0; i<NCharges;++i)
			{
				if (Charge_track[j]==Charge_indices[i] && i==j)
				{i_values[i]=j;}
			}
		}

		for(i=0;i<NCharges;++i)
		{	Charge_indices[i]=Charge_track[i];}

        //Calculate Initial Coulomb Energy of System//
        E_initial=0;
        if(Ncharges>1){
            for (i=0; i<num_of_particles; ++i)
            {
                for (j=i+1; j<num_of_particles; ++j)
                {
                    if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
                    {
                        Vector3D_t NID = getNID(i,j);
                        rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
                        r = sqrt(rr) ;
                        E_initial += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Final Energy
                    }
                }
            }
        }

	for (k=0; k<NCharges; ++k) // Do MC_Move for each charge
	{

		cindex = Charge_track[k] ;
		rbix = rx[cindex] ;
		rbiy = ry[cindex] ;
		rbiz = rz[cindex] ;
		for(i=0;i<num_of_p_particles;++i)
		{
			Charge[i]=0;
		}

		for(i=0;i<num_of_p_particles;++i)
		{
			for(j=0;j<NCharges;++j)
			{
				if(i==Charge_track[j])
				Charge[i]=1;

			}
		}

		//printf("mm=%f\n",mm);
		mm=ran1(idum);
		cnumber=Charge_indices[k]; //find the next charge on the shuffled list
		test=0;// test for am I close?
		test2=0;

		for(j = 0; j<num_of_p_particles; ++j) //for loop to check distance between cnumber and every other particle
					{
						if (j!=cnumber && j!=cnumber+1 && j!=cnumber-1)
						{
							Vector3D_t NID = getNID(cnumber,j);
							rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
								if(rr<4.41)
								{
					if (mm<.33) //Move to the right
					{
						cnumber2=cnumber+1; //right hand neighbor
						//printf("t=%d RR cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
						for(i=0;i<NCharges;++i)
						{

							if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1 || Charge_indices[i]==cnumber2) // don't check if right hand neighbor is charged
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("RR Skip\n");
							continue;
						}
						for(i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber )
							{
								Charge_indices[i]=cnumber2;
								//printf("RR Swap\n");
								test2=0;
							}
						}

						test=1;//possibility of hopping to NOT neighbor occurred

					}
					if(mm>0.66) //move to the left - if random # is less than .5
					{
						cnumber2=cnumber-1; //left hand neighbor
						for (i=0;i<NCharges;++i)
						{
							if(cnumber2<0 || cnumber==0 || cnumber%Nb==0 || Charge_indices[i]==cnumber2)// don't attempt to go left at bead 0
							{
								test2=1;

							}
						}
						if (test2!=0)
						{
							//printf("LL Skip\n");
								continue;
						}
						else{
							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{	Charge_indices[i]=cnumber2;
									//printf("LL Swap\n");
									test2=0;
								}
							}
							test=1;//possibility of hopping to NOT neighbor occurred
						}
					}
					else if(mm<=0.66 && mm >=0.33) ///swap with j particle that's close to cnumber
					{
						cnumber2=j;
						//printf("t=%d P cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
						for(i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber2)
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("Pop Skip\n");
								continue;
						}
						else{
							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber)
								{
									//printf("Pop\n");  //Pop=non-adjacent hop
									Charge_indices[i]=cnumber2;
									test2=0;
								}
							}
							test=2;
						}
					}
					break;
				}
			}
		}

		//////////// Perform swap if NOT close enough, only to a neighbor bead///////////////
			if(test==0){
					if (mm>.5) //Move to the right
					{

					cnumber2=cnumber+1; //right hand neighbor
					//printf("t=%d R cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

						for(i=0;i<NCharges;++i)
						{

							if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1 || Charge_indices[i]==cnumber2) // don't check if right hand neighbor is charged
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("R Skip\n");
								continue;
						}
						else{
							for (i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
									{
										Charge_indices[i]=cnumber2;
										//printf("R Swap\n");
										test2=0;
									}
							}
						}
					}
					else //move to the left - if random # is less than .5
					{

						cnumber2=cnumber-1; //left hand neighbor
						//printf("t=%d L cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
						for (i=0;i<NCharges;++i)
						{
							if(cnumber2<0 || cnumber==0 || cnumber%Nb==0 || Charge_indices[i]==cnumber2)// don't attempt to go left at bead 0
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("L Skip\n");
								continue;
						}
						else{
							for (i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber)
								{
									Charge_indices[i]=cnumber2;
									//printf("L Swap\n");
									test2=0;
									test=0;

								}
							}
						}
					}
				}

				// Update the charge positions to test new energy
				for(i=0;i<num_of_p_particles;++i)
				{
					Charge[i]=0;
				}

				for(i=0;i<num_of_p_particles;++i)
				{
					for(j=0;j<NCharges;++j)
					{
						if(i==Charge_indices[j])
						Charge[i]=1;

					}
				}

				//Calculate New Energy if Multiple Charges Exist
				E_final=0;
				if(Ncharges>1){
					for (i=0; i<num_of_particles; ++i)
					{
						for (j=i+1; j<num_of_particles; ++j)
								{
									if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
									{
										Vector3D_t NID = getNID(i,j);
										rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;
										r = sqrt(rr) ;
										E_final += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Final Energy
									}
								}
					}
				}
				// printf("this is Efinal=");
				// printf("%lf\n",E_final);
				 //printf("this delta E=(%lf)\n", (E_final-E_initial));
				if (test==0 || test==1) // Perform MC but under neighbor conditions - barrier
				{
				//Actual Monte Carlo Test
					if (1<=exp((-barrier-0.5*(E_final - E_initial))))
					{
						count=0;
										for(i = 0; i<NCharges; ++i)
										{
												if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
												{ 	count+=1;
													Charge_old[i]=cnumber;
													Charge_indices[i]=cnumber2;
														Charge_track[i]=cnumber2;
										// 				Output2 = fopen(str2, "a");
										// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
										// // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
										// fclose(Output2);
												}
										}
										counter1++ ;
										cindex = Charge_track[k] ;
										//printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
										dx = rx[cindex] - rbix ;
										dy = ry[cindex] - rbiy ;
										dz = rz[cindex] - rbiz ;
										//printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
										dx -= L2xp*round(dy/L2yp);
										dy -= L2yp*round(dy/L2yp);
										dx -= L1xp*round(dx/L1xp);
										dz -= L*round(dz/L);

										dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
										dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
										dzt[k*printprops+tcount] += dz;
					}
					else
					{
							dice_roll=ran1(idum);
							if (dice_roll<=exp((-barrier-0.5*(E_final - E_initial)))) // if greater, accept swap
							{	count=0;
								for(i = 0; i<NCharges; ++i)
								{
										if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
										{
											count+=1;
											Charge_old[i]=cnumber;
											Charge_indices[i]=cnumber2;
												Charge_track[i]=cnumber2;
							// 				Output2 = fopen(str2, "a");
							// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
							// fclose(Output2);
							// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
										}
								}
								counter2++ ;
								cindex = Charge_track[k] ;
								//printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
								dx = rx[cindex] - rbix ;
								dy = ry[cindex] - rbiy ;
								dz = rz[cindex] - rbiz ;
								//printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
								dx -= L2xp*round(dy/L2yp);
								dy -= L2yp*round(dy/L2yp);
								dx -= L1xp*round(dx/L1xp);
								dz -= L*round(dz/L);

								dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
								dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
								dzt[k*printprops+tcount] += dz;
							}
							else
							{	count=0;
								for(i=0;i<NCharges; ++i)
								{
									if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
									{	count+=1;
										// printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
										Charge_indices[i]=cnumber;
										Charge_old[i]=cnumber;
										}
								}
								counter3++ ;
							}
						}
					}

					///////// If there's a possibility of a bead (not neighbor) close by to hop a charge
					else if(test==2)
					{
						if (1<=exp((-barrier2-0.5*(E_final - E_initial))))
						{
							count=0;
							for(i = 0; i<NCharges; ++i)
							{
									if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
									{ 	count+=1;
										Charge_old[i]=cnumber;
										Charge_indices[i]=cnumber2;
											Charge_track[i]=cnumber2;
							// 				Output2 = fopen(str2, "a");
							// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
							// // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
							// fclose(Output2);
									}
							}
							counter4++ ;
							cindex = Charge_track[k] ;
							//printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
							dx = rx[cindex] - rbix ;
							dy = ry[cindex] - rbiy ;
							dz = rz[cindex] - rbiz ;
							//printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
							dx -= L2xp*round(dy/L2yp);
							dy -= L2yp*round(dy/L2yp);
							dx -= L1xp*round(dx/L1xp);
							dz -= L*round(dz/L);

							dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
							dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
							dzt[k*printprops+tcount] += dz;
						}
						else{
							dice_roll=ran1(idum);
							if (dice_roll<=exp((-barrier2-0.5*(E_final - E_initial)))) // if greater, accept swap
							{	count=0;
								for(i = 0; i<NCharges; ++i)
								{
										if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
										{
											count+=1;
											Charge_old[i]=cnumber;
											Charge_indices[i]=cnumber2;
												Charge_track[i]=cnumber2;
							// 				Output2 = fopen(str2, "a");
							// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
							// fclose(Output2);
							// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
										}
								}
								counter5++ ;
								cindex = Charge_track[k] ;
								//printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
								dx = rx[cindex] - rbix ;
								dy = ry[cindex] - rbiy ;
								dz = rz[cindex] - rbiz ;
								//printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
								dx -= L2xp*round(dy/L2yp);
								dy -= L2yp*round(dy/L2yp);
								dx -= L1xp*round(dx/L1xp);
								dz -= L*round(dz/L);

								dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
								dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
								dzt[k*printprops+tcount] += dz;
							}
							else
							{	count=0;
								for(i=0;i<NCharges; ++i)
								{
									if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
									{	count+=1;
										// printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
										Charge_indices[i]=cnumber;
										Charge_old[i]=cnumber;
										}
								}
								counter6++ ;
							}
						}
					}

		}//end of k loop
	}//end of MC step

	for(i=0;i<num_of_p_particles;++i)
	{
		Charge[i]=0;
	}

	for(i=0;i<num_of_p_particles;++i)
	{
		for(j=0;j<NCharges;++j)
		{
			if(i==Charge_track[j])
			Charge[i]=1;
		}
	}

}//end of function
void chargeHop_nonint(){

	int i, j, k,c,l, test, test2, cnumber, cnumber2, count, NCharges,num_of_particles,cindex;
	double E_initial, E_final, dx, dy, dz, r, rr , mm, dice_roll;
	double rxtempi, rxtempf,rytempi,rytempf,rztempi,rztempf,rbix,rbiy,rbiz;
	rxtempi = 0, rxtempf = 0, rytempi = 0, rytempf = 0, rztempi = 0, rztempf = 0;
	for(i=0;i<Ncharges;i++){
		check_ch[3*i]=0;
		check_ch[3*i+1]=0;
		check_ch[3*i+2]=0;
	}

	NCharges = Ncharges;
	num_of_particles = N ;
    double theta = atan(L1[1]/L1[0]);

	if (t%MStep==0) // Do Monte Carlo Move every MSteps
	{
		for (j=0;j<NCharges;++j)
		{
			for (i=0; i<NCharges;++i)
			{
				if (Charge_track[j]==Charge_indices[i] && i==j)
				{i_values[i]=j;}
			}
		}

		for(i=0;i<NCharges;++i)
		{	Charge_indices[i]=Charge_track[i];}

	for (k=0; k<NCharges; ++k) // Do MC_Move for each charge
	{		
            cindex = Charge_track[k] ;
            rbix = rx[cindex] ;
            rbiy = ry[cindex] ;
            rbiz = rz[cindex] ;
            
            mm=ran1(idum);
			//printf("mm=%f\n",mm);

			cnumber=Charge_indices[k]; //find the next charge on the shuffled list
			test=0;// test for am I close?
			test2=0;
			for(j = 0; j<num_of_particles; ++j) //for loop to check distance between cnumber and every other particle
						{
							if (j!=cnumber && j!=cnumber+1 && j!=cnumber-1)
							{
								Vector3D_t NID = getNID(cnumber,j);
							    rr = NID.x*NID.x + NID.y*NID.y + NID.z*NID.z;


                                if(rr<4.41)
                                {
						if (mm<.33) //Move to the right
						{
							cnumber2=cnumber+1; //right hand neighbor
							//printf("t=%d RR cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
							if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1) // don't check if right hand neighbor is charged
							{
								continue;
							}

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{
									Charge_indices[i]=cnumber2;
									//printf("RR Swap\n");
								}
							}

							test=1;//possibility of hopping to NOT neighbor occurred

						}
						if(mm>0.66) //move to the left - if random # is less than .5
						{
							cnumber2=cnumber-1; //left hand neighbor

							if(cnumber==0 || cnumber%Nb==0)// don't attempt to go left at bead 0
							{
								continue;
							}

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{	Charge_indices[i]=cnumber2;
									//printf("LL Swap\n");
								}
							}
							test=1;//possibility of hopping to NOT neighbor occurred
						}
						else if(mm<=0.66 && mm >=0.33) ///swap with j particle that's close to cnumber
						{
							cnumber2=j;
							//printf("t=%d P cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber)
								{
									//printf("Pop\n");  //Pop=non-adjacent hop
									Charge_indices[i]=cnumber2;
								}
							}
							test=2;
						}
						break;
					}
				}
			}


//////////// Perform swap if NOT close enough, only to a neighbor bead///////////////
		if(test==0)
		{
			if (mm>.5) //Move to the right
			{

			cnumber2=cnumber+1; //right hand neighbor
			//printf("t=%d R cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

				if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1) // don't check if right hand neighbor is charged
				{
					continue;
				}

				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber )
						{
							Charge_indices[i]=cnumber2;
							//printf("R Swap\n");
						}
				}
			}
			else //move to the left - if random # is less than .5
			{

				cnumber2=cnumber-1; //left hand neighbor
				//printf("t=%d L cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

				if(cnumber==0 || cnumber%Nb==0)// don't attempt to go left at bead 0
				{
					continue;
				}

				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber)
					{
						Charge_indices[i]=cnumber2;
						//printf("L Swap\n");
					}
				}
			}
		}

		if (test==0 || test==1) // Perform MC but under neighbor conditions - barrier
		{


			//Actual Monte Carlo Test
			E_final = 0;
			E_initial = 0;
			if (1<=exp((-barrier-0.5*(E_final - E_initial)))) //if energy change is greater than 1 accept change
			{	count=0;
                for(i = 0; i<NCharges; ++i)
                {
                        if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
                        { 	count+=1;
                            Charge_old[i]=cnumber;
                            Charge_indices[i]=cnumber2;
                                Charge_track[i]=cnumber2;
                // 				Output2 = fopen(str2, "a");
                // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
                // // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
                // fclose(Output2);
                        }
                }
                counter1++ ;
                cindex = Charge_track[k] ;
                //printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
                dx = rx[cindex] - rbix ;
                dy = ry[cindex] - rbiy ;
                dz = rz[cindex] - rbiz ;
                //printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
                dx -= L2xp*round(dy/L2yp);
                dy -= L2yp*round(dy/L2yp);
                dx -= L1xp*round(dx/L1xp);
                dz -= L*round(dz/L);

                dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
                dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
                dzt[k*printprops+tcount] += dz;
			}

			else
			{
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier-0.5*(E_final - E_initial)))) // if greater, accept swap
				{	count=0;
                    for(i = 0; i<NCharges; ++i)
                    {
                            if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
                            {
                                count+=1;
                                Charge_old[i]=cnumber;
                                Charge_indices[i]=cnumber2;
                                Charge_track[i]=cnumber2;
                // 				Output2 = fopen(str2, "a");
                // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
                // fclose(Output2);
                // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
                            }
                    }
                    counter2++ ;
                    cindex = Charge_track[k] ;
                    //printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
                    dx = rx[cindex] - rbix ;
                    dy = ry[cindex] - rbiy ;
                    dz = rz[cindex] - rbiz ;
                    //printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
                    dx -= L2xp*round(dy/L2yp);
                    dy -= L2yp*round(dy/L2yp);
                    dx -= L1xp*round(dx/L1xp);
                    dz -= L*round(dz/L);

                    dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
                    dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
                    dzt[k*printprops+tcount] += dz;
                }
				else
				{	count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{	count+=1;
							// printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
							Charge_indices[i]=cnumber;
							}
					}
					counter3++ ;
				}
			}
		}

		///////// If there's a possibility of a bead (not neighbor) close by to hop a charge


		else if(test==2) ///Perform MC under NOT neighbor conditions - barrier2
		{
			//Actual Monte Carlo Test
			if (1<=exp((-barrier2-0.5*(E_final - E_initial)))) //if energy change is greater than 1 accept change
			{
				count=0;
                for(i = 0; i<NCharges; ++i)
                {
                        if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
                        {	count+=1;

                            Charge_old[i]=cnumber;
                                Charge_track[i]=cnumber2;
                                Charge_indices[i]=cnumber2;
                // 				Output2 = fopen(str2, "a");
                // 	fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
                // fclose(Output2);
                // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);

                        }
                }
                counter4++ ;
                cindex = Charge_track[k] ;
                //printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
                dx = rx[cindex] - rbix ;
                dy = ry[cindex] - rbiy ;
                dz = rz[cindex] - rbiz ;
                //printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
                dx -= L2xp*round(dy/L2yp);
                dy -= L2yp*round(dy/L2yp);
                dx -= L1xp*round(dx/L1xp);
                dz -= L*round(dz/L);

                dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
                dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
                dzt[k*printprops+tcount] += dz;
			}

			else
			{
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier2-0.5*(E_final - E_initial)))) // if greater, accept swap
				{
					count=0;
                    for(i = 0; i<NCharges; ++i)
                    {
                            if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
                            {count+=1;
                                Charge_old[i]=cnumber;
                                    Charge_track[i]=cnumber2;
                                    Charge_indices[i]=cnumber2;
                // 					Output2 = fopen(str2, "a");
                // 					fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
                // fclose(Output2);
                // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
                            }
                    }
                    counter5++ ;
                    cindex = Charge_track[k] ;
                    //printf("rxf: %lf rxi: %lf\n",rb[3*cindex],rbix);
                    dx = rx[cindex] - rbix ;
                    dy = ry[cindex] - rbiy ;
                    dz = rz[cindex] - rbiz ;
                    //printf("dx: %lf dy: %lf dz: %lf\n",dx,dy,dz);
                    dx -= L2xp*round(dy/L2yp);
                    dy -= L2yp*round(dy/L2yp);
                    dx -= L1xp*round(dx/L1xp);
                    dz -= L*round(dz/L);

                    dxt[k*printprops+tcount] += dx*cos(theta) - dy*sin(theta);
                    dyt[k*printprops+tcount] += dx*sin(theta) + dy*cos(theta);
                    dzt[k*printprops+tcount] += dz;
				}
				else
				{counter6++ ; //reject and swap back

					count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{count+=1;
							// printf("t= %d Charge_indices=%d\n",t,Charge_indices[i]);
						Charge_indices[i]=cnumber;
						}
					}

				}
			}
		}


	}//end of k loop
}//end of MC step


}//end of function
void calcDisFD()
{
	double dcomx,dcomy,dcomz,dxi,dyi,dzi;
	long double dx,dy,dz ;
	int i,j;
	double theta = atan(L1[1]/L1[0]);

	//CalcDis Chrg//
    double drelxi,drelyi;
	int chain ; 
    drelxi = 0; drelyi = 0;
	//if(t==0){printf("Warning - only hopping displacement considered\n");}
	tcount++ ;
	for(i=0;i<Ncharges;i++){
		j = Charge_track[i];
		dx = 0.0; dy = 0.0; dz = 0.0;

        //Calculate relative displacement between unwrapped chain + C.O.M
		chain = (int) floor(j/Nb) ; 
		//printf("j: %d chain: %d\n",j,chain);
        drelxi = comx[chain]-px[j] ; 
        drelyi = comy[chain]-py[j] ; 
		//printf("drelxi: %lf\n",drelxi);
        //Update flow force relative to unwrapped positions + C.O.M
        dx = dt*(fx[j] + flowrate*drelxi) + pc*R[3*j] ;
		dy = dt*(fy[j] - flowrate*drelyi) + pc*R[3*j+1];
		// dx = dt*(fx[i] + flowrate*(px[i]-comx[0])); //+ pc*R[3*j] ;
		// dy = dt*(fy[i] - flowrate*(py[i]-comy[0])); //+ pc*R[3*j+1];
		dz = dt*fz[j] + pc*R[3*j+2];
		////////////////////////////////////

		//Calculate Displacement Directly From Forces//
		// dx = dt*(fx[j] + flowrate*rx[j]); //+ pc*R[3*j] ;
		// dy = dt*(fy[j] - flowrate*ry[j]); //+ pc*R[3*j+1];
		// dz = dt*fz[j] ;  //+ pc*R[3*j+2];

		//Calculate Displacement Directly From Forces Relative to C.O.M//
		// dx = dt*(fx[j] + flowrate*(px[j]-comx[j])); //+ pc*R[3*j] ;
		// dy = dt*(fy[j] - flowrate*(py[j]-comy[j])); //+ pc*R[3*j+1];
		// dz = dt*fz[j] ;

		//printf("total_c: %.9Lf %.9Lf %.9Lf\n",dx,dy,dz);
		// printf("force_c: %.9lf %.9lf %.9lf\n",dt*fx[j],dt*fy[j],dt*fz[j]);
		// printf("flow_c: %.9lf %.9lf %.9lf\n",dt*(flowrate*rx[j]),dt*(flowrate*ry[j]),0.0);

		//printf("force_c: %.9lf %.9lf %.9lf\n",dt*fx[j],dt*fy[j],dt*fz[j]);
		//printf("flow_c: %.9lf %.9lf %.9lf\n",dt*flowrate*(rx[j]-comx[0]),dt*flowrate*(ry[j]-comy[0]),0.0);

		//printf("pos_c: %lu %.9lf %.9lf %.9lf\n",t,rx[j],ry[j],rz[j]);
		//printf("map_c: %lu %.9lf %.9lf %.9lf\n",t,sx[j],sy[j],sz[j]);

		dxt[i*printprops+tcount] = dx;
		dyt[i*printprops+tcount] = dy;
		dzt[i*printprops+tcount] = dz;

		//printf("drel: %lu %lf %lf\n",t,drelxi,drelyi);

		// dxt[i*printprops+tcount] = dx*cos(theta) + dy*sin(theta);
	  // dyt[i*printprops+tcount] = -dx*sin(theta) + dy*cos(theta);
	  // dzt[i*printprops+tcount] = dz ;

		//Take Displacement Relative to C.O.M

		// dxi = dx - dcomx;
	  // dyi = dy - dcomy;
	  // dzi = dz - dcomz;

		// dxt[i*printprops+tcount] = dxi*cos(theta) - dyi*sin(theta);
	  // dyt[i*printprops+tcount] = dxi*sin(theta) + dyi*cos(theta);
	  // dzt[i*printprops+tcount] = dzi ;
	}
}
void initPos(){
	int i,j;


	rbi[0] = comx[0];
	rbi[1] = comy[0];
	rbi[2] = comz[0];

}
void calcCH(){
	int i,j,ind2;
	long double dx,dy,dz;

	chfile = fopen(str2,"a");
	//printf("t: %lu\n",t);
	if(t==0){
		dx=0;
		dy=0;
		dz=0;
		for(j=0;j<Ncharges;j++){
			fprintf(chfile,"%lu %d %d %d %Lf %Lf %Lf\n",t,j,Charge_old[j],Charge_track[j],dx,dy,dz);
		}
		//fprintf(chfile,"%lu %lf %lf %lf\n",t,dx,dy,dz);
	}
	else{
		for(j=0;j<Ncharges;j++){
			dx=0;
			dy=0;
			dz=0;
			for(i=1;i<printprops+1;i++){
					ind2 = j*printprops+i ;
					//printf("j: %d t: %d dxt[%d]: %lf\n",j,i,ind2,dxt[ind2]);
					dx += dxt[ind2] ;
					dy += dyt[ind2] ;
					dz += dzt[ind2] ;
				}
			//printf("t: %lu j: %d dx[0]: %lf\n",t,j,dx);
			fprintf(chfile,"%lu %d %d %d %Lf %Lf %Lf\n",t,j,Charge_old[j],Charge_track[j],dx,dy,dz);
		}
		// for(i=1;i<printprops+1;i++){
		// 		//printf("j: %d t: %d dxt[%d]: %lf\n",j,i,ind2,dxt[ind2]);
		// 		dx += dxt[i] ;
		// 		dy += dyt[i] ;
		// 		dz += dzt[i] ;
		// 	}
		// //printf("t: %lu j: %d dx[0]: %lf\n",t,j,dx);
		// fprintf(chfile,"%lu %lf %lf %lf\n",t,dx,dy,dz);
	}
	fclose(chfile);
}

void calcMSD()
{
  int nsamp,ind1,ind2,count,m,l, wmax,Tmax,q,k,ctemp,n_it , n_nc,nc,ttemp1, i, j, jcount, nb1, n_nb1, nb2, n_nb2,nedot,n_nedot,n_nit,nit, dint,dint2;
  double dcomx,dcomy,dcomz,dx,dy,dz,rxtemp,rytemp,rztemp,MSDtemp,avgmsd;

  Tmax = tmax/printprops;
	wmax = Tmax+1;
	nsamp = Ncharges*wmax;


  double *msd = calloc(nsamp, sizeof(double));
  double *allmsd = calloc(nsamp, sizeof(double));
  //double *avgmsd = calloc(nsamp, sizeof(double));
  double *msdcount = calloc(nsamp, sizeof(double));
	double *vx = calloc(nsamp, sizeof(double));
	double *vy = calloc(nsamp, sizeof(double));
	double *vz = calloc(nsamp, sizeof(double));
	double *wx = calloc(nsamp, sizeof(double));
	double *wy = calloc(nsamp, sizeof(double));
	double *wz = calloc(nsamp, sizeof(double));
	double *Ree2 = calloc(Ncharges, sizeof(double));
	double *Ree2x = calloc(Ncharges, sizeof(double));
	double *Ree2y = calloc(Ncharges, sizeof(double));
	double *Ree2z = calloc(Ncharges, sizeof(double));
	double *Rcount = calloc(Ncharges, sizeof(double));
    double *avg = calloc(wmax, sizeof(double));
    double *avgx = calloc(wmax, sizeof(double));
    double *avgy = calloc(wmax, sizeof(double));
    double *avgz = calloc(wmax, sizeof(double));
	FILE *datafile;

    //Read in displacements from CH File///
	datafile = fopen(str2,"r");
	if(!datafile)
	{
			printf("Error - datafile not found, check data string %s\n",str2);
			exit(1);
	}

	for(i=0;i<nsamp;i++)
	{
			fscanf(datafile, "%d %d %d %d %lf %lf %lf\n", &ttemp, &ctemp,&dint2, &dint, &rxtemp, &rytemp, &rztemp);

			wx[i] = rxtemp;

			wy[i] = rytemp;

			wz[i] = rztemp;

	}
	fclose(datafile);

	// printf("nsamp: %d\n",nsamp);
	// for(i=0;i<nsamp;i++){
	// 	printf("wx[%d]=%lf\n",i,wx[i]); //correct value for displacemnt (agrees with rf - ri)
	// }


    //Sort displacements into correct order [t0,...,tmax,t0,...,tmax]
	for(j=0;j<Ncharges;j++){
		for(i=1;i<Tmax+1;i++){
			ind1 = j*Tmax+i;
			if(j>0){
				ind1 = j*(Tmax+1)+i;;
			}
			ind2 = i*Ncharges + j ;
			//printf("ind1: %d ind2: %d\n",ind1,ind2);
			vx[ind1] = wx[ind2] ;
			vy[ind1] = wy[ind2] ;
			vz[ind1] = wz[ind2] ;
		}
	}

	// printf("nsamp: %d\n",nsamp);
	// for(i=0;i<nsamp;i++){
	// 	printf("vx_new[%d]=%lf\n",i,vx[i]); //correct value for displacemnt (agrees with rf - ri)
	// }

    //Set Timestamp to Start MSD Calculation (in Tmax/tau)//
    
    //int MSDstart = 10000 ; 
    //printf("MSDstart: %d\n",MSDstart);
    //Calculate MSD//
    int nq ; 
    if(Ncharges>=25){
        nq=25 ;
    }
    else{
        nq = Ncharges;
    }
    #pragma omp parallel
    {
        double *Ree2 = calloc(nq, sizeof(double));
        double *Ree2x = calloc(nq, sizeof(double));
        double *Ree2y = calloc(nq, sizeof(double));
        double *Ree2z = calloc(nq, sizeof(double));
        double *Rcount = calloc(nq, sizeof(double));

        #pragma omp for

            for(i = 1; i<wmax-MSDstart; ++i) //i is shift
            {
                for(l=0;l<nq;l++){
                    Ree2[l] = 0 ;
                    Ree2x[l] = 0 ;
                    Ree2y[l] = 0 ;
                    Ree2z[l] = 0 ;
                    Rcount[l] = 0;
                }
                avg[i] = 0 ;
                avgx[i] = 0 ;
                avgy[i] = 0 ;
                avgz[i] = 0 ;
                for(q=0;q<nq;q++){
                    for(j = MSDstart; j<wmax-i; ++j)
                    {
                        dx = 0.0; dy = 0.0; dz = 0.0;
                        for(k = 0; k<i; ++k)
                        {
                            dx += vx[j+k+q*Tmax];
                            dy += vy[j+k+q*Tmax];
                            dz += vz[j+k+q*Tmax];
                        }
                        Ree2[q] += dx*dx+dy*dy+dz*dz;
                        Ree2x[q] += dx*dx;
                        Ree2y[q] += dy*dy;
                        Ree2z[q] += dz*dz;
                        Rcount[q]++;
                    }
                    avg[i]+=Ree2[q]/(double)Rcount[q];
                    avgx[i]+=Ree2x[q]/(double)Rcount[q];
                    avgy[i]+=Ree2y[q]/(double)Rcount[q];
                    avgz[i]+=Ree2z[q]/(double)Rcount[q];
                }
            }
    }

    //Write to File
    outputfile=fopen(outp1,"w");
    for(i=1;i<wmax-MSDstart;i++){
        fprintf(outputfile,"%f,%lf,%lf,%lf,%lf\n",i*dt*printprops,avg[i]/nq,avgx[i]/nq,avgy[i]/nq,avgz[i]/nq);
    }
    fclose(outputfile);


    //Calculate MSD - OG//

	// outputfile=fopen(outp1,"w");
	// for(i = 1; i<wmax-1; ++i) //i is shift
	// {
	// 	for(l=0;l<Ncharges;l++){
	// 		Ree2[l] = 0 ;
	// 		Ree2x[l] = 0 ;
	// 		Ree2y[l] = 0 ;
	// 		Ree2z[l] = 0 ;
	// 		Rcount[l] = 0;
	// 	}
	// 	avg = 0 ;
	// 	avgx = 0 ;
	// 	avgy = 0 ;
	// 	avgz = 0 ;
	// 	for(q=0;q<Ncharges;q++){
	// 		for(j = 1; j<wmax-i; ++j)
	// 		{
	// 			dx = 0.0; dy = 0.0; dz = 0.0;
	// 			for(k = 0; k<i; ++k)
	// 			{
	// 				dx += vx[j+k+q*Tmax];
	// 				dy += vy[j+k+q*Tmax];
	// 				dz += vz[j+k+q*Tmax];
	// 			}
	// 			Ree2[q] += dx*dx + dy*dy + dz*dz;
	// 			Ree2x[q] += dx*dx;
	// 			Ree2y[q] += dy*dy;
	// 			Ree2z[q] += dz*dz;
	// 			Rcount[q]++;
	// 		}
	// 		avg+=Ree2[q]/(double)Rcount[q];
	// 		avgx+=Ree2x[q]/(double)Rcount[q];
	// 		avgy+=Ree2y[q]/(double)Rcount[q];
	// 		avgz+=Ree2z[q]/(double)Rcount[q];
	// 	}
	// 	//average MSD over charges//
	// 	//print average to file//
	// 	fprintf(outputfile,"%f,%lf,%lf,%lf,%lf\n",i*dt*printprops,avg/Ncharges,avgx/Ncharges,avgy/Ncharges,avgz/Ncharges);

	// }
	// fclose(outputfile);

}
void calcCOM(){
	int r,c,j,ind;
	double dx, dy, dz,COMx,COMy,COMz,dxo,dyo,dzo,rxi,ryi,det,invDet,sxi,syi,szi;

	//Calculate C.O.M of Unwrapped Chain (in unit cube frame)//
	double comxi,comyi,comzi ; 
	for(i=0;i<Nc;i++){
		for(j=0;j<Nb;j++){
			COMx += px[Nb*i+j] ;
			COMy += py[Nb*i+j] ;
			COMz += pz[Nb*i+j] ;
		}
		COMx /= Nb ;
		COMy /= Nb ;
		COMz /= Nb ; 

		comxi = COMx;
		comyi = COMy;
		comzi = COMz ; 

        comx[i] = comxi;
		comy[i] = comyi;
		comz[i] = comzi ; 

	}

}
void calcMSDcom()
{
  int nsamp,ind1,ind2,count,m,l, wmax,Tmax,q,k,ctemp,n_it , n_nc,nc,ttemp1, i, j, jcount, nb1, n_nb1, nb2, n_nb2,nedot,n_nedot,n_nit,nit, dint,dint2;
  double dcomx,dcomy,dcomz,avg,avgx,avgy,avgz,dx,dy,dz,rxtemp,rytemp,rztemp,MSDtemp,avgmsd;

	Ncharges = 1;
  Tmax = tmax/printprops;
	wmax = Tmax+1;
	nsamp = Ncharges*wmax;


  double *msd = calloc(nsamp, sizeof(double));
  double *allmsd = calloc(nsamp, sizeof(double));
  //double *avgmsd = calloc(nsamp, sizeof(double));
  double *msdcount = calloc(nsamp, sizeof(double));
	double *vx = calloc(nsamp, sizeof(double));
	double *vy = calloc(nsamp, sizeof(double));
	double *vz = calloc(nsamp, sizeof(double));
	double *wx = calloc(nsamp, sizeof(double));
	double *wy = calloc(nsamp, sizeof(double));
	double *wz = calloc(nsamp, sizeof(double));
	double *Ree2 = calloc(Ncharges, sizeof(double));
	double *Ree2x = calloc(Ncharges, sizeof(double));
	double *Ree2y = calloc(Ncharges, sizeof(double));
	double *Ree2z = calloc(Ncharges, sizeof(double));
	double *Rcount = calloc(Ncharges, sizeof(double));
	FILE *datafile;

	datafile = fopen(str2,"r");
	if(!datafile)
	{
			printf("Error - datafile not found, check data string %s\n",str2);
			exit(1);
	}

	for(i=0;i<nsamp;i++)
	{
			fscanf(datafile, "%d %lf %lf %lf\n", &ttemp, &rxtemp, &rytemp, &rztemp);

			wx[i] = rxtemp;

			wy[i] = rytemp;

			wz[i] = rztemp;

	}
	fclose(datafile);

	// printf("nsamp: %d\n",nsamp);
	// for(i=0;i<nsamp;i++){
	// 	printf("wx[%d]=%lf\n",i,wx[i]); //correct value for displacemnt (agrees with rf - ri)
	// }

	for(j=0;j<Ncharges;j++){
		for(i=1;i<Tmax+1;i++){
			ind1 = j*Tmax+i;
			if(j>0){
				ind1 = j*(Tmax+1)+i;;
			}
			ind2 = i*Ncharges + j ;
			//printf("ind1: %d ind2: %d\n",ind1,ind2);
			vx[ind1] = wx[ind2] ;
			vy[ind1] = wy[ind2] ;
			vz[ind1] = wz[ind2] ;
		}
	}

	// printf("nsamp: %d\n",nsamp);
	// for(i=0;i<nsamp;i++){
	// 	printf("vx_new[%d]=%lf\n",i,vx[i]); //correct value for displacemnt (agrees with rf - ri)
	// }



  //Calculate MSD//

	outputfile=fopen(outp1,"w");
	for(i = 1; i<wmax-1; ++i) //i is shift
	{
		for(l=0;l<Ncharges;l++){
			Ree2[l] = 0 ;
			Ree2x[l] = 0 ;
			Ree2y[l] = 0 ;
			Ree2z[l] = 0 ;
			Rcount[l] = 0;
		}
		avg = 0 ;
		avgx = 0 ;
		avgy = 0 ;
		avgz = 0 ;
		for(q=0;q<Ncharges;q++){
			for(j = 1; j<wmax-i; ++j)
			{
				dx = 0.0; dy = 0.0; dz = 0.0;
				for(k = 0; k<i; ++k)
				{
					dx += vx[j+k+q*Tmax];
					dy += vy[j+k+q*Tmax];
					dz += vz[j+k+q*Tmax];
				}
				Ree2[q] += dx*dx + dy*dy + dz*dz;
				Ree2x[q] += dx*dx;
				Ree2y[q] += dy*dy;
				Ree2z[q] += dz*dz;
				Rcount[q]++;
			}
			avg+=Ree2[q]/(double)Rcount[q];
			avgx+=Ree2x[q]/(double)Rcount[q];
			avgy+=Ree2y[q]/(double)Rcount[q];
			avgz+=Ree2z[q]/(double)Rcount[q];
		}
		//average MSD over charges//
		//print average to file//
		fprintf(outputfile,"%f,%lf,%lf,%lf,%lf\n",i*dt*printprops,avg,avgx/Ncharges,avgy/Ncharges,avgz/Ncharges);

	}
	fclose(outputfile);

}
