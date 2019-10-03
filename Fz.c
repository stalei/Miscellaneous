//Author : Shahram Talei @ 2017 
// This code reads r and creates theta and fi for exponential disk. Saves
// the disk and then creates prob points above disk, calculates the force 
//for every point by direct sum and them saves the force inside another binary file.
// You need MPI to compile and run the code.
/*
 * 
 * */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_sf_bessel.h>
#include<string.h>

#define G 6.67e-8 // gravity constant in cgs units

// Fg calculates direct force between two points
double Fg(double r0, double fi0, double z0,double r1, double fi1, double z1 )
{
	double f,lsq,l;
	// convert from kpc to cm
	//r0*=3.086e21;
	//r1*=3.086e21;
	//z0*=3.086e21;
	//z1*=3.086e21;
	lsq=r1*r1+r0*r0-2*r1*r0*cos(fi0+fi1)+(z1-z0)*(z1-z0);
	l=pow(lsq+0.0005,1.5); // adding small softening to make smooth function
	f=(z1-z0)/l; 
	// ignoring G because we are calculating F/G[Sigma0],
	//we need to multiply by cos (theta) to project the force in z direction f= G*(z1-z0)/l
	return f;
}

//this is used to define disk points
struct Mass
{
	double R;
	double Fi;
	double Z;
};
//this is used to define probe points
struct Point
{
	double r;
	double fi;
	double z;
	double fz;
};

///////////////////main//////////
int main()
{
	MPI_Init(NULL,NULL);
	/////////////////////////////// Numbers to modify///////////////////
	///////////////////////////////////////////////////////////////////
	//make sure these numbers are the same as MakeR.py
	const int N=200000; //number of particles in disk
	const int M=14; // number of probe points in length
	const int C=50; // core index, just to check how code distributes particles on every core
	double R0=2.8,Z0=0.48; // disk parameters
	float radprobe,zprobe; // probe size
	radprobe=10;
	zprobe=40;
	
	char diskaddress []="Disk.out"; // this is the file we save position of particles in disk
	char faddress []="PF.bin"; // this is the file we save force in every probe point
	char RFile[]="R.bin";// this is input file we read interpolated Rs from python code
	
	////////////////////////////////////////////////
	const double PI=3.14159265359;
	////////////////////////////////////////////////
	struct Mass m[N];
	struct Point p[M][M][M]; //Probe points
	double R[N],Fi;
	long int NCore[C];
    MPI_File DiskFile, DevFile1, DevFile2, ForceFile;
    MPI_Status status;
	//measure time
	double t1, t2; 
	t1 = MPI_Wtime(); 
   srand(time(NULL));
   int i=0,j=0,k=0,pp=0;
	double FZ,Fz;
	double zmax,rmax;
	double zrand,zcheck;
	double fii,ri,zi;
	
	const int binnumber=20;
	float binR[binnumber], binRSum[binnumber];
	
	int CoreNum;
	MPI_Comm_rank(MPI_COMM_WORLD, & CoreNum);
	int CoreTot;
	MPI_Comm_size(MPI_COMM_WORLD, & CoreTot);
	rmax=R0*radprobe;
	zmax=Z0*zprobe;
	
	for(i=0;i<C;i++)
		NCore[i]=0;
	//print initial interface
	if(CoreTot>1)
      {
		  if(CoreNum==0)
		  {
		  printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@/Fz\\@@@@@@@@@@@@@@@@@@@@@@@@@@@");
		  printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
		  printf("\n\tSimple code to Calculate vertical force above exponential Disk");
		  printf("\n This code, creates an exponential disk and finds forces ");
		  printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
		  printf("\nsimulation parameters are:");
		  printf("\n\tNumber of Disk Particles:%d\n\tNumber of force probe points:%d * %d * %d= %d points", N,M,M,M, M*M*M);
		  printf("\n\tDisk parameters{Rd:%g kpc,Zd:%g kpc}",R0,Z0);
		  printf("\n\tProbe radius is:%g X Rd around disk",radprobe);
		  printf("\n\t---------------------------------------------------\n Creating Disk ...");
		  }
		  ///////////////////////////Create disk ////////////////////////////////////
		  //read R's
		  FILE *Rfile;
		  Rfile=fopen(RFile,"rb");
		  MPI_File_open(MPI_COMM_WORLD,diskaddress,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&DiskFile);
		 for(i=0;i<binnumber;i++)
		 {
			 binR[i]=0;
		 }
		 
		 for(pp=CoreNum;pp<N;pp+=CoreTot)
		  {
			  fread(&R[pp],sizeof(double),1,Rfile);
			  Fi=((double)rand()/RAND_MAX)*2*PI;
			  zrand=(double)rand()/RAND_MAX;
			  zcheck=zrand*10000;
			  if(((int)zcheck)%2==0)
				m[pp].Z=-Z0*log(zrand);
			  else
				  m[pp].Z=Z0*log(zrand);
			  m[pp].R=R[pp];
			  m[pp].Fi=Fi;
			  NCore[CoreNum]++;
			  binR[((int)((R[pp]/(radprobe*R0))*binnumber))]++; // print hist at the end, after force to keep synchronization
		  }	
		fclose(Rfile);
		  //Save Disk output file
		if(CoreNum==0)
			  printf("\n Disk is Dine, Please wait...");
		  //{
			  ///////////Check distribution
			for(i=0;i<binnumber;i++)
			{
				MPI_Reduce(&binR[i],&binRSum[i],1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
			}
			  for(pp=CoreNum;pp<N;pp+=CoreTot)
			  {
				MPI_File_write(DiskFile,&m[pp],3,MPI_DOUBLE,MPI_STATUS_IGNORE);  
			  }			  
		  //}		  
		  MPI_File_close( &DiskFile );
		  if(CoreNum==0)
			  printf("\nDisk is saved !\n Calculating Force ...");
		  
		  ///////////////////////////////////////////////////////////Calculate point to point force //////
		  //find forces on a grid
		  for(i=0;i<M;i++)
		  {
			  for(j=0;j<M;j++)
			  {
				for(k=0;k<M;k++)
				{
					p[i][j][k].r=j*(rmax/M)+0.5;
					p[i][j][k].fi=i*(2*PI/M);
					p[i][j][k].z=k*(zmax/M)+0.1;
					for(pp=CoreNum;pp<N;pp+=CoreTot)
					{
						Fz+=2.0*PI*R0*R0*Fg(m[pp].R,m[pp].Fi,m[pp].Z,p[i][j][k].r,p[i][j][k].fi,p[i][j][k].z)/N;
					}
					MPI_Reduce(&Fz,&FZ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					if(FZ>1.0e-100)
					  p[i][j][k].fz=FZ;
					else
						p[i][j][k].fz=0;
					printf("Force is:%g\n",p[i][j][k].fz);
					Fz=0;
				}
			  }
		  }
		  if(CoreNum==0)
			  printf("\nForce is Done !");
		 
	   if(CoreNum==0)
	   {
			FILE *PF_file;
			PF_file=fopen(faddress,"wb");
			if (!PF_file)
				{
					printf("Unable to open file!");
				}
			else
			{
			for(i=0;i<M;i++)
				for(j=0;j<M;j++)
					for(k=0;k<M;k++)
						fwrite(&p[i][j][k], sizeof(struct Point), 1, PF_file);
					
			printf("\nForce file is saved !");		
			}
			fclose(PF_file);
			//simple output to check distributions and some numbers
			printf("\n Numbers of particles/bin:\n\t ___________total number - OutR _______________________\n");
			double scale;
			double surface,surfden,sigma0=0;
			scale=radprobe*R0/binnumber;
			for(i=0;i<binnumber;i++)
			{
				printf("\t|");
				for(j=0;j<(int)(binRSum[i]*180/N);j++)
					printf("#");
				printf(" %g - %g kpc\n",binRSum[i],(i+1)*scale);
			}
			
			printf("\n Check Constant Surface density Coefficient:\n\t ______________Surface - Sigma0-ish - OutR _________________________________\n");
			for(i=0;i<binnumber;i++)
			{
				surface=PI*((2*i+1)*scale*scale);
				surfden=binRSum[i]*exp((i+1)*scale/R0)/surface;
				sigma0+=surfden;
				printf("\t|");
				for(j=0;j<(int)(surfden*100/N);j++)
					printf("#");
				printf(" %g - %g - %g Rd \n",surface,surfden,(i+1)*scale/R0);
			}
			printf("sum of surfacedensity=%g",sigma0);
			//print bin diagram for 
			printf("\n Surface Number density (Number of particles per bin area):\n\t ______________ N[bin]/Area - OutR _____________________________\n");
			for(i=0;i<binnumber;i++)
			{
				surface=PI*((2*i+1)*scale*scale);
				//printf("\nsurf:%g",)
				surfden=binRSum[i]/surface;
				printf("\t|");
				for(j=0;j<(int)(surfden*120/N);j++)
					printf("#");
				printf(" %g - %g Rd \n",surfden,(i+1)*scale/R0);
			}
		}
	  /////////////////////
		t2 = MPI_Wtime();
		
		printf("\nNumber of particles in core %d is:%ld\n",CoreNum,NCore[CoreNum]);
		if(CoreNum==0)
		{
			printf("\nDone !");
			printf( "\nElapsed time is %f sec\n", t2 - t1 );
		}
		  			
         MPI_Finalize();
	}
	else
	{
		printf("\nNot enough Cores");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
return 0;
   
}
