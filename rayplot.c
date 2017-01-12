/*
	Program to plot rays
	eik dec 2016

*/
char *helpText = "\n\
rayplot options\n\
	Write files that can be used to plot raypaths.\n\
	Options and [defaults] are:\n\
	-P 	velfile	[plot.vel]	Name of velocity file\n\
	-d	dz			Resample velocity file to dz\n\
	-s	dz			Resample using cubic splines do dz\n\
	-n	nz	[100]		Number of resampled velocoity values\n\
	-b	zBottom [7.45]		Depth were ray turns\n\
	-r	rayFile	[ray.points]	\n\
	-h				print instructions\n\
\n" ;
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define _GNU_SOURCE
#include <fenv.h>
#include <unistd.h>
#include <signal.h>
#include "ray.h"

extern DepthPoint *rayPoints ;
extern int nRaypoints ;
VelModel mp ;
double zSource ;
double dz  ;
double zBottom = 7.45 ;
char *rayFile = "ray.points" ;
int nz = 100 ;
char *pModel = "plot.vel" ;
int splineFlag ;

void doIt()
{
	FILE *ff ;
	double z,x,t,p,xMax ;
	int il,i ;
	VelModel m ;
	DepthPoint *rp ;
	if( dz > 0.0 ) {
		if (splineFlag ) m = resampleVelModel(&mp,dz,nz) ;
		else m = resampleVelModel2(&mp,dz,nz) ;
		m.nVel = nz ;
		printVelModel(&m) ;
	}
	else    m = mp ;
	ff = fopen(rayFile,"w") ;
	rp = rayPoints ;
	p = 1.0/velZ(zBottom,&m,&il) ;
	x = traceUD(SURFACE,p,zSource,&m,&t) ;
	for ( i =  0 ; i < nRaypoints ; i++ ) {
		fprintf(ff,"%10.3f %10.3f\n",rp->x,rp->z)  ;
		rp++ ;
	}
	rp-- ;
	xMax = 2.0 * rp->x ;
	for( i = 1 ; i < nRaypoints ; i++) {
		rp-- ;
		fprintf(ff,"%10.3f %10.3f\n",xMax - rp->x,rp->z)  ;
	}
}
void init()
{
	int i,j,size ;
	long long index,i2 ;
	Station *sp ;
	initVelModel(20,&mp ) ;
	readVelModel(pModel,&mp) ;
	size = 100 + 2 * mp.nVel ;
	rayPoints = ( DepthPoint *) malloc(size* sizeof(DepthPoint)) ;
	
}
int main(int ac, char **av) {
	extern char *optarg ;
	extern int optind ;
	int cc,n ;
	feenableexcept(FE_INVALID) ;  
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"d:n:b:P:s:r:hH?"))) {
	    switch(cc) {
		case 'P' : pModel = optarg ; break ;
		case 'd' : dz = atof(optarg) ; break ;
		case 's' : dz = atof(optarg) ; splineFlag++ ; break ;
		case 'n' : nz = atoi(optarg) ; break ;
		case 'b' : zBottom = atof(optarg) ; break ;
		case 'r' : rayFile = optarg ; break ;
		default : printf("%s\n",helpText ) ;
	}}
	init() ;
	doIt() ;
	return 0 ;
}
