/*
	Routines to solve for eq locations.
	eik dec 2016

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ray.h"

VelModel mp,ms ;
void printModel( char *text, double *x ) 
{
/* 	printf("Model %s : %9.4f %8.4f %8.4f %8.4f\n",text,x[0],x[1]*100.0,x[2]*100.0,x[3]) ; */
	printf("Model %s : %9.5f %9.5f %9.5f %9.5f\n",text,x[0],x[1],x[2],x[3]) ;
}
void insertRow( double *x, double *a, int i, int m, int n, double weight )
{
	while(n--) { a[i] = weight * *x++ ; i += m ; }
}
int  locate( Solution *sol, Phase *pp )
{
#define MAXPHASES 30 
	double x[4],x0[4],b[MAXPHASES],a[4*MAXPHASES] ;
	double dx[4], dxtest[4] ;
	double tt,dlon,distance,residual,sum ;
	double rayp,dtdx,dxdp,km2lat,km2lon ;
	double dz,ttz,dtdz,dtdla,dtdlo ;
	double dl,dtdlax,dtdlox,damp, weight ;
	VelModel *vm ;
	Phase *p ;
	Station *s ;
	int np,i,j,iter ;
	if (sol->index != pp->index ) {
		rLog(2,"Index in locate does not match",NULL) ;
		return 0 ;
	}
	p = pp ;
	while( p->index == sol->index ) p++ ;
	np = p-pp ;
	if( np > MAXPHASES ) rLog(1,"locate: more than %d phases", (void*) MAXPHASES ) ;
	printf("locate: %d phases\n",np) ;
	x0[0] = 0.0 ;         /* origin time, sec */
	x0[1] = sol->lat ;    /* latitude, degrees */
	x0[2] = sol->lon ;	/* longitude, degrees */
	x0[3] = sol->depth ;	/* depth, km */
	dz = 0.01 ;
	km2lat = 6391*M_PI/180.0 ;
	km2lon = km2lat * cos(sol->lat * M_PI/180.0 ) ;
	printModel("x0",x0) ;
	for ( iter = 0 ; iter < 12 ; iter++) {
	  sum = 0.0 ;
	  for( i = 0 ; i < np ; i++ ) {
		p = pp+i ;
		s = p->statP ;
		dlon = x0[2] - s->lon ;
		if( p->type == 'P' ) vm = &mp ; else vm = &ms ;
		if( p->type == 'P' ) weight = 1.0 ; else weight = 2.7 ;
		weight = 1.0 ;
		distance = gDistance(x0[1],s->lat,dlon) ;
		tt = timeFromDist(vm,distance,x0[3],&rayp,&dtdx,&dxdp) + x0[0] ;
		ttz = timeFromDist(vm,distance,x0[3]+dz,&rayp,&dtdx,&dxdp) +x0[0] ;
		dtdz = (ttz-tt)/dz ;
		dtdla = dtdx * ( x0[1] - s->lat )*km2lat*km2lat /  distance  ;
		dtdlo = dtdx * ( x0[2] - s->lon )*km2lon*km2lon /  distance  ;
		dx[0] = 1.0 ; dx[1] = dtdla ; dx[2] = dtdlo ; dx[3] = dtdz ;
		insertRow(dx,a,i,np,4,weight) ;
		b[i] = weight * ( p->pTime - tt ) ;

		dl = 0.0001 ;
		distance = gDistance(x0[1]+dl,s->lat,dlon) ;
		ttz = timeFromDist(vm,distance,x0[3],&rayp,&dtdx,&dxdp) +x0[0] ;
		dtdlax = (ttz-tt)/dl ;
		distance = gDistance(x0[1],s->lat,dlon+dl) ;
		ttz = timeFromDist(vm,distance,x0[3],&rayp,&dtdx,&dxdp) +x0[0] ;
		dtdlox = (ttz-tt)/dl ;

		residual = tt - p->pTime ;
		sum += residual*residual ;
		printf("%s %c %10.6f %10.6f %10.6f %10.6f %10.6f ",s->name,p->type,p->pTime,tt,residual,distance,x0[3]) ;
		printf("%8.4f  %8.4f %8.4f ",dtdlax,dtdlox,1.0/dtdx) ;
		printModel(" dx ",dx) ;
	   }
	   golubC(a,x,b,np,4) ;
	   printModel(" x ",x) ;
	   damp = 0.55 ;
	   if( iter > 7 ) damp = 1.0 ;
	   for( j = 0 ; j < 4 ; j++) x0[j] += damp*x[j] ;
	   printModel(" x0",x0) ;
	   printf("std_dev =%10.6f damp=%7.2f\n",sqrt(sum/np),damp) ;
	}
}

#ifdef TEST

char *pModel = "silp.vel" ;
char *sModel = "sils.vel" ;
char *phaseFile = "../geysir/phase.dat" ;
char *solFile = "../geysir/ctloc2" ;

doit()
{
	Phase *phases, *ip ;
	Solution *location, *lp ;
	int nPhases,nLoc ,i,j ;
	long long index,i2 ;
	Station *sp ;
	initVelModel(20,&mp ) ;
	initVelModel(20,&ms ) ;
	readVelModel(pModel,&mp) ;
	readVelModel(sModel,&ms) ;
	nPhases = readPhases(phaseFile,&phases ) ;
	nLoc = readCtloc(solFile,&location) ;
	lp = location + 5 ;
	index = lp->index ;
	printf("index=%ld\n",index) ;
	ip = phases ;
	while( ip->index < index ) ip++ ;
	j = locate( lp, ip ) ;
	while( ip->index == index ) {
		sp = ip->statP ;
		printf("%ld %s %10.6f %10.6f\n",ip->index,sp->name,ip->pTime,ip->weight ) ;
		ip++ ;
	}
	printf("%d phases\n",ip-phases ) ;
}
int main(int ac, char **av) {
	int cc,n ;
/*	feenableexcept(FE_INVALID) ; */
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"tpcr"))) {
	    switch(cc) {

	}}
	doit() ;
	return 0 ;
}
#endif
