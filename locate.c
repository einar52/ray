/*
	Routines to solve for eq locations.
	eik dec 2016

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ray.h"

int  locate( Solution *sol, Phase *pp )
{
#define MAXPHASES 30 
	double x[4],x0[4],b[MAXPHASES],a[4*MAXPHASES] ;
	Phase *p ;
	int np,i ;
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
	for( i = 0 ; i < np ) {
	}
}

#ifdef TEST

char *sModel = "sils.vel" ;
char *pModel = "silp.vel" ;
char *phaseFile = "../geysir/phase.dat" ;
char *solFile = "../geysir/ctloc2" ;

doit()
{
	VelModel mp, mP ;
	Phase *phases, *ip ;
	Solution *location, *lp ;
	int nPhases,nLoc ,i,j ;
	long long index,i2 ;
	Station *sp ;
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
