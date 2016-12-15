/*
	Routines to solve for eq locations.
	eik dec 2016

*/
char *helpText = "\n\
reltest options\n\
	Trace rays to check location list. One of -r and -c must be given\n\
	options include:\n\
	-r  file 			Name of file listing relative locations	\n\
	-c  file 			Name of file listing catalog locations	\n\
	-p  file	phase.dat	Name of file listing phases\n\
	-P  file	silp.vel	P velocity model \n\
	-S  file	sils.vel	S velocity model \n\
	-v				List more information, may be repeated\n\
	-h				print instructions\n\
\n\
\n" ;
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ray.h"

VelModel mp,ms ;


char *pModel = "silp.vel" ;
char *sModel = "sils.vel" ;
char *phaseFile = "phase.dat" ;

/*char *solFile = "../geysir/ctloc2" ; 
char *solFile = "../geysir/reloc2" ; */
Phase *phases ;
Solution *location ;
int nPhases,nLoc ;
int phaseList ;

void checkPhases() 
{
	int iphase,il,np ;
	Phase *pp, *p1, *p ;
	Station *s ;
	Solution *sol ;
	VelModel *vm ;
	double rayp,dtdx,dxdp ;
	double dlon,distance,tt,dtt ;
	int numberP,numberS,numberPe,numberSe ;
	double sumP,sumS, sumPe,sumSe ;
	sol = location ;
	pp = phases ;
	printf("%d phases, %d solutions\n",nPhases,nLoc) ;
	numberP = 0 ; numberS = 0 ; sumS = 0.0 ; sumP = 0.0 ;
	for( il = 0 ; il < nLoc ; il++) {
		sol = location + il  ;
		while( pp->index < sol->index) pp++ ;
		p1 = pp ;
		while( pp->index == sol->index) pp++ ;
		np = pp - p1 ;
		if( phaseList > 1 )printf( " %3d %3d %3d \n",il,pp-p1,pp-phases) ;
		numberPe = numberP ; numberSe = numberS; sumPe = sumP ; sumSe = sumS ;
		for ( p = p1 ; p < pp ; p++) {
			s = p->statP ;
			if( p->type == 'P' ) vm = &mp ; else vm = &ms ;
			dlon = sol->lon - s->lon ;
			distance = gDistance(sol->lat,s->lat,dlon) ;
			tt = timeFromDist(vm,distance,sol->depth,&rayp,&dtdx,&dxdp) ;
			dtt = p->pTime-tt ;
			if( phaseList ) printf("%ld %s %c %8.4f %8.2f %8.4f %8.4f\n",
				p->index,s->name,p->type,p->pTime,distance,tt,dtt ) ;
			if( p->type == 'P') 
				{ numberP++, sumP += dtt * dtt ; }
			else 
				{ numberS++, sumS += dtt * dtt ; }
		}
		sumPe = sumP - sumPe ; 
		sumSe = sumS - sumSe ;
		numberPe = numberP - numberPe ; 
		numberSe = numberS - numberSe ;
		if( phaseList ) {
			if( numberPe ) printf("%4d P phases, RMS residual: %6.3f\n",numberPe, sqrt(sumPe/numberPe)) ;
			if( numberSe ) printf("%4d S phases, RMS residual: %6.3f\n",numberSe, sqrt(sumSe/numberSe)) ;
		}
	}
	printf("%4d P phases, RMS residual: %6.3f\n",numberP, sqrt(sumP/numberP)) ;
	printf("%4d S phases, RMS residual: %6.3f\n",numberS, sqrt(sumS/numberS)) ;

/*
	
	for (iphase = 0 ; iphase < nPhases-100 ; iphase++ ) {
		printf("%4d %4d %ld %ld \n", iphase,sp-location,sp->index, pp->index ) ;
		while( sp->index < pp->index ) sp++ ;
		pp++ ;
	}
*/
}

void rData()
{
	int i,j ;
	long long index,i2 ;
	Station *sp ;
	initVelModel(20,&mp ) ;
	initVelModel(20,&ms ) ;
	readVelModel(pModel,&mp) ;
	readVelModel(sModel,&ms) ;
	nPhases = readPhases(phaseFile,&phases ) ;
/*	nLoc = readCtloc(solFile,&location) ; 
	nLoc = readReloc(solFile,&location) ; */
/*
	lp = location + 6 ;
	index = lp->index ;
	printf("index=%ld\n",index) ;
	ip = phases ;
	while( ip->index < index ) ip++ ;
	while( ip->index == index ) {
		sp = ip->statP ;
		printf("%ld %s %10.6f %10.6f\n",ip->index,sp->name,ip->pTime,ip->weight ) ;
		ip++ ;
	}
	printf("%d phases\n",ip-phases ) ;
*/
}
int main(int ac, char **av) {
	extern char *optarg ;
	extern int optind ;
	int cc,n ;
/*	feenableexcept(FE_INVALID) ; */
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"vS:P:p:c:r:hH?"))) {
	    switch(cc) {
		case 'c' : nLoc = readCtloc(optarg,&location) ; break ;
		case 'r' : nLoc = readReloc(optarg,&location) ; break ;
		case 'p' : phaseFile = optarg ; break ;
		case 'P' : pModel = optarg ; break ;
		case 'S' : sModel = optarg ; break ;
		case 'v' : phaseList++ ; break ;
		default : printf("%s\n",helpText ) ;
	}}
	rData() ;
	checkPhases() ;
	return 0 ;
}
