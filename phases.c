
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <fenv.h>
#include "ray.h"

VelModel mS, mP ;
off_t testFileSize(char *fileName)
{
	int fd;
	off_t size ;
	fd = open(fileName,O_RDONLY) ;
	size = lseek(fd,0,SEEK_END) ;
	(void) close(fd) ;
	return size ;
}
int readPhases(char *fileName, Phase **phases )
{
#define PHASELENGTH 31
	int size,nPhase,i,nf ;
	long long sIndex ;
	char *line[300] ;
	FILE *ffd ;
	Phase *pp ;
	size = testFileSize(fileName) ;
	nPhase = size/PHASELENGTH ;
	*phases = (Phase*) malloc(sizeof(Phase)*nPhase) ;
	pp = *phases ;
	ffd = fopen(fileName,"r") ;
	if( ffd == NULL) rLog(1,"cannot open %s", (void *)fileName) ;
	while( nf = readData(ffd,line)) {
		if( nf > 4 ) {
			sIndex = strtoll(line[14],NULL,10) ;
		} else {
			pp->index = sIndex ;
			pp->statP = lookUpStation(line[0]) ;
			pp->type = *line[3] ;
			pp->pTime = atof(line[1]) ;
			pp->weight = atof(line[2]) ;
			pp++ ;
		}
	}
	return pp-*phases ;
}
void printPhases( Phase *p, int n  ) 
{
	int i ;
	printf("n = %d\n",n) ;
	for ( i = 0 ; i< n ; i++) {
/*		printf("%s %c %ld\n",p->statP->name,p->type,p->index) ; */
		printf("%d %ld %s %8.3f %6.3f %c\n",i,p->index,p->statP->name,p->pTime,p->weight,p->type) ;
		p++ ;
	}
}
void printSol( Solution *p, int n )
{
	int i ;
	for(i = 0 ; i < n ; i++) {
		printf("%ld %10.6f %10.6f %8.3f %10.4f\n",
			p->index,p->lat,p->lon,p->depth,p->timeShift) ;
		p++ ;
	}
}
int readReloc(char *fileName, Solution **solutions)
#define RELLINELENGTH 183
{
	int size ,nSol, i,ii, nf ;
	char *line[300] ;
	FILE *ffd ;
	Solution *sp ;
	double ts,ti,dt ;
	size = testFileSize(fileName) ;
	nSol = size/RELLINELENGTH ;
	ffd = fopen(fileName,"r") ;
	if( ffd == NULL) rLog(1,"cannot open %s", (void *)fileName) ;
	*solutions = (Solution*)  malloc(sizeof(Solution) * (nSol+5) ) ;
	sp = *solutions ;
	while( nf = readData(ffd,line)) {
		sp->index = strtoll(line[0],NULL,10) ;
		sp->lat = atof(line[1]) ;
		sp->lon = atof(line[2]) ;
		sp->depth = atof(line[3] ) ;
		ts = atof(line[15] ) ;
		ti = (sp->index % 100000 ) * 0.001 ;
		dt = ts - ti ;
		if( dt >  30.0 ) dt -= 60.0 ;
		if( dt < -30.0 ) dt += 60.0 ;
		sp->timeShift = dt ;
		sp++ ;
/*		printf("%s %s %8.4f %8.4f %8.4f\n",line[0],line[1], ts,ti,dt) ; */
	}
	return sp-*solutions ;
}
int readCtloc(char *fileName, Solution **solutions )
{
#define CTLINELENGTH 133 
	int size ,nSol, i,ii, nf ;
	char *line[300] ;
	FILE *ffd ;
	Solution *sp ;
	size = testFileSize(fileName) ;
	nSol = size/CTLINELENGTH ;
	ffd = fopen(fileName,"r") ;
	if( ffd == NULL) rLog(1,"cannot open %s", (void *)fileName) ;
	*solutions = (Solution*)  malloc(sizeof(Solution) * (nSol+5) ) ;
	sp = *solutions ;
	while( nf = readData(ffd,line) ) {
		sp->index = strtoll(line[0],NULL,10) ;
		sp->lat = atof(line[3]) ;
		sp->lon = atof(line[4]) ;
		sp->depth = atof(line[5]) ;
		sp++ ; 
	}
	return sp-*solutions ;
/*
	ii = 0 ;
	for( i = 0 ; i < nSol; i++) {
		sp = solutions + i ;
		printf("%3d %ld %19ld %9.5f %9.5f %8.3f\n",
			i,sp->index,sp->index-solutions[ii].index,sp->lat,sp->lon,sp->depth) ;
		ii = i ;
	}
*/
}
void testPhases( char *filename)
{
	FILE *fd ;
	int nf ;
	double lat,lon,depth,dist,tim ;
	double timModel, p,dtdx,dxdp,z,dtt ;
	char *line[300], phase ;
	VelModel *mp ;
	Station *sta ;
	initVelModel(200,&mS) ;
	initVelModel(200,&mP) ;
	readVelModel("silp.vel",&mP) ;
	readVelModel("sils.vel",&mS) ;
	fd = fopen(filename,"r") ;
	if( NULL == fd ) rLog( 0,"unable to open phase file %s",(void *) filename) ;	
	while(2<(nf = readData(fd,line ) )) {
		if( nf > 4 ) {
			lat = atof(line[7]) ;
			lon = atof(line[8]) ;
			depth = atof(line[9]) ;
	 		printf("EV %d %8.4f %8.4f %8.2f %s\n",nf,lat,lon,depth,line[14]) ;
		} else {
			sta = lookUpStation(line[0]) ;
			phase = *line[3] ;
			if( phase == 'P') mp = &mP ;
			else  mp = &mS ;
			tim = atof(line[1]) ;
			dist = gDistance(lat,sta->lat,lon-sta->lon) ;
			timModel = timeFromDist(mp,dist,depth, &p,&dtdx,&dxdp) ;
			dtt = tim-timModel ;
			printf("%s %c %8.4f %8.4f %8.3f %8.3f %6.3f %8.3f\n",
				sta->name,phase,sta->lat,sta->lon,dist,tim,timModel,dtt) ;
		}
	}
}
#ifdef TEST
int main(int ac, char **av) {
	Solution  *sol ;
	Phase *phases ;
	int cc,n ;
/*	feenableexcept(FE_INVALID) ; */
	shLogLevel = 2 ;
	while( EOF != ( cc = getopt(ac,av,"tpcr"))) {
	    switch(cc) {
		case 't' : testPhases("../geysir/phase.dat") ;  break ;
		case 'c' : readCtloc("../geysir/ctloc2",&sol) ; break ;
		case 'r' : n = readReloc("../geysir/reloc2",&sol) ; 
				printSol(sol,n) ;
				break ;
		case 'p' :  n = readPhases("../geysir/phase.dat",&phases) ;
				printPhases(phases,n) ;
				printf("n=%d\n",n) ;
			break ;
	}}
	return 0 ;
}
#endif
