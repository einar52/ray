
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
int readPhase(char *fileName, Phase *phases )
{
#define PHASELENGTH 31
	int size,nPhase,i,nf ;
	long long sIndex ;
	char *line[300] ;
	FILE *ffd ;
	Phase *pp ;
	size = testFileSize(fileName) ;
	nPhase = size/PHASELENGTH ;
	phases = (Phase*) malloc(sizeof(Phase)*nPhase) ;
	pp = phases ;
	ffd = fopen(fileName,"r") ;
	while( nf = readData(ffd,line)) {
		if( nf > 4 ) {
			sIndex = strtoll(line[14],NULL,10) ;
		} else {
			pp->index = sIndex ;
			pp->statP = lookUpStation(line[0]) ;
			pp->type = *line[1] ;
			pp->pTime = atof(line[1]) ;
			pp->weight = atof(line[2]) ;
		}
		pp++ ;
	}
	return phases-pp ;
}
int readCtloc(char *fileName, Solution *solutions )
{
#define CTLINELENGTH 133 
	int size ,nSol, i,ii, nf ;
	char *line[300] ;
	FILE *ffd ;
	Solution *sp ;
	size = testFileSize(fileName) ;
	nSol = size/CTLINELENGTH ;
	ffd = fopen(fileName,"r") ;
	solutions = (Solution*)  malloc(sizeof(Solution) * (nSol+5) ) ;
	sp = solutions ;
	while( nf = readData(ffd,line) ) {
		sp->index = strtoll(line[0],NULL,10) ;
		sp->lat = atof(line[3]) ;
		sp->lon = atof(line[4]) ;
		sp->depth = atof(line[5]) ;
		sp++ ; 
	}
	return sp-solutions ;
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
/*	feenableexcept(FE_INVALID) ; */
	shLogLevel = 2 ;
	testPhases("../geysir/phase.dat") ; 
	readCtloc("../geysir/ctloc2",sol) ;
	return 0 ;
}
#endif
