
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
int readCtloc(char *fileName, Solution *solutions )
{
#define CTLINELENGTH 133 
	int size ,nSol, i,ii, nf ;
	char *line[300] ;
	FILE *ffd ;
	Solution *sp ;
/*	char *fields[200] ; */
	size = testFileSize(fileName) ;
	nSol = size/CTLINELENGTH ;
	solutions = (Solution*)  malloc(sizeof(Solution) * (nSol+5) ) ;
	sp = solutions ;
	ffd = fopen(fileName,"r") ;
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
	while(nf = readData(fd,line ) ) {
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

int main(int ac, char **av) {
	Solution  *sol ;
/*	feenableexcept(FE_INVALID) ; */
	shLogLevel = 2 ;
	readCtloc("../geysir/ctloc2",sol) ;
	testPhases("../geysir/phase.dat") ; 
	return 0 ;
}
