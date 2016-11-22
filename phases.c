
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ray.h"
VelModel mS, mP ;

void readPhases( char *filename)
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
	shLogLevel = 2 ;
	readPhases("../geysir/phase.dat") ;
	return 0 ;
}
