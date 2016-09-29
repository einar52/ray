/* 
	Program to print out traveltime data for 1d linear gradeint layered models.
	Einar Kjartransson, 2016 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ray.h"

char *velFile = "silp.vel" ;
double sourceDepth = 4.5 ;
double bottomDepth = 15 ;
double velReduce = 7.0 ;
int nPoint = 10 ;

char outputName[300] ;
VelModel m ;

void doIt() 
{
	double pMax, pBottom, p , dp, x,t, reduce ;
	int il,i, mmode ;
	FILE *ofd ;
	ofd = fopen(outputName,"w") ;
	if( velReduce > 0.0 ) reduce = 1.0/velReduce ; else reduce = 0.0 ;
	initVelModel(200, &m  ) ;
	readVelModel(velFile, &m) ;
	if( shLogLevel > 3 ) printVelModel( &m) ;
	pMax = 1.0/velZ(sourceDepth,&m,&il ) ;
	mmode = RayDown ;
	if( sourceDepth <= 0.0 ) mmode = SURFACE ; 
	else for( i = 0 ; i < nPoint ; i++) {
		p = (i + 0.9) *pMax / nPoint ;
		x = traceUD( RayUP, p, sourceDepth, &m, &t ) ;
		fprintf(ofd,"%10.5f %10.5f %6d\n",x,t-x*reduce,i) ;
	}
	pBottom = 1.0/velZ(bottomDepth,&m,&il) ;
	dp = pMax - pBottom ;
	for( i = 0 ; i < nPoint ; i++) {
/*		p = pMax - (i+0.1)*dp/nPoint ;  */
		p = 1.0/velZ( sourceDepth + (i+0.5)*(bottomDepth-sourceDepth)/nPoint,&m,&il);   
		x = traceUD( mmode , p, sourceDepth, &m, &t ) ; 
		fprintf(ofd,"%10.5f %10.5f %6d\n",x,t-x*reduce,i) ;
	} 
	fclose(ofd) ;
}
void makeOutputName( int i , char *a ) 
{
	char  *b ;
	b = outputName ; 
	if(*a == '.') a++ ;
	if(*a == '/') a++ ;
	b=strcpy(b,"PL_") ; b += 3 ;
	while ( i > 0 ) {
		*b = *a++ ;
		if( *b == 0 ) { *b = '_' ; i-- ; }
		if( *b == '-' ) b-- ;
		b++ ;
	} 
	*(--b) = 0 ;
}
int main(int ac , char **av ) 
{
	extern char *optarg ;
	extern int optind ;
	int cc ;
	while( EOF != (cc = getopt(ac,av,"f:l:d:n:b:v:o:hH?"))) {
		switch(cc) {
		case 'f':	velFile=optarg ; break ;
		case 'l' :	shLogLevel = atoi(optarg) ; break ;
		case 'd' :      sourceDepth = atof(optarg) ; break ;
		case 'n' : 	nPoint = atoi(optarg) ; break ;
		case 'b' : 	bottomDepth = atof(optarg) ; break ;
		case 'v' :	velReduce = atof(optarg) ; break ;
		case 'o' :	strcpy( outputName,optarg) ; break ;
	}}
	if( 0 == *outputName )  makeOutputName(ac,*av) ; 
	fprintf(stderr,"_%s_\n",outputName) ;
	doIt() ;
	return 0 ; 
}
