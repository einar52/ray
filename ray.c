#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int shLogLevel = 4 ;
typedef struct {
	int nVel ;
	double *z ; /* depth  */
	double *v ; /* velocity */
} VelModel ;
#define MaxLayer 1000
void  rLog( int level, char *s1 , void *p )
{
        char buff[200] ;
        time_t t ;
        FILE *lf ;
        if( level > shLogLevel ) return ;
        time(&t) ;
        (void) strftime(buff,30,"%d/%m %H:%M:%S",gmtime(&t)) ;
/*        lf = fopen(logFile,"a") ;  */
	lf = stderr ;
        fprintf(lf,"rayt, %s l=%d: ",buff,level) ;
        fprintf(lf,s1,p) ;
        fprintf(lf,"\n") ;
        fclose(lf) ;
        if( level < 0 ) abort() ;
        if( level < 1 ) exit(-1) ;
}
void initVelModel( int nVel , VelModel *m )
{
	m->nVel = nVel ;
	m->z = ( double *) calloc( sizeof(double ) , nVel ) ;
	m->v = ( double *) calloc( sizeof(double ) , nVel ) ;
}
void printVelModel( VelModel *m )
{
	int i ;
	printf("nVel = %d\n",m->nVel ) ;
	for( i = 0 ; i < m->nVel ; i++) 
		printf("%10.4f %10.4f %6d \n",*(m->z+i),*(m->v+i),i ) ;
}
void readVelModel( char *inputFile, VelModel *m )
/*
	Velocity model is defined by velocity at specified depth. First depth
	value is taken to be zero. Velocity below max depth givine is linearly
	extrapolated using gradient between last two given depths.
*/
{
	FILE *ifil ;
	char line[300] ;
	double v[MaxLayer], z[MaxLayer], vold, zold ;
	int i,j, more ;
	ifil = fopen(inputFile,"r") ;
	if( NULL == ifil ) rLog(0,"Cannot open %s",( void*) inputFile) ;
	i = 0 ; more = 1 ;
	zold = -1.0 ; vold = 0.0 ;
	do {
		if( line != fgets(line,300,ifil)) 
			rLog(0,"End of file reading from %s", (void*)inputFile);
		j = sscanf(line,"%lf %lf ", z + i,v + i  );
/*		printf("j=%d\n",j) ; */
		if( j != 2 ) rLog(0,"Missing data reading: %s", (void*) line ) ;
/*		printf("%10.4f %10.4f %6d \n",v[i],z[i],i ) ; */
		more = (z[i] > zold) ;
		if( vold == v[i]) rLog(0,"Constant velocity is not supported",NULL) ;
		vold = v[i]; zold = z[i++] ; 
		if( i >=  MaxLayer ) more = 0 ;
	} while  ( more ) ;
	initVelModel( i, m ) ;
	for( j = 0 ; j < i ; j++ ) {
		*(m->z+j) = z[j] ;
		*(m->v+j) = v[j] ;
	}
}


double rtrace( double v1,double v2, double z, double p, double *x, double *t )

/* trace ray with rayparameter p through a layer of thickness z,
	v1	velocity at top
	v2	velocity at bottom
	z	thickness
	p	rayparameter
	x	horizontal distance
	t	traveltime
	returns	vertical distance from bottom where ray turns
		zero if the ray does not turn.
		returns -1.0 if velocity is to high for ray.

		Telford WM et. al. page 272-273
*/
{
	double si1, si2, co1, co2 ;
	double z0, r, b, b1 ;
	double zin, v2in,zturn ;
	si1 = p*v1 ;
	if( si1 >= 1.0 ) return( -1.0 ) ;
	si2 = p*v2 ;
	zturn = 0.0 ;
	if ( si2 >= 1.0 ) { 
		v2in = v2 ;
		v2 = 1.0 / p ;
		si2 = 1.0 ;
		zin = z ;
		z = zin * ( v2 - v1 ) / ( v2in - v1 ) ;
		zturn = zin - z ; 
	}
	co1 = sqrt( 1.0 - si1 * si1 ) ;
	co2 = sqrt( 1.0 - si2 * si2 ) ;
	b = ( v2 - v1 ) / z ;
	b1 = 1.0 / b ;
	z0 = v1 * b1 ;
	r = z0 / si1 ;
	*x = r * ( co1 - co2 ) ;
/*	*t = b1 * log( v1 * (1.0-co2)/(v2*(1.0-co1))) ;      */
	*t = b1 * log( si2 * (1.0+co1)/(si1*(1.0+co2))) ;
	return( zturn) ;
}
double traceModel( double p , double z0, VelModel *m, double *time ) 
{
	int i ;
	double x, v1,v2, z,zz,dz, xsum,t, tsum,zold ;
	xsum = 0.0 ; tsum = 0.0 ;
	i = 0 ;
	zold = 0.0 ;
	do {
		z = m->z[i] ;
		v1 = m->v[i] ;
		v2 = m->v[i+1] ;
		dz = z - zold ;
		if( i == m->nVel ) {
			
		}
		zz = rtrace(v1,v2,dz,p,&x,&t) ;
		xsum += x ;  tsum += t ;
	} while (zz == 0.0) ;
	*time = tsum+tsum ;
	return xsum+xsum ;
}

void testVel()
{
	VelModel m ;
	readVelModel("test.vel",&m) ;
	printVelModel(&m) ;
}
void testRay()
{
	int i ;
	double v1,v2,z,p,x,t,zturn ;
	for ( i = 60 ; i < 85 ; i++ ) {
		p = 0.006*(i+1) ;
		zturn = rtrace(2.0,2.5,6.0,p,&x,&t) ;
		printf(" %10.4f %10.4f %10.4f %10.4f\n",p,x,t,zturn) ;
	}
}

int main(int ac, char **av)
{
	int i ;
	testVel() ;
	testRay() ;
	return(0) ;
}

