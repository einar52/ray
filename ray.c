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
		vold = v[i]; zold = z[i] ; 
		if( i >=  MaxLayer ) more = 0 ;
		if ( more ) i++ ;
	} while  ( more ) ;
	initVelModel( i, m ) ;
	for( j = 0 ; j < i ; j++ ) {
		*(m->z+j) = z[j] ;
		*(m->v+j) = v[j] ;
	}
	*(m->z) = 0.0 ;  /* depth of first velocity value must be 0.0 */
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
	printf("v1=%10.3f v2=%10.3f dz=%10.3f p=%10.3f x=%10.3f t=%10.3f\n",v1,v2,z,p,*x,*t) ;
	return( zturn) ;
}

double traceModel( double p , double zSource, VelModel *m, double *time ) 
{
	int i,more ;
	double x, v, z,zMax,dz,zz, xsum,t, tsum,vold,zold, vBottom ;
	double tSource,xSource,vSource, xtotal ;
	tSource = 0.0 ; xSource = 0.0 ;
	vBottom = 1.0/p ;
	zMax = m->z[(m->nVel)-1] ;
	if( zSource <= 0.0 ) zSource = zMax+1.0 ;
	xsum = 0.0 ; tsum = 0.0 ;
	zold = m->z[0] ;
	vold = m->v[0] ;
	i = 1 ;
	zold = 0.0 ;
	do {
		z = m->z[i] ;
		v = m->v[i] ;
		if( z >= zSource ) {
			vSource = vold + (zSource - zold)*(v-vold)/(z-zold) ;
			zz = rtrace(vold,vSource,zSource-zold,p,&x,&t ) ;
		printf("t =%12.5f zSource =%12.5f zz =%12.5f x =%12.5f \n",t,zSource,zz,x) ;
			xSource = x + xsum ;
			tSource = t + tsum ;
			zSource = zMax+1.0 ;
		}
		if( (i+1) == (m->nVel )) { /* if ray goes below last layer, extrapolate */
			if( vBottom > v ) {
				z = zold + (vBottom-vold)*(z-zold)/(v-vold) ;
				v = vBottom ;
			}
		}
		dz = z - zold ;
		zz = rtrace(vold,v,dz,p,&x,&t) ;
		printf("t =%12.5f z =%12.5f zz =%12.5f x =%12.5f \n",t,z,zz,x) ;
		vold = v ;
		zold = z ;
		xsum += x ;  tsum += t ;
		i++ ;
		more = (i < (m->nVel)) && ( zz == 0.0 ) ;
	} while (more) ;
	*time = tsum+tsum-tSource ;
	xtotal = xsum+xsum - xSource ;
	printf("tsum=%12.6f xsum=%12.6f\n",*time,xtotal) ;
	return xtotal ;
}
double velZ( double z , VelModel *m, int *iLayer)
{
	int i ;
	double v0,v1,z0,z1 ;
	for ( i = 0 ; i < m->nVel ; ) {
		if ( z < m->z[i] ) { i++ ; break ; } 
		i++ ;
	}
	i-- ;
	*iLayer = i ;
	v0 = m->v[i-1] ;
	v1 = m->v[i] ;
	z0 = m->z[i-1] ;
	z1 = m->z[i] ;
	printf(" %g %g %g %g \n",z0,z1,v0,v1) ;
	return v0 + (z-z0)*( v1-v0)/(z1-z0) ;
}
double traceUp( double p, double zSource, VelModel *m, double *tTime)
{
	int iLayer,i ;
	double vSource,zold,vold,z,v,xResult,x,t,zz ;
	vSource = velZ( zSource, m, &iLayer) ;
	xResult = 0.0 ; *tTime = 0.0 ;
	zold = zSource ;
	vold = vSource ;
	for( i = iLayer ; i > 0 ; ) {
		i-- ;
		z = m->z[i];
		v = m->v[i];
		zz = rtrace( v,vold,zold-z,p,&x,&t ) ;
		*tTime += t ;
		xResult += x ;
		zold = z ; vold = v ;
	}
	return( xResult ) ;
}
void testZ( VelModel *m)
{
	double z,v ;
	int i ;
	z = 3850 ;
	z = 5.99 ;
	z = 0.4 ;
	v = velZ(z, m,&i ) ;
	printf("z=%12.6f v = %12.6f i=%d\n",z,v,i ) ;
}
void testVel()
{
	VelModel m ;
	int n ;
	double vmax,p,t,x,z ;
	readVelModel("test.vel",&m) ;
	printVelModel(&m) ;
	n = m.nVel ;
	vmax = m.v[n-1] ;
	p = 1.01/vmax ;
	p = 1.00/7.4 ;
	z = 5.99 ;
	z = 0.4 ;
        x = traceModel(p,z,&m,&t) ;
	x = traceUp(p,z,&m,&t) ;
	testZ( &m ) ;
}
void testRay()
{
	int i ;
	double v1,v2,z,p,x,t,zturn ;
	for ( i = 75 ; i < 80 ; i++ ) {
		p = 0.006*(i+1) ;
		zturn = rtrace(2.0,2.5,6.0,p,&x,&t) ;
		printf(" %10.4f %10.4f %10.4f %10.4f\n",p,x,t,zturn) ;
	}
}
int main(int ac, char **av)
{
	int i ;
	testVel() ;
/*	testRay() ;  */
	return(0) ;
}

