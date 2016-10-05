#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ray.h"

double zBottom ; 
int shLogLevel = 3 ;
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
	FILE *out ;
	out = fopen("dump.vel","w") ;
	for( i = 0 ; i < m->nVel ; i++) 
		fprintf(out,"%10.4f %10.4f %6d \n",*(m->z+i),*(m->v+i),i ) ;
	fclose(out) ;
}
double spline2(double x, double y, double yy )
{
	double x2,x3,res ;
	x2 = x*x ;
	x3 = x2*x ;
	res = y * ( 1.0 - 3.0 * x2 + 2.0 * x3 ) ;
	res += yy * (  x - x2 - x2 + x3 ) ;
	return(res) ;
}
void splineTest()
{
	double x ;
	int n,i ;
	n = 200 ;
	for(i = 0 ; i <= n ; i++) {
		x = (i * 1.0) / n ;
		printf("%6d %10.4f %10.4f  %10.4f \n",i,x, spline2(x,1.0,0.0), spline2(x,0.0,1.0) );
	}
}
VelModel resampleVelModel( VelModel *mIn, double dz, int nz ) 
{
	VelModel m ;
	initVelModel(nz+mIn->nVel,&m) ;
	double z0,z1,z2,z3, v0,v1,v2,v3,z ;
	double x,s1,s2 ;
	int i,j ;
	z1 = mIn->z[0] ;
	z2 = mIn->z[1] ;
	z3 = mIn->z[2] ;
	v1 = mIn->v[0] ;
	v2 = mIn->v[1] ;
	v3 = mIn->v[2] ;
	z0 = z1+z1-z2 ;
	v0 = v1+v1-v2 ;
	s1 = (z2-z1)*(v2-v0)/(z2-z0) ;
	s2 = (z2-z1)*(v1-v3)/(z3-z1) ;
	z = z1 ;
	j = 3 ;
	for( i = 0 ; i < nz ; i++) {
		x = (z-z1)/(z2-z1) ;
		m.z[i] = z ;
		m.v[i] = spline2(x,v1,s1) + spline2(1.0-x,v2,s2) ;
		z += dz ;
		while( z >= z2 ) {
			z0 = z1 ; v0 = v1 ;
			z1 = z2 ; v1 = v2 ;
			z2 = z3 ; v2 = v3 ;
			z3 = mIn->z[j] ; v3 = mIn->v[j++] ;
			s1 = (z2-z1)*(v2-v0)/(z2-z0) ;
			s2 = (z2-z1)*(v1-v3)/(z3-z1) ;
			if( j >= mIn->nVel ) {
				m.nVel = i ;
				return(m) ;
			}
		}
	}
	j -= 2 ;
	while ( j < mIn->nVel ) {
		m.z[i] = mIn->z[j] ;
		m.v[i++] = mIn->v[j++] ;
	}
	m.nVel = i ; 
	return m ;
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
	if( z <= 0.0 ) 	{ *t = 0.0, *x =0.0 ; return(0.0) ; }
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
	if( shLogLevel < 4 ) return( zturn ) ;
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
	return v0 + (z-z0)*( v1-v0)/(z1-z0) ;
}
double traceUD( int mode, double p, double zSource, VelModel *m, double *tTime)
{
	int iLayer,i ;
	double vSource,zold,vold,z,v,xResult,x,t,zz ;
	double timeUp,xUp, timeDown,xDown ;
	zBottom = 0.0 ;
	xUp = 0.0 ; timeUp = 0.0 ;
	if( mode != SURFACE ) {
		vSource = velZ( zSource, m, &iLayer) ;
		zold = zSource ;
		vold = vSource ;
		for( i = iLayer ; i > 0 ; ) {
			i-- ;
			z = m->z[i];
			v = m->v[i];
			zz = rtrace( v,vold,zold-z,p,&x,&t ) ;
			zold = z ; vold = v ;
			timeUp += t ;
			xUp += x ;
		}
		if( mode == RayUP ) {
			*tTime = timeUp ;
			return( xUp ) ;
		}
		zold = zSource ;
		vold = vSource ;
	} else {
		zold = 0.0 ;
		vold = m->v[0] ;
		iLayer = 1 ;
	}
	timeDown = 0.0 ; xDown = 0.0 ;
	for( i = iLayer ; i < m->nVel ; i++) {
/*		printf("%d ",i) ; */
		z = m->z[i];
		v = m->v[i];
		zz = rtrace( vold,v,z-zold,p,&x,&t ) ;
		timeDown += t ;
		xDown    += x ;
		zold = z ; vold = v ;
		if( zz > 0.0 ) break ;
	}
	zBottom = zold - zz ;
	*tTime = timeUp + timeDown + timeDown ;
	return( xUp + xDown + xDown ) ;
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
	p = 1.00/5.4 ;
	z = 0.4 ;
	z = 5.99 ;
        x = traceModel(p,z,&m,&t) ;
/*	x = traceUp(p,z,&m,&t) ;
	x = traceDown(p,z,&m,&t) ;  */
	x = traceUD(RayUP,p,z,&m,&t) ;
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
void test2()
{
	VelModel m ;
	int i,n ;
	double pMax, depth,p,x,t,xd,td,xs,ts ;
	readVelModel("test.vel", &m) ;
	if( shLogLevel > 3 )	printVelModel(&m) ; 
	depth = 3.8 ;
	n = 5 ;
	pMax = 1.0/velZ(depth,&m,&i) ;
	for( i = 1 ; i < 2*n ; i++) {
		p = i*pMax/(2*n) ;
		x = traceUD(RayUP,p,depth,&m,&t) ;
		xd = traceUD(RayDown,p,depth,&m,&td) ;
		xs = traceUD(SURFACE,p,depth,&m,&ts) ;
		printf("p=%10.4f x=%8.4f %8.4f %8.4f t=%8.4f %8.4f %8.4f\n",
			p,x,xd,x+xd-xs,t,td,t+td-ts) ;
	}
}

#ifdef DEBUG
int main(int ac, char **av)
{
	int i ;
	shLogLevel = 3 ;
	test2() ;
/*	testVel() ; */
/*	testRay() ;  */
	return(0) ;
}
#endif 
