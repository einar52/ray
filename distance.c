/*  compute distance between two points on a sphere (sDistance) or an
    ellipsiod (gDistance) .
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <geodesic.h>

#define A 6378.137
#define F 1/298.257223563 /* WGS84 */
#define deg2rad M_PI/180.0
#define earthRadius 6391.57

struct geod_geodesic distance_g  ;
int distance_I = 1 ;
double gDistance(double la1, double la2, double dlon)
{
	double dist,az1,az2 ;
	if( distance_I ) { geod_init(&distance_g,A,F) ;
				distance_I = 0 ;
	}
	geod_inverse(&distance_g,la1,0.0,la2,dlon,&dist,&az1,&az2) ;
	return dist ;
}
double sDistance(double la1, double la2, double dlon)
{
    double b1,b2,dlo,x1,z1,co2,x2,y2,z2,xp,yp,zp,sind,cosd,dist ;
    b1 = la1 * deg2rad ;
    b2 = la2 * deg2rad ;
    dlo = dlon * deg2rad ;
    x1 = cos(b1)  ;
    z1 = sin(b1) ;
    co2 = cos(b2)  ;
    x2 = co2*cos(dlo) ;
    y2 = co2*sin(dlo) ;
    z2 = sin(b2) ;
    xp = -z1*y2 ;
    yp = z1*x2 - x1*z2 ;
    zp = x1*y2  ;
    sind = sqrt(xp*xp + yp*yp + zp*zp) ;
    cosd = x1*x2  +  z1*z2 ;
    dist = earthRadius * atan2(sind,cosd)  ;
    return(dist) ;
}

#ifdef TEST
speedTest()
{
	double  la1,la2,dist,dlon ;
	int i,j ;
	dlon = 100 ;
   	for( i = 0 ; i < 2000 ; i++ ) {
		la1 = 10 + 0.01*i ;
		for( j = 0 ; j < 5000 ; j++ ) {
			la2 = 70 + 0.01*j ;
			dist = gDistance(la1,la2,dlon) ;
		}
	} 
	/*  time for 10M calculations is about 7 seconds for gDistance and
		1.25 seconds for sDistance */
}
main()
{
	float la1,la2,dlon ;
	speedTest() ;
	la1 = 64 ; la2 = 65 ; dlon = 2.3 ;
	la1 = 64 ; la2 = 66 ; dlon = 4.6 ;
	float distg,dists,distx, a1,a2 ;
	distg = gDistance(la1,la1,dlon) ;
	dists = gDistance(la1,la2,0) ;
	distx = gDistance(la1,la2,dlon) ;
	printf("distg=%10.3f dists=%10.3f distx=%10.3f\n",distg,dists,distx) ;
	distg = sDistance(la1,la1,dlon) ;
	dists = sDistance(la1,la2,0) ;
	distx = sDistance(la1,la2,dlon) ;
	printf("distg=%10.3f dists=%10.3f distx=%10.3f\n",distg,dists,distx) ;
}
#endif
