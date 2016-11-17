#!/usr/bin/python

#print "hello world"
station = []
lat = []
lon = []
def readStations() :
    ff = open("../hau99/net.dat")
    fflines = ff.readlines()
    for l in fflines :
        ll = l.split()
        if( len(ll) == 4 ) :
            station.append(ll[0])
            lat.append( float(ll[1])) 
            lon.append( float(ll[2]))
    ff.close()
import math
earthRadius = 6366.8
earthRadius = 6399.0
earthRadius = 6394.0

deg2rad = math.pi/180.0
def distSphere( la1,la2,dlon ) :
    b1 = la1 * deg2rad
    b2 = la2 * deg2rad
    dlo = dlon * deg2rad
    x1 = math.cos(b1) 
#   y1 = 0.0
    z1 = math.sin(b1)
    co2 = math.cos(b2) 
    x2 = co2*math.cos(dlo)
    y2 = co2*math.sin(dlo)
    z2 = math.sin(b2)
    xp = -z1*y2
    yp = z1*x2 - x1*z2
    zp = x1*y2 
    sind = math.sqrt(xp*xp + yp*yp + zp*zp)
    cosd = x1*x2  +  z1*z2
    dist = earthRadius * math.atan2(sind,cosd) 
    return(dist)
def main() :
    file = open("../hau99/phase_good_1999.tbl")
    lines = file.readlines()
    n = 0
    for x in lines :
       l= x.split()
       phase = l[2]
       if( n > 1 ) :
         time = float(l[5]) - float(l[15])
    #     if( time < 0 ) :
         if( l[4] != l[14] ) :
            time = time + 60.0 
	 sta = l[1] 
	 latq = float(l[17])
	 lonq = float(l[19])
         j = station.index(sta)
	 distq = distSphere(latq,lat[j],lonq-lon[j] )
         dist = float(l[8])
         res = float(l[6]) 
         depth = float(l[24]) 
         if ( (phase == 'p') & ( depth > 2.5) ) :
#            print  '%6.1f %6.1f %6.3f %6.2f' % ( distq,depth,time, res),l[0],l[1],l[2]
            print  '%6.1f %6.1f %6.3f %6.2f' % ( distq,depth,time, res),l[0],l[1],l[2],'%6.1f %4.1f'%(distq,dist-distq)
       n = n + 1
#print "done"
readStations()
main()
