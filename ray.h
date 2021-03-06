
extern int shLogLevel ;
extern double zBottom ;
typedef struct {
	int nVel ;
	double *z ; /* depth  */
	double *v ; /* velocity */
} VelModel ;

typedef struct {
	char name[4] ;
	double lat,lon,depth ;
} Station ;

typedef struct {
	long long index ;
	Station *statP ;
	double pTime ;
	double weight ;
	char type ;
} Phase ;
typedef struct {
	long long index ;
	double lat,lon,depth, timeShift ;
} Solution ;
typedef struct { double x,z ; } DepthPoint ;

#define SURFACE	1
#define RayUP	2
#define RayDown 3
#define MaxLayer 1000
void  rLog( int level, char *s1 , void *p );
void initVelModel( int nVel , VelModel *m );
void printVelModel( VelModel *m );
void readVelModel( char *inputFile, VelModel *m );
double rtrace( double v1,double v2, double z, double p, double *x, double *t );
double traceModel( double p , double zSource, VelModel *m, double *time ) ;
double velZ( double z , VelModel *m, int *iLayer);
double traceUD( int mode, double p, double zSource, VelModel *m, double *tTime);
VelModel resampleVelModel( VelModel *mIn, double dz, int nz )  ;
VelModel resampleVelModel2( VelModel *mIn, double dz, int nz )  ;
/* resample velocity model using 2 point interpolation */
void splineTest() ;
double timeFromDist( VelModel *m, double x, double z, double *p, double *dtdx, double *dxdp) ;

/*stations.c */

Station *lookUpStation( char *name )  ;

/* distance.c */
double gDistance(double la1, double la2, double dlon);

double sDistance(double la1, double la2, double dlon);

/* readdata.c */
int readData( FILE *fd, char **fields) ;
/* phases.c */
int readPhases(char *fileName, Phase **phases )  ;
int readReloc(char *fileName, Solution **solutions) ;
int readCtloc(char *fileName, Solution **solutions ) ;
void printPhases( Phase *p, int n  ) ;
void printSol( Solution *p, int n ) ;

void  golubC(double *a, double *x, double *b, int m, int n) ;
