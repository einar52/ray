
extern int shLogLevel ;
extern double zBottom ;
typedef struct {
	int nVel ;
	double *z ; /* depth  */
	double *v ; /* velocity */
} VelModel ;
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
void splineTest() ;
double timeFromDist( VelModel *m, double x, double z, double *p, double *dtdx, double *dxdp) ;

