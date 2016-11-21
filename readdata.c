
#include <stdlib.h> 
#include <stdio.h>
#include <ctype.h>

#define SIZE 510
char buffer[SIZE] ;
char *fields[SIZE] ;
int readData( FILE *fd, char **fields) 
/*	Read line from fd and put a list of non-blank fields into array fields. 
	Returns number of fields found
*/
{
	int n,nn,i,type,otype ;
	char *cp ;
	nn = SIZE ;
	cp = fgets(buffer,SIZE,fd) ;
	otype = isspace(' ') ;
	i = 0 ;
	while( *cp ) {
		type = isspace(*cp) ;
		if(( type == 0) && otype) { fields[i++] = cp ;	}
		otype = type ;
		cp++ ;
	}
	return  i ;
}
#ifdef TEST
int main()
{
	FILE *ff ;
	int n ;
	char *fields[200] ;
	double lat ;
	ff = fopen("../geysir/phase.dat","r") ;
	n = readData(ff,(char **) &fields) ;
	lat = atof(fields[7]) ;
	printf("n=%d lat=%lf _%s_\n",n,lat,fields[7]) ;

}
#endif
