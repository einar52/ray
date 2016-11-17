
CFLAGS=-g

top : t14

tt : travelt 
	./travelt -d 3 -x 15

T = travelt.o ray.o

$T : ray.h

travelt : $T
	cc ${CFLAGS}  -o travelt $T -lm

t : ray
	./ray

ray : ray.c  ray.h
	cc -g -DDEBUG -o ray ray.c -lm

t22 : travelt
	travelt -n 0 -l 9 -d 7.3 -x 22.1 > jt22  2>&1

t14 : travelt
	travelt -n 0 -l 9  -x 14.6 -d 5.6 > jtest  2>&1

clean : 
	rm *.o travelt ray


