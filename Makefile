
CFLAGS=-g

tt : travelt 
	./travelt -d 3

T = travelt.o ray.o

$T : ray.h

travelt : $T
	cc ${CFLAGS}  -o travelt $T -lm

t : ray
	./ray

ray : ray.c  ray.h
	cc -g -DDEBUG -o ray ray.c -lm

clean : 
	rm *.o travelt ray


