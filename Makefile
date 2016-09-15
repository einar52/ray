
CFLAGS=-g

t : ray
	./ray

ray : ray.o
	cc -g -o ray ray.c -lm

clean : 
	rm ray.o ray


