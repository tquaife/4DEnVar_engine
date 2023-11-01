CC = gcc
CFLAGS = -Wall -g -std=c99
LIBS = -lgsl -lgslcblas -lm 
INCLS = -I./ -I/usr/include/gsl


4DEnVar:	4DEnVar.o 4DEnVar_engine.o matrixio.o
		${CC} ${CFLAGS} $@.o -o $@ matrixio.o 4DEnVar_engine.o ${INCLS} ${LIBS}

4DEnVar_surf:	4DEnVar_surf.o 4DEnVar_engine.o matrixio.o
		${CC} ${CFLAGS} $@.o -o $@ matrixio.o 4DEnVar_engine.o ${INCLS} ${LIBS}

.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<

clean:
		\rm -f *.o *~ *%
