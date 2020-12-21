.ONESHELL:
CC=gcc
CFLAGS= -O2 -lm -lgsl
DEPS=
OBJS=main.o
OUT=main

$(OUT): $(OBJS)
	$(CC) -o $@  $^ $(CFLAGS)

.c.o:
	$(CC) -c -o $@ $< $(CFLAGS)

run: $(OUT)
	./$(OUT)

all:
	make clean
	make run

memcheck: $(OUT)
	valgrind --leak-check=full ./$(OUT)
plot: $(OUT)
	gnuplot gnuplot.sh

clean:
	rm -f *.o
	rm -f *.dat
	rm -f *~
	rm -f *.eps
	rm -f $(OUT)

.PHONY: clean plot run
