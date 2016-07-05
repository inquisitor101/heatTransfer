CC = mpicc
LD = mpicc
CFLAGS = -std=c99 -Wall -Werror -O0 # change -O0 to -O3
LDFLAGS = -lm -lmpi
RM = /bin/rm -f
OBJS = main.o funcs.o util.o
EXECUTABLE = heat3D

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS)

util.o: header.h util.h util.c
	$(CC) $(CFLAGS) -c util.c

funcs.o: header.h funcs.h funcs.c
	$(CC) $(CFLAGS) -c funcs.c

main.o: main.c funcs.h header.h util.h
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)
