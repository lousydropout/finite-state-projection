CC=g++-4.8
EIGEN=/home/lousydropout/Software/eigen/
CFLAGS=-Wall -Wextra -O3 -std=c++0x
all: fsp

# fsp_package.o: fsp_package.h
# 	$(CC) $(CFLAGS) -I $(EIGEN) -c fsp_package.cpp
fsp: fsp.cpp
	$(CC) $(CFLAGS) -I $(EIGEN) fsp.cpp -o fsp

clean:
	rm -f *~ fsp *.o

