CC=g++-4.8
EIGEN=/home/lousydropout/Software/eigen/
CFLAGS=-Wall -O3 -std=c++0x
all: feedback fsp

feedback: feedback.cpp
	$(CC) $(CFLAGS) -I $(EIGEN) feedback.cpp -o feedback

fsp: fsp.cpp
	$(CC) $(CFLAGS) -I $(EIGEN) fsp.cpp -o fsp

clean:
	rm -f *~ fsp feedback *.o

