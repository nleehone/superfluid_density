CC = g++
CFLAGS = -O3 -fPIC
all: rhosf logging
	gcc -shared MMA.cc rhosf.o logging.o -o librhosf.so -I/usr/local/Wolfram/Mathematica/10.3/SystemFiles/IncludeFiles/C/ -fPIC -lgsl -lgslcblas
	cp librhosf.so ../

rhosf: rhosf.cpp logging.o
	$(CC) $(CFLAGS) rhosf.cpp -c -std=c++11

logging: logging.cpp
	$(CC) $(CFLAGS) logging.cpp -c -std=c++11

.PHONY : clean
clean:
	rm *.o

