CC=g++-5
CFLAGS= -std=c++11 -g -O3 -larmadillo -llapack -lblas 
OBJ =vib.o 
DEPS = mass.h 

all: vib.exe

vib.exe: $(OBJ)
	$(CC) -o vib.exe $^ $(CFLAGS)

%.o: %.cpp $(DEP)
	$(CC) -c $< $(CFLAGS)

clean:
	rm *.o *.exe
