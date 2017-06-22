OBJ = mycorr.o sacio.o
mycorr : $(OBJ)
	cc -o mycorr $(OBJ) -lm -lfftw3

$(OBJ) : sacio.h

clean : 
	rm -f mycorr $(OBJ)
