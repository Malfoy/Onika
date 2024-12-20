CC = g++

CFLAGS= -O3 -g -std=c++17 -Wall  -w -fopenmp -flto
CFLAGS_DELTA= -O3 -std=c++17 -Wall -w -flto
#LDFLAGS = -LTurboPFor-Integer-Compression -l:libic.a -pthread  -lpthread -lz -fopenmp -flto
LDFLAGS = -pthread -lpthread -lz -fopenmp -flto


#PREXEC = turboPFor_compile
EXEC = onika
OBJ = utils.o onika_index.o 
OBJ_INDEX = onika.o onika_index.o utils.o #TurboPFor-Integer-Compression/libic.a

all :  $(EXEC)

# TurboPFor-Integer-Compression/libic.a:
# 	make -C TurboPFor-Integer-Compression/

onika : $(OBJ_INDEX)
	$(CC) -o onika $^ $(LDFLAGS)

onika.o: src/onika.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: src/utils.cpp headers/utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

onika_index.o: src/onika_index.cpp headers/onika_index.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) $(OBJ) $(OBJ_INDEX) .vscode onikaOutput.gz

rebuild: clean $(EXEC)