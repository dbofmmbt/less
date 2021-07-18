CC = mpicc
SRC = src/*.c

check:
	${CC} ${SRC} -fsyntax-only -Wall

release:
	mkdir -p build/release
	${CC} ${SRC} -Wall -o build/release/main

run: release
	mpirun -n $(p) ./build/release/main