CC = gcc
SRC = src/main.c

check:
	${CC} ${SRC} -fsyntax-only -Wall

release:
	mkdir -p build/release
	${CC} ${SRC} -Wall -o build/release/main

run: release
	./build/release/main