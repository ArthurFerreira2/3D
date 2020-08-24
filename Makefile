LIBS = -lmingw32 -lSDL2main -lSDL2
FLAGS = -std=c99 -pedantic -Wpedantic -Wall -Werror -Wl,-subsystem,windows -O3
CC = gcc

3Drenderer: 3Drenderer.c
	$(CC) $^ $(FLAGS) $(LIBS) -o $@

all: 3Drenderer
