CC = gcc
CFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic -pthread
LDLIBS = -pthread
all: cache_lab
cache_lab: main.c
	$(CC) $(CFLAGS) main.c -o $@ $(LDLIBS)
quick: cache_lab
	mkdir -p results plots
	./cache_lab --quick
	python3 plot_results.py
run: cache_lab
	mkdir -p results plots
	./cache_lab
	python3 plot_results.py
clean:
	rm -f cache_lab
.PHONY: all quick run clean
