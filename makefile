CFLAGS = -O3 -Wall -Winline -Wshadow -std=c++11

mgsolve: ./mgsolver.cpp
	g++ $(CFLAGS) ./mgsolver.cpp -o mgsolve

clean:
	rm -rf *.o solution.txt mgsolve

.PHONY: clean
