CC=g++
CFLAGS=-Wall -Wextra
LDFLAGS=-L. -lplotter -lGL -lglut

all: libplotter.a test_2d test_3d test_points

libplotter.a: plotter.o
	ar rcs $@ $^

plotter.o: plotter.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_2d: test_2d.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test_3d: test_3d.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test_points: test_points.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test_2d.o: test_2d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_3d.o: test_3d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_points.o: test_points.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f *.o *.a test_2d test_3d test_points
