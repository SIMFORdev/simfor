CC=g++
CFLAGS=-Wall -Wextra
LDFLAGS=-L. -lplotter -lGL -lglut

# Цели
all: libplotter.a test_2d test_3d test_points

# Библиотека
libplotter.a: plotter.o
	ar rcs $@ $^

# Объектные файлы
plotter.o: src/plotter.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_2d.o: examples/test_2d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_3d.o: examples/test_3d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

test_points.o: examples/test_points.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

# Исполняемые файлы
test_2d: test_2d.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test_3d: test_3d.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test_points: test_points.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Очистка
clean:
	rm -f *.o *.a test_2d test_3d test_points
