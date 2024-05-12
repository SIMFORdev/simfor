CC=g++
CFLAGS=-Wall -Wextra
LDFLAGS=-L. -lplotter -lGL -lglut
SRCDIR=src
OBJDIR=build
EXAMDIR=examples

# Цели
all: $(OBJDIR) libplotter.a examples/test_2d examples/test_3d examples/test_points

# Каталог build
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Библиотека
libplotter.a: $(OBJDIR)/plotter.o
	ar rcs $@ $^

# Объектные файлы
$(OBJDIR)/plotter.o: $(SRCDIR)/plotter.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

$(OBJDIR)/test_2d.o: $(EXAMDIR)/test_2d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

$(OBJDIR)/test_3d.o: $(EXAMDIR)/test_3d.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

$(OBJDIR)/test_points.o: $(EXAMDIR)/test_points.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

# Исполняемые файлы
$(EXAMDIR)/test_2d: $(OBJDIR)/test_2d.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(EXAMDIR)/test_3d: $(OBJDIR)/test_3d.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(EXAMDIR)/test_points: $(OBJDIR)/test_points.o libplotter.a
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Очистка
clean:
	rm -f *.a $(EXAMDIR)/test_2d $(EXAMDIR)/test_3d $(EXAMDIR)/test_points
	rm -rf $(OBJDIR)
