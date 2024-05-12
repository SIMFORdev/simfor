# SIMFOR_plotter_c

## Установка
1. Установите требуемые зависимости:
```
sudo apt-get install freeglut3 freeglut3-dev
```

2. Соберите библиотеку и тесты:

* Makefile
```
make
```

* или CMake
```
mkdir build
```
cd build
```
```
cmake ..
```
```
make
```

* или сборка статической библиотеки:
```
g++ -c plotter.cpp -o plotter.o
```
```
ar rcs libplotter.a plotter.o
```
Компиляция вашей программы test:
```
g++ test.cpp -o test -L. -lplotter -lGL -lglut
```

Makefile собирает статическую библиотеку libplotter.a и примеры test_2d, test_3d, test_points.

**Управление (зависит от раскладки и регистра)**
* "o" - циклический выбор отображаемого графика;
* "i" - вкл/выкл показа всех графиков;
* "g" - циклическая смена режима отображения графиков:
Для графика: линия -> точка -> точки и линии;
Для поверхности: точки-> линии с точками -> линии -> заливка с слабозаметными линиями;

* "x","y","z" - поворот вокруг соответсвующей оси;
* "X","Y","Z" - поворот в обратную сторону;

* "w","s" - перемещение по оси У;
* "a","d" - перемещение по оси X;
* "q","e" - перемещение по оси Z;

* "l" - показать легенду;

* "n" - сетка;
* "m" - миллиметровый режим.


