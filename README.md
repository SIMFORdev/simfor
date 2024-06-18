# SIMFOR

## Зависимости
1) [Boost 1.84.0](https://www.boost.org/users/history/version_1_84_0.html) (нужно собрать вместе с модулем [MPI](https://www.boost.org/doc/libs/master/doc/html/mpi/getting_started.html))
2) OpenMP
3) MPI
4) GLUT 
5) OpenGL
6) \*для выполнения установки советую использовать утилиту checkinstall

## Добавление в проект
```cmake
cmake_minimum_required(VERSION 3.23)
project(TestProject)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(simfor)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC simfor)
```

## Запуск примеров
```shell
mkdir build
cd build
cmake ..
make
# Пусть наша папка в которой мы работаем называется examples/SomeClassTest.
cd examples
```

## Установка
Установку библиотеки можно выполнить с помощью скрипта `install-proj`, либо вручную.
```shell
mkdir build && cd build 
cmake ..
make -j 8
sudo checkinstall
```
После установки добавление в проект
```cmake
cmake_minimum_required(VERSION 3.23)
project(TestProject)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost 1.84.0 COMPONENTS mpi REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(simfor REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC
        simfor
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARY}
        OpenMP::OpenMP_CXX
)
```