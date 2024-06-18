# SIMFOR

## Зависимости
1) [Boost 1.84.0](https://www.boost.org/users/history/version_1_84_0.html) (нужно собрать вместе с модулем [MPI](https://www.boost.org/doc/libs/master/doc/html/mpi/getting_started.html))
2) OpenMP
3) MPI
4) GLUT 
5) OpenGL
6) \*для выполнения установки советую использовать утилиту checkinstall 

## Установка всех зависимостей
```shell
sudo apt update
sudo apt upgrade
sudo apt install wget make cmake build-essential g++ checkinstall libopenmpi-dev mesa-utils freeglut3-dev

wget https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.gz
tar -xvf boost_1_84_0.tar.gz
cd boost_1_84_0/
./bootstrap.sh
echo "using mpi ;" >> project-config.jam
sudo ./b2 install
```

## Установка SIMFOR
Установку SIMFOR можно выполнить с помощью скрипта `install-proj`, либо вручную.
После установки её можно добавить в проект следующим образом:
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

## Добавление в проект как сабмодуль
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
cd examples # тут хранятся примеры по папкам, в консоли их можно позапускать
```
