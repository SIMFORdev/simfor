# SIMFOR

## Добавление нового модуля
1) Создать хедер и поместить его в include/simfor 
2) Создать файл сорцов и поместить его в src 
3) В CMakeLists.txt прописать путь до хедера в `SIMFOR_INCLUDES` и 
путь до исходного файла в `SIMFOR_SOURCES`.
4) В `examples` необходимо сделать свою папку, в которой будет лежать 
ваш исполняемый файл. Пусть папка будет называться `example`.
5) Там необходимо создать одноименный `.cpp` файл, например `example.cpp`,
и `CMakeLists.txt` файл. Его содержание скопировать от сюда:
```cmake
cmake_minimum_required(VERSION 3.23)
project(example)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} example.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC simfor)
```

Пример:

- Создали `SomeClass.hpp` и добавили его в `include/simfor`. Теперь его путь 
`include/simfor/SomeClass.hpp`
- Создали `SomeClass.cpp` и добавили его в `src`. Теперь его путь
  `include/simfor/SomeClass.cpp`.
- Обновили `CMakeLists.txt`.
```cmake
set(SIMFOR_INCLUDES
        # Public API includes
        include/simfor/SomeClass.hpp # <- Добавили тут
        # Private API includes
        include/simfor/internal/types.hpp
)

set(SIMFOR_SOURCES
        ${SIMFOR_INCLUDES}
        src/SomeClass.cpp            # <- Добавили тут
)
```
- Создаем папку и исполняемый файл в examples. Папку назовем `SomeClassTest`,
а файл соответственно `SomeClassExe.cpp`.
```c++
#include <iostream>
#include "simfor/SomeClass.hpp"

int main(int argc, char **argv) {
	std::cout << "Hello examples in SIMFOR\n";
	SomeClass testclass;
	testclass.HelloWorld();
	return 0;
}
```
- Создаем `CMakeLists.txt` и заполняем его, не забывая изменить `project` и 
`add_executable`.
```cmake
cmake_minimum_required(VERSION 3.23)
project(SomeClassTest)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} SomeClassTest.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC simfor)
```

В итоге получается следующая структура папок:

```text
examples/
    ...
    SomeClassTest/
        CMakeLists.txt         <- CMakeLists.txt от SomeClassTest
        SomeClassTest.cpp
    ...
include/
    simfor/
        ...
        SomeClass.hpp
        ...
src/
    ...
    SomeClass.cpp
    ...
CMakeLists.txt                 <- корневой CMakeLists.txt
```