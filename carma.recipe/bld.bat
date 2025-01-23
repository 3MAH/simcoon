cmake -S . -B build ^
      -D CMAKE_BUILD_TYPE:STRING=Release ^
      -D CMAKE_INSTALL_PREFIX=%PREFIX%/Library

cmake --build build --config Release
cmake --install build
