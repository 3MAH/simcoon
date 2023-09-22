cmake -G"Visual Studio 17 2022" ^
      -D CMAKE_BUILD_TYPE:STRING=Release ^
      -D CMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -D CARMA_INSTALL_LIB=ON ^
      -S . -B build

cmake --build build --config Release
cmake --install build