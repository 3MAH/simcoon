cmake -G Ninja ^
      -D CMAKE_BUILD_TYPE:STRING=Release ^
      -D CMAKE_INSTALL_PREFIX:PATH=%PREFIX%
      -D CARMA_INSTALL_LIB=ON
      -S . -B build ^      

cmake --build build --target ALL_BUILD --config Release
cmake --install build