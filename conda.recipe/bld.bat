:: Simcoon build - C++ library via CMake, Python via scikit-build-core
cd %SRC_DIR%

:: Configure C++ library and tests (Python bindings built separately via pip)
cmake -G"Visual Studio 17 2022" ^
      -S . -B build ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DSIMCOON_BUILD_PYTHON_BINDINGS:BOOL=OFF ^
      -DSIMCOON_BUILD_TESTS:BOOL=ON
if errorlevel 1 exit 1

:: Build C++ library
cmake --build build --config Release
if errorlevel 1 exit 1

:: Install C++ library and headers to conda prefix
cmake --install build --config Release
if errorlevel 1 exit 1

:: Install Python package via scikit-build-core
cd %SRC_DIR%
python -m pip install . --no-deps --no-build-isolation -vv
if errorlevel 1 exit 1
