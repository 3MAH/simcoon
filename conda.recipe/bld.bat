:: Simcoon unified build with Python bindings
cd %SRC_DIR%

:: Set pybind11 path for CMake
python -c "import pybind11; print(pybind11.get_cmake_dir())" > %TEMP%\_pb11.txt
set /p PYBIND11_DIR=<%TEMP%\_pb11.txt
del %TEMP%\_pb11.txt
if "%PYBIND11_DIR%"=="" (
  echo ERROR: Could not find pybind11 CMake directory
  exit 1
)
echo Using pybind11 from: %PYBIND11_DIR%

:: Configure with Python bindings enabled
:: CONDA_BUILD env var is detected by CMakeLists.txt to enable SIMCOON_BUILD_PYTHON
:: Use Ninja (already a build requirement) with the conda-activated MSVC (CC/CXX=cl.exe) rather
:: than a hardcoded Visual Studio generator. The VS generator fails on CI when the activated VS
:: edition does not match the hardcoded version ("could not find any instance of Visual Studio").
cmake -G Ninja ^
      -S . -B build ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DSIMCOON_BUILD_TESTS:BOOL=OFF ^
      -DPython3_EXECUTABLE:FILEPATH="%PYTHON%" ^
      -Dpybind11_DIR:PATH="%PYBIND11_DIR%" ^
      -DCMAKE_PREFIX_PATH="%PYBIND11_DIR%"
if errorlevel 1 exit 1

:: Build everything (C++ library and Python bindings)
cmake --build build --config Release
if errorlevel 1 exit 1

:: Install C++ library to %PREFIX%/Library (lib, include, cmake)
:: and Python bindings to site-packages
cmake --install build --config Release
if errorlevel 1 exit 1
