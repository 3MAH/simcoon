:: Simcoon unified build with Python bindings
cd %SRC_DIR%

:: Set pybind11 path for CMake
for /f "delims=" %%i in ('python -c "import pybind11; print(pybind11.get_cmake_dir())"') do set PYBIND11_DIR=%%i
if "%PYBIND11_DIR%"=="" (
  echo ERROR: Could not find pybind11 CMake directory
  exit 1
)
echo Using pybind11 from: %PYBIND11_DIR%

:: Configure with Python bindings enabled
cmake -G"Visual Studio 17 2022" ^
      -S . -B build ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DSIMCOON_BUILD_PYTHON_BINDINGS:BOOL=ON ^
      -DSIMCOON_BUILD_TESTS:BOOL=ON ^
      -DPython3_EXECUTABLE:FILEPATH="%PYTHON%" ^
      -Dpybind11_DIR:PATH="%PYBIND11_DIR%"
if errorlevel 1 exit 1

:: Build everything (C++ library and Python bindings)
cmake --build build --config Release
if errorlevel 1 exit 1

:: Install C++ library to conda prefix
cmake --install build --config Release
if errorlevel 1 exit 1

:: Install Python package components to build/python-package
cmake --install build --config Release --component python
if errorlevel 1 exit 1

:: Install Python package from the generated package directory
python -m pip install %SRC_DIR%/build/python-package --no-deps -vv
if errorlevel 1 exit 1
