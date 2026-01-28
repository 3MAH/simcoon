:: conda install -c conda-forge git cmake python=3.9 -y

set PREFIX=%CONDA_PREFIX%
set SRC_DIR=%cd%
set SP_DIR=%PREFIX%\Lib\site-packages

:: Armadillo
git clone https://gitlab.com/conradsnicta/armadillo-code.git armadillo
cd armadillo
copy /y .\examples\lib_win64\libopenblas.dll %CONDA_PREFIX%\Library\bin\
copy /y .\examples\lib_win64\libopenblas.lib %CONDA_PREFIX%\Library\lib\
cmake -S . -B build -G "Visual Studio 17 2022" -DBLAS_LIBRARY:FILEPATH="%CONDA_PREFIX%\Library\lib\libopenblas.lib" -DLAPACK_LIBRARY:FILEPATH="%CONDA_PREFIX%\Library\lib\libopenblas.lib"
cmake --build build --config RELEASE
cmake --install build --prefix %CONDA_PREFIX%\Library
if errorlevel 1 exit 1

:: Simcoon
cd %SRC_DIR%
cmake -S . -B build ^
      -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -Wno-dev
cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

:: Arma2numpy
cd arma2numpy-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev 
cmake --build build --target ALL_BUILD --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cd %SRC_DIR%\simcoon-python-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev
cmake --build build --target ALL_BUILD --config Release
cmake --install build
if errorlevel 1 exit 1

cd %SP_DIR%

:: Install simcoon python 
mkdir simcoon
cd simcoon
type NUL > __init__.py
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\lib\Release\_core.pyd %SP_DIR%\simcoon
if errorlevel 1 exit 1

REM xcopy /s /i %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\arma2numpy-builder\build\bin\Release\arma2numpy.dll %PREFIX%\Library\bin
if errorlevel 1 exit 1