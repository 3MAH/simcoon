:: gmp for windows with channel -c salilab
:: conda install ninja
:: conda install -c conda-forge m2w64-gcc
:: conda install -c conda-forge fortran-compiler 
:: conda install -c conda-forge vs2019_win-64

:: Armadillo
git clone https://gitlab.com/conradsnicta/armadillo-code.git armadillo
if errorlevel 1 exit 1
cd armadillo
cmake -S . -B build -G "Visual Studio 17 2022"
if errorlevel 1 exit 1
cmake --build build --config RELEASE
if errorlevel 1 exit 1
cd ..

:: -DCMAKE_CXX_FLAGS="-DARMA_DONT_USE_WRAPPER"

cmake -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX:PATH="%PREFIX%" ^
      -Wno-dev ^
      -S . -B build
if errorlevel 1 exit 1

cmake --build build --target simcoon --config Release
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1

cd simcoon-python-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -Wno-dev ^
      -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 exit 1
cmake --build build --target ALL_BUILD --config Release
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1

cd %SP_DIR%
mkdir simcoon
cd simcoon
type NUL > __init__.py
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\lib\Release\simmit.pyd %SP_DIR%\simcoon
if errorlevel 1 exit 1

xcopy /s /i %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\bin\Release\arma2numpy.dll %PREFIX%\Library\bin
if errorlevel 1 exit 1