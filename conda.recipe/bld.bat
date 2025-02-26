:: Simcoon
conda install set3mah::carma

cd %SRC_DIR%
cmake -G"Visual Studio 17 2022" ^
      -S . -B build ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_EXECUTABLE=%PREFIX%/python.exe ^
      -Wno-dev

cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cmake -G"Visual Studio 17 2022" ^
      -S simcoon-python-builder -B simcoon-python-builder/build ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DUSE_OPENMP=False ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -DPython3_EXECUTABLE=%PREFIX%/python.exe ^
      -Wno-dev

cmake --build simcoon-python-builder/build --target simmit --config Release
cmake --install simcoon-python-builder/build
if errorlevel 1 exit 1

copy %SRC_DIR%\simcoon-python-builder\build\lib\simmit.pyd %SRC_DIR%\python-setup\simcoon\
python -m pip install %SRC_DIR%\python-setup
if errorlevel 1 exit 1
