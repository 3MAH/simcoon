:: Simcoon
cd %SRC_DIR%
cmake -S . -B build ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^      
      -Wno-dev

cmake --build build --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cmake -S simcoon-python-builder -B simcoon-python-builder/build ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
      -Wno-dev

cmake --build simcoon-python-builder/build --config Release
cmake --install simcoon-python-builder/build
if errorlevel 1 exit 1

copy %SRC_DIR%/simcoon-python-builder/build/lib/simmit.pyd %SRC_DIR%/python-setup/simcoon/simmit.pyd
python -m pip install %SRC_DIR%/python-setup
if errorlevel 1 exit 1
