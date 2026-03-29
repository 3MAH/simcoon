
Installation
============

Quick install
-------------

**Conda (recommended)**

.. code-block:: bash

    conda install -c conda-forge -c set3mah simcoon

**pip**

.. code-block:: bash

    pip install simcoon

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (arm64),
and Windows (x64) on Python 3.10--3.14.


BLAS/LAPACK and OpenMP
----------------------

Simcoon uses `Armadillo <http://arma.sourceforge.net>`_ for linear algebra,
which in turn relies on a BLAS/LAPACK implementation. The threading model of
the BLAS library matters because it can conflict with OpenMP if both are loaded
in the same process.

OpenMP is enabled on **Linux and macOS** for parallel batch operations.
On **Windows**, OpenMP is disabled in the Python bindings to avoid conflicts
between different OpenMP runtimes (e.g. ``vcomp140.dll`` from MSVC and
runtimes from other packages). Batch operations on Windows run sequentially
while BLAS handles internal threading.

.. list-table:: BLAS/OpenMP summary
   :header-rows: 1
   :widths: 20 25 25 30

   * - Platform
     - BLAS
     - OpenMP
     - Conflict risk
   * - macOS
     - Apple Accelerate
     - ON (libomp, bundled in wheel)
     - None
   * - Linux
     - System OpenBLAS
     - ON (libgomp)
     - None
   * - Windows
     - vcpkg OpenBLAS (pip) / netlib+MKL (conda)
     - OFF
     - None

**Using MKL with conda on Linux**

If you prefer Intel MKL for performance, switch the BLAS backend and set the
threading layer to avoid conflicts between ``libiomp5`` (Intel) and
``libgomp`` (GCC):

.. code-block:: bash

    conda install libblas=*=*mkl mkl
    export MKL_THREADING_LAYER=GNU


Developer installation
----------------------

Prerequisites (conda)
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    conda create --name simcoon_dev
    conda activate simcoon_dev

**Linux:**

.. code-block:: bash

    conda env update -f environment.yml

**macOS (Apple Silicon):**

.. code-block:: bash

    conda env update -f environment_arm64.yml

**Windows:**

.. code-block:: bash

    conda env update -f environment_win.yml

Prerequisites (system packages)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Debian/Ubuntu:**

  .. code-block:: bash

      sudo apt-get install libarmadillo-dev libopenblas-dev liblapack-dev \
          libgtest-dev ninja-build cmake

- **macOS (Homebrew):**

  .. code-block:: bash

      brew install armadillo ninja cmake

- **Windows (vcpkg):**

  .. code-block:: powershell

      vcpkg install armadillo:x64-windows openblas:x64-windows


Building from source
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/3MAH/simcoon.git
    cd simcoon
    pip install -e . --no-build-isolation

This builds the C++ library and Python bindings in one step using
scikit-build-core.

**Enabling OpenMP (optional, for conda environments):**

.. code-block:: bash

    pip install -e . --no-build-isolation \
        --config-settings=cmake.define.SIMCOON_USE_OPENMP=ON


Running tests
~~~~~~~~~~~~~

**Python tests:**

.. code-block:: bash

    pytest

**C++ tests:**

.. code-block:: bash

    mkdir build && cd build
    cmake .. -DSIMCOON_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
    cmake --build .
    ctest --output-on-failure


Using simcoon with fedoo
~~~~~~~~~~~~~~~~~~~~~~~~~

Simcoon is designed to work with `fedoo <https://github.com/3MAH/fedoo>`_
for finite-element simulations. Both packages can be installed together:

.. code-block:: bash

    # conda
    conda install -c conda-forge -c set3mah simcoon fedoo

    # pip
    pip install simcoon fedoo

No special configuration is needed -- the BLAS/OpenMP setup described above
ensures that both libraries coexist without runtime conflicts.
