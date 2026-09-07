
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
     - Duplicate libomp when mixed with conda packages (see below)
   * - Linux
     - System OpenBLAS
     - ON (libgomp)
     - None
   * - Windows
     - vcpkg OpenBLAS (pip) / netlib+MKL (conda)
     - OFF
     - None

Duplicate OpenMP runtimes on macOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**The rule: one environment = one OpenMP runtime.** In a conda environment that
runtime is conda-forge's ``llvm-openmp`` (``libomp.dylib``), and every native
package must share it. Two setups break the rule:

- a **PyPI wheel** that bundles its own ``libomp`` (PyTorch wheels do; the
  simcoon wheel does, via delocate) next to conda numpy/scipy/simcoon;
- a **source build** of simcoon linked against an Armadillo that does not come
  from the environment (a system or Homebrew copy): its OpenBLAS drags a second
  ``libomp`` into the process.

Symptoms range from the explicit abort message (``OMP: Error #15: Initializing
libomp.dylib, but found libomp.dylib already initialized``) to silent crashes
(``SIGSEGV`` inside ``libomp`` worker threads under threaded runs).

What to do:

1. **conda-forge for everything native**: simcoon, numpy, scipy, armadillo and
   PyTorch (``conda install -c conda-forge pytorch``, not ``pip install torch``).
   Do not mix in Homebrew or system libraries.
2. **Source builds in a conda environment** pick the environment's Armadillo
   automatically (CMake detects ``CONDA_PREFIX``) and warn when they do not.
   Install it first and purge any stale build cache when switching:

   .. code-block:: bash

       conda install -c conda-forge armadillo
       rm -rf build/           # a cached non-conda Armadillo path would be kept otherwise
       pip install -e . --no-build-isolation

3. **Check** which runtimes a process really loads:

   .. code-block:: bash

       python -m simcoon.doctor

   It imports numpy, scipy, simcoon and torch (if present), lists the OpenMP
   runtimes each one brings and exits with 1 when there are several.

``KMP_DUPLICATE_LIB_OK=TRUE`` is **not** a remedy: it silences the guard and
leaves two runtimes fighting (crashes or silently wrong results remain
possible, this is Intel's own warning).

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

- **macOS:** use the conda environment (``environment_arm64.yml`` above:
  conda-forge Armadillo, Accelerate BLAS, cmake, ninja). Do not install the
  dependencies with Homebrew: a Homebrew Armadillo links Homebrew's OpenBLAS
  and ``libomp``, i.e. a second OpenMP runtime next to the environment's (see
  *Duplicate OpenMP runtimes on macOS* above).

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

Keep both packages (and numpy/scipy) on the **same channel** — all-conda or
all-PyPI — and the BLAS/OpenMP setup described above ensures they coexist
without runtime conflicts (on macOS, mixing channels can load a second
OpenMP runtime; see *Duplicate OpenMP runtimes on macOS*).
