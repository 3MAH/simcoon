Install simcoon on Conda environnement
===========================

Installation of the environement
--------------------------------

#Installation of an environment called "scientific" (python 3.7 ou 3.8)
conda create --name scientific
#activate the environment
. /path/to/anaconda3/bin/activate && conda activate /path/to/anaconda3/envs/scientific; 
#Install packages

conda install -c conda-forge armadillo
conda install -c conda-forge numpy


Installation of simcoon (easy)
-----------------------

Open Install.sh, fill the anaconda path at the beginning of the file (full path) and run:
sh Install.sh

Then just follow the instructions

Installation of simcoon
-----------------------
cd build (or mkdir build and then cd build if not already there)

cmake .. -DCMAKE_INCLUDE_PATH=/path/to/anaconda3/envs/scientific/include -DCMAKE_LIBRARY_PATH=/path/to/anaconda3/envs/scientific/lib
make
make test

To finally install it, it is best to copy the simcoon file that is in the /include folder in the /path/to/anaconda3/envs/scientific/include, and the simcoon library (libsimcoon.so) in the /path/to/anaconda3/envs/scientific/lib folder manually. This last library should be situated in the /build/lib folder (or /build/lib/Relase, or /build/lib/Debug)
You can also specify to cmake these installation folders and run 'make install' command


Installation of simcoon python builder
-----------------------
Go to the simcoon-python-builder and cd build (or mkdir build and then cd build if not already there)

cmake .. -DCMAKE_INCLUDE_PATH=/path/to/anaconda3/envs/scientific/include -DCMAKE_LIBRARY_PATH=/path/to/anaconda3/envs/scientific/lib
make
make test

To finally it, it is best to copy the additional simcoon include folders (arma2numpy and python_wrappers) files that are in the simcoon-python-builder/include folder in the /path/to/anaconda3/envs/scientific/include, and the arma2numpy library (libarma2numpy.so) in the /path/to/anaconda3/envs/scientific/lib folder manually.
The last library should be situated in the simcoon-python-builder/build/lib folder (or /build/lib/Relase, or /build/lib/Debug)

Finally you should copy the '_core.so' file situated in the simcoon-python-builder/build/lib folder (or /build/lib/Relase, or /build/lib/Debug) in the folder 'python-setup/simcoon'

Go to he folder 'python-setup/simcoon' and finally run:
'python setup.py install'


Installation of the Jupyter notebook
--------------------------------


conda install ipython
conda install ipykernel
python -m ipykernel install --user --name=scientific


Possible issues
===============

If you get an error such as:

Error processing line 1 of /Users/ychemisky/opt/anaconda3/envs/test/lib/python3.8/site-packages/matplotlib-3.1.1-py3.8-nspkg.pth:

Just suppress the file 
/Users/ychemisky/opt/anaconda3/envs/test/lib/python3.8/site-packages/matplotlib-3.1.1-py3.8-nspkg.pth:

