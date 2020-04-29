from setuptools import setup
#from distutils.core import setup

setup(name='simcoon',
      version='1.0',
      description='Simulation in Mechanics and Materials: Interactive Tools',
      author='Yves Chemisky',
      author_email='yves.chemisky@gmail.com',
      #url=
      packages=['simcoon',],
      package_data={'simcoon': ['simmit.so']},
      
      include_package_data=True,
      license='GPL'
      )
