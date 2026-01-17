==============
Homogenization
==============

This module provides mean-field homogenization methods for predicting 
effective properties of composite materials.

Eshelby Tensor
==============

The foundation of mean-field homogenization is the Eshelby tensor, which 
describes the strain field inside an ellipsoidal inclusion embedded in an 
infinite matrix.

.. doxygengroup:: eshelby
   :project: simcoon
   :content-only:

Homogenization Schemes
======================

Functions for computing effective properties using Mori-Tanaka, 
self-consistent, and other mean-field schemes.

.. doxygengroup:: homogenization
   :project: simcoon
   :content-only:
