Python and machine-learning constitutive laws
-----------------------------------------------

Constitutive laws written in Python and driven by the C++ material-point solver
through the ``PYEXT`` callback (see :doc:`/simulation/python_umat`), and the
recurrent neural network models of :mod:`simcoon.ml` (see :doc:`/simulation/ml_lstm`).

- **pyumat_numpy_j2** - J2 plasticity written in numpy, compared with ``EPICP``
- **plot_lstm_epicp** - a stress LSTM trained on EPICP paths and run in the solver
  under strain and stress control (requires PyTorch: conda-forge ``pytorch`` in a conda
  environment, ``pip install simcoon[ml]`` otherwise)
- **plot_lmsc_selfconsistency** - the same data through a gated LSTM and through the
  linearized minimal state cell, comparing stationarity and self-consistency
