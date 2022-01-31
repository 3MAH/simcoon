The Constitutive Library
========================

.. default-domain:: python

.. function:: mat Ireal()

    Provides the fourth order identity tensor written in Voigt notation :math:`I_{real}`, where :

    .. math::

        I_{real} = \left( \begin{array}{cccccc}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0.5 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0.5 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0.5 \end{array} \right)

    .. code-block:: cpp

        mat Ir = Ireal();

