The Eshelby tensor library
===================

.. default-domain:: cpp

.. function:: mat Eshelby_sphere(double)

    Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity in the Simcoon formalism. Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`nu`
    
    .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{matrix}\right)

    .. code-block:: cpp

        mat S = Eshelby_sphere(nu);

.. function:: mat Eshelby_cylinder(double)

    Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity in the SMART+ formalism, as a function of the Poisson ratio \(\nu\)The cylinder is oriented such as the longitudinal axis is the axis :math:`1`. Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor
  
      .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{matrix}\right)

  
  <pre>mat S = Eshelby_cylinder(nu);</pre> return the Eshelby tensor as a mat, according to the conventions of a localisation tensor \[\boldsymbol{S}=\left(\begin{matrix} 0 & 0 & 0 & 0 & 0 & 0 \\ \frac{\nu}{2(1-\nu)} & \frac{5-4\nu}{8(1-\nu)} & \frac{4\nu-1}{8(1-\nu)} & 0 & 0 & 0 \\ \frac{\nu}{2(1-\nu)} & \frac{4\nu-1}{8(1-\nu)} & \frac{5-4\nu}{8(1-\nu)} & 0 & 0 & 0 \\ 0 & 0 & 0 & 2\frac{1}{4} & 0 & 0 \\ 0 & 0 & 0 & 0 & 2\frac{1}{4} & 0 \\ 0 & 0 & 0 & 0 & 0 & 2\frac{2(3-4\nu)}{8(1-\nu)} \end{matrix}\right)\]
</div>