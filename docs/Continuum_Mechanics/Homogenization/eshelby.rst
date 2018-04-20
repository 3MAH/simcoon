The Eshelby tensor library
===================

.. default-domain:: cpp

.. function:: mat Eshelby_sphere(double)

    Provides the Eshelby tensor of a spherical inclusion for isotropic linear elasticity in the Simcoon formalism. Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor, as a function of the Poisson ratio :math:`\nu`
    
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

    Provides the Eshelby tensor of a cylindrical inclusion for isotropic linear elasticity in the Simcoon formalism, as a function of the Poisson ratio :math:`\nu`. The cylinder is oriented such as the longitudinal axis is the axis :math:`1`. Returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor.
  
      .. math::

        \boldsymbol{S}=\left(\begin{matrix}
        \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & 0 & 0 & 0 \\
        \frac{5\nu-1}{15(1-\nu)} & \frac{5\nu-1}{15(1-\nu)} & \frac{7-5\nu}{15(1-\nu)} & 0 & 0 & 0 \\
        0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 & 0 \\
        0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} & 0 \\
        0 & 0 & 0 & 0 & 0 & 2\frac{4-5\nu}{15(1-\nu)} \end{matrix}\right)

    .. code-block:: cpp

        mat S = Eshelby_cylinder(nu);
  
.. function:: mat Eshelby_prolate(double,double)


    Provides the Eshelby tensor of a prolate inclusion for isotropic linear elasticity in the Simcoon formalism, as a function of the Poisson ratio :math:`\nu` and the aspect ratio :math:`a_r = frac{a1}{a2} = frac{a1}{a3}`. The prolate inclusion is oriented such as the axis of rotation is the axis :math:`1`.
  
      .. math::

        \boldsymbol{S}=\left(\begin{matrix} S_{11} & S_{12} & S_{12} & 0 & 0 & 0 \\ 
        S_{21} & S_{22} & S_{23} & 0 & 0 & 0 \\
        S_{21} & S_{23} & S_{22} & 0 & 0 & 0 \\
        0 & 0 & 0 & S_{44} & 0 & 0 \\
        0 & 0 & 0 & 0 & S_{44} & 0 \\
        0 & 0 & 0 & 0 & 0 & S_{66} \end{matrix}\right)
        
    with the following components:
    
      .. math::            
        
        S_{11} &= \frac{1}{2(1-\nu)}\left(1-2\nu+\frac{3a_r^2-1}{a_r^2-1}-g\left(1-2\nu+\frac{3a_r^2}{a_r^2-1}\right)\right) \\
        S_{12} &= \frac{-1}{2(1-\nu)}\left(1-2\nu+\frac{1}{a_r^2-1}+g\left(1-2\nu+\frac{3}{a_r^2-1}\right)\right) \\
        S_{21} &= \frac{-a_r^2}{2(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(\frac{3a_r^2}{a_r^2-1}-\left(1-2\nu\right)\right) \\
        S_{22} &= \frac{3a_r^2}{8(1-\nu)}\left(a_r^2-1\right)+\frac{g}{4\left(1-\nu\right)}\left(1-2\nu-\frac{9}{4\left(a_r^2-1\right)}\right) \\
        S_{23} &= \frac{1}{4(1-\nu)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}-g\left(1-2\nu+\frac{3}{4\left(a_r^2-1\right)}\right)\right) \\
        S_{44} &= \frac{2}{4\left(1-\nu\right)}\left(1-2\nu-\frac{a_r^2+1}{a_r^2-1}-\frac{g}{2}\left(1-2\nu-\frac{3a_r^2+1}{a_r^2-1}\right)\right) \\
        S_{66} &= \frac{2}{4\left(1-\nu\right)}\left(\frac{a_r^2}{2\left(a_r^2-1\right)}+g\left(1-2\nu-\frac{3}{4\left(a_r^2-1\right(}\right)\right) 
        
        
     with :math:`g = a_r\frac{a_r\sqrt{a_r^2-1}}{\left(a_r^2-1\right)^{\frac{3}{2}}} - acos(a_r)`

