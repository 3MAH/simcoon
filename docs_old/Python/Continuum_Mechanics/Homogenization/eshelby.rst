The Eshelby tensor library
===================

The Eshelby Library provides various estimations of the Eshelby tensor and the Hill interaction tensor (also called polarisation tensor in some references). In particular, this library offers an analytical expression for special cases, in the framework on linear elasticity. Also, it provides an numerical estimation of the Eshelby tensor in the framework of an anisotropic linear behavior, for a general ellipsoidal inclusion shape.

.. default-domain:: cpp

    .. code-block:: cpp

        #include <simcoon/Continuum_Mechanics/Homogenization/eshelby.hpp>

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
     
    .. code-block:: cpp

        mat S = Eshelby_prolate(nu,a_r);
        
.. function:: mat Eshelby_oblate(double,double)


    Provides the Eshelby tensor of a oblate inclusion for isotropic linear elasticity in the Simcoon formalism, as a function of the Poisson ratio :math:`\nu` and the aspect ratio :math:`a_r = frac{a1}{a2} = frac{a1}{a3}`. The oblate inclusion is oriented such as the axis of rotation is the axis :math:`1`.
  
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
             
     with :math:`g = a_r\frac{-a_r\sqrt{1-a_r^2}}{\left(1-a_r^2\right)^{\frac{3}{2}}} - acos(a_r)`
     
    .. code-block:: cpp

        mat S = Eshelby_oblate(nu,a_r);        
          
.. function:: mat Eshelby(mat, double, double, double, vec, vec, vec, vec, int, int)

    Provides the numerical estimation of the Eshelby tensor of an ellispoid in the general case of anisotropic media, as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction :math:`1`,:math:`2` and :math:`3`, respectively. It also requires the list of integration and their respective weight for the numerical integration, as well as the number of integration points in the :math:`1` and :math:`2` directions. The points and weights are calculated using the points_  function.

    .. code-block:: cpp
        
        mat S = Eshelby(L, a1, a2, a3, x, wx, y, wy, mp, np);
        
    *L* is the stiffness tensor of the media; *a1* is the semi-axis of the ellispoid length in the direction :math:`1`; *a2* is the semi-axis of the ellispoid length in the direction :math:`2`; *a3* is the semi-axis of the ellipsoid length in the direction :math:`3`; *x* is the vector of points in the direction :math:`1`; *wx* is the vector of the weights of points in the direction :math:`1`; *y* is the vector of points in the direction :math:`2`; *wx* is the vector of the weights of points in the direction :math:`2`; *mp* is the number of points in the direction :math:`1`; *np* is the number of points in the direction :math:`2`;

    The function returns the Eshelby tensor as a mat, according to the conventions of a localisation tensor

.. function:: mat T_II(mat, double, double, double, vec, vec, vec, vec, int, int)

    Provides the numerical estimation of the Hill interaction tensor of an ellispoid in the general case of anisotropic media, as a function of the stiffness tensor, and the three semi-axis length of the ellipsoid in the direction :math:`1`,:math:`2` and :math:`3`, respectively. It also requires the list of integration and their respective weight for the numerical integration, as well as the number of integration points in the :math:`1` and :math:`2` directions. The points and weights are calculated using the points_  function.

    .. code-block:: cpp
        
        mat S = T_II(L, a1, a2, a3, x, wx, y, wy, mp, np)
                
    *L* is the stiffness tensor of the media; *a1* is the semi-axis of the ellispoid length in the direction :math:`1`; *a2* is the semi-axis of the ellispoid length in the direction :math:`2`; *a3* is the semi-axis of the ellipsoid length in the direction :math:`3`; *x* is the vector of points in the direction :math:`1`; *wx* is the vector of the weights of points in the direction :math:`1`; *y* is the vector of points in the direction :math:`2`; *wx* is the vector of the weights of points in the direction :math:`2`; *mp* is the number of points in the direction :math:`1`; *np* is the number of points in the direction :math:`2`;

    The function returns the Hill interaction tensor as a mat, according to the conventions of a localisation tensor

.. function:: void points(mat, double, double, double, vec, vec, vec, vec, int, int)

    This methods computes the list of integration and their respective weight for the numerical integration, as a function of the number of integration points in the 1 and 2 directions. 
  
    .. code-block:: cpp
      
        vec x(mp);
        vec wx(mp);
        vec y(np);
        vec wy(np);
        points(x, wx, y, wy, mp, np);
       
    *x* is the vector of points in the direction :math:`1`; *wx* is the vector of the weights of points in the direction :math:`1`; *y* is the vector of points in the direction :math:`2`; *wx* is the vector of the weights of points in the direction :math:`2`; *mp* is the number of points in the direction :math:`1`; *np* is the number of points in the direction :math:`2`.
    Update *x*, *wx*, *y* and *wy* according to *mp* and *np*. Note that *x*, *wx*, *y*, *wy* have to be initialized first with the size of *mp* and *np*, respectively.
    
