The Micromechanics schemes
===================

The Eshelby Library provides various estimations of the Eshelby tensor and the Hill interaction tensor (also called polarisation tensor in some references). In particular, this library offers an analytical expression for special cases, in the framework on linear elasticity. Also, it provides an numerical estimation of the Eshelby tensor in the framework of an anisotropic linear behavior, for a general ellipsoidal inclusion shape.

.. default-domain:: cpp

    .. code-block:: cpp

        #include <simcoon/Continuum_Mechanics/Homogenization/schemes.hpp>

.. function:: void umat_multi(phase_characteristics &, const mat &, const double &, const double &, const int &, const int &, const bool &, double &, const int &)

The procedure umat_multi takes care of the constitutive response of a composite material that possesses :math:`N` distinct phases. 
In this procedure, the *phase_characteristics* object is being updated, with the decomposition of the total strain :math:`\Delta \mathbf{\varepsilon}` and temperature :math`\Delta T` increments.

The Mori Tanaka scheme
----------------------------------

This Library provides the macroscopic response of a composite with N phases, using the Mori Tanaka method. The algorithm  requires the following information: 


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





# The Micromechanics libraries

<div id="Mori_Tanaka_Library">
  <h2>
    The Mori Tanaka Library (Mori_Tanaka.hpp)
  </h2> This Library provides the macroscopic response of a composite with N phases, using the Mori Tanaka method. The algorithm umat_MT_N requires the following information: 
  
  <ul>
    <li>
      Etot (vec) : The total macroscopic strain at the beginning of the increment.
    </li>
    <li>
      DEtot (vec) : The increment of the total macroscopic strain.
    </li>
    <li>
      sigma (vec) : The macroscopic stress (initially at the beginning of the increment, updated at the end).
    </li>
    <li>
      Lt (mat) : The macroscopic tangent stiffness tensor.
    </li>
    <li>
      DR (mat) : The rotation increment matrix.
    </li>
    <li>
      nprops (int) : The number of constants associated with the composite and each phase.
    </li>
    <li>
      *props (double) : A table of material properties: props[0] defines the number of phases, props[1] is the value (X) giving the file number containing phases properties to homogenize (this file is called "NphasesX.dat"), while props[2] and props[3] are the number of integration points in the two directions for the computation of the Eshelby tensors. The rest of the material properties are associated with each phase.
    </li>
    <li>
      nstatev (int) : The number of state variables stored for all the phases.
    </li>
    <li>
      *statev (double) : A table of state variables. At each material phase: the first 6 variables store the total strain, the next 6 the increment of the total strain, the next 6 the stress, the next 36 the elastic stiffness tensor and the next 36 the tangent stiffness tensor of the phase (all these in the global coordinate system). The rest of the statev are related with the constitutive law of the phase (plastic strains, viscous strains etc).
    </li>
    <li>
      T (double) : The macroscopic temperature at the beginning of the increment.
    </li>
    <li>
      DT (double) : The increment of the macroscopic temperature.
    </li>
    <li>
      Time (double): The time at the beginning of the increment.
    </li>
    <li>
      DTime (double): The increment of time.
    </li>
    <li>
      sse (double): The specific elastic strain energy of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      spd (double): The specific plastic dissipation of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      ndi (int): Number of direct stress components used in the analysis.
    </li>
    <li>
      nshr (int): Number of engineering shear stress components used in the analysis.
    </li>
    <li>
      start (bool): It is related with the initialization of the algorithm.
    </li>
  </ul> The algorithm reads the material properties of all the phases from the file "Nphases.dat", which is included in the folder "data". At the end of the computations, the umat_MT_N returns the updated values of the macroscopic stress, the macroscopic tangent stiffness tensor and the statev of each phase.
</div>

* * *

<div id="Self_Consistent_Library">
  <h2>
    The Self Consistent Library (Self_Consistent.hpp)
  </h2> This Library provides the macroscopic response of a composite with N phases, using the self consistent method. The algorithm umat_SC_N requires the following information: 
  
  <ul>
    <li>
      Etot (vec) : The total macroscopic strain at the beginning of the increment.
    </li>
    <li>
      DEtot (vec) : The increment of the total macroscopic strain.
    </li>
    <li>
      sigma (vec) : The macroscopic stress (initially at the beginning of the increment, updated at the end).
    </li>
    <li>
      Lt (mat) : The macroscopic tangent stiffness tensor.
    </li>
    <li>
      DR (mat) : The rotation increment matrix.
    </li>
    <li>
      nprops (int) : The number of constants associated with the composite and each phase.
    </li>
    <li>
      *props (double) : A table of material properties: props[0] defines the number of phases, props[1] is the value (X) giving the file number containing phases properties to homogenize (this file is called "NphasesX.dat"), while props[2] and props[3] are the number of integration points in the two directions for the computation of the Eshelby tensors. The rest of the material properties are associated with each phase.
    </li>
    <li>
      nstatev (int) : The number of state variables stored for all the phases.
    </li>
    <li>
      *statev (double) : A table of state variables. At each material phase: the first 6 variables store the total strain, the next 6 the increment of the total strain, the next 6 the stress, the next 36 the elastic stiffness tensor and the next 36 the tangent stiffness tensor of the phase (all these in the global coordinate system). The rest of the statev are related with the constitutive law of the phase (plastic strains, viscous strains etc).
    </li>
    <li>
      T (double) : The macroscopic temperature at the beginning of the increment.
    </li>
    <li>
      DT (double) : The increment of the macroscopic temperature.
    </li>
    <li>
      Time (double): The time at the beginning of the increment.
    </li>
    <li>
      DTime (double): The increment of time.
    </li>
    <li>
      sse (double): The specific elastic strain energy of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      spd (double): The specific plastic dissipation of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      ndi (int): Number of direct stress components used in the analysis.
    </li>
    <li>
      nshr (int): Number of engineering shear stress components used in the analysis.
    </li>
    <li>
      start (bool): It is related with the initialization of the algorithm.
    </li>
  </ul> The algorithm reads the material properties of all the phases from the file "Nphases.dat", which is included in the folder "data". At the end of the computations, the umat_SC_N returns the updated values of the macroscopic stress, the macroscopic tangent stiffness tensor and the statev of each phase.
</div>

* * *

<div id="Periodic_Layer_Library">
  <h2>
    The Periodic Layers Library (Periodic_Layer.hpp)
  </h2> This Library provides the macroscopic response of a multilayered composite with N layers, using the periodic homogenization method. The algorithm umat_PL_N requires the following information: 
  
  <ul>
    <li>
      Etot (vec) : The total macroscopic strain at the beginning of the increment.
    </li>
    <li>
      DEtot (vec) : The increment of the total macroscopic strain.
    </li>
    <li>
      sigma (vec) : The macroscopic stress (initially at the beginning of the increment, updated at the end).
    </li>
    <li>
      Lt (mat) : The macroscopic tangent stiffness tensor.
    </li>
    <li>
      DR (mat) : The rotation increment matrix.
    </li>
    <li>
      nprops (int) : The number of constants associated with the composite and each phase.
    </li>
    <li>
      *props (double) : A table of material properties: props[0] defines the number of phases, while the rest of the material properties are associated with each phase.
    </li>
    <li>
      nstatev (int) : The number of state variables stored for all the phases.
    </li>
    <li>
      *statev (double) : A table of state variables. At each material phase: the first 6 variables store the total strain, the next 6 the increment of the total strain, the next 6 the stress, the next 36 the elastic stiffness tensor and the next 36 the tangent stiffness tensor of the phase (all these in the global coordinate system). The rest of the statev are related with the constitutive law of the phase (plastic strains, viscous strains etc).
    </li>
    <li>
      T (double) : The macroscopic temperature at the beginning of the increment.
    </li>
    <li>
      DT (double) : The increment of the macroscopic temperature.
    </li>
    <li>
      Time (double): The time at the beginning of the increment.
    </li>
    <li>
      DTime (double): The increment of time.
    </li>
    <li>
      sse (double): The specific elastic strain energy of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      spd (double): The specific plastic dissipation of the composite. Given at the beginning of the increment, updated at the end (unused in this version of SMART+).
    </li>
    <li>
      ndi (int): Number of direct stress components used in the analysis.
    </li>
    <li>
      nshr (int): Number of engineering shear stress components used in the analysis.
    </li>
    <li>
      start (bool): It is related with the initialization of the algorithm.
    </li>
  </ul> The algorithm reads the material properties of all the phases from the file "Nlayers.dat", which is included in the folder "data". At the end of the computations, the umat_PL_N returns the updated values of the macroscopic stress, the macroscopic tangent stiffness tensor and the statev of each phase.
</div>


