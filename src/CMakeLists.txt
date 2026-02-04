# Source files for simcoon library
# Listed explicitly to ensure CMake detects changes

target_sources(simcoon PRIVATE
  # Continuum Mechanics - Functions
  Continuum_mechanics/Functions/constitutive.cpp
  Continuum_mechanics/Functions/contimech.cpp
  Continuum_mechanics/Functions/criteria.cpp
  Continuum_mechanics/Functions/damage.cpp
  Continuum_mechanics/Functions/derivatives.cpp
  Continuum_mechanics/Functions/func_N.cpp
  Continuum_mechanics/Functions/hyperelastic.cpp
  Continuum_mechanics/Functions/kinematics.cpp
  Continuum_mechanics/Functions/natural_basis.cpp
  Continuum_mechanics/Functions/objective_rates.cpp
  Continuum_mechanics/Functions/recovery_props.cpp
  Continuum_mechanics/Functions/stress.cpp
  Continuum_mechanics/Functions/transfer.cpp

  # Continuum Mechanics - Homogenization
  Continuum_mechanics/Homogenization/cylinder_multi.cpp
  Continuum_mechanics/Homogenization/ellipsoid_multi.cpp
  Continuum_mechanics/Homogenization/eshelby.cpp
  Continuum_mechanics/Homogenization/layer_multi.cpp
  Continuum_mechanics/Homogenization/phase_multi.cpp

  # Continuum Mechanics - Material
  Continuum_mechanics/Material/crystallo.cpp
  Continuum_mechanics/Material/ODF.cpp
  Continuum_mechanics/Material/ODF2Nphases.cpp
  Continuum_mechanics/Material/PDF.cpp
  Continuum_mechanics/Material/PDF2Nphases.cpp
  Continuum_mechanics/Material/peak.cpp
  Continuum_mechanics/Material/read.cpp
  Continuum_mechanics/Material/variant.cpp

  # Continuum Mechanics - Micromechanics
  Continuum_mechanics/Micromechanics/multiphase.cpp
  Continuum_mechanics/Micromechanics/schemes.cpp

  # Continuum Mechanics - UMAT Finite
  Continuum_mechanics/Umat/Finite/generic_hyper_invariants.cpp
  Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.cpp
  Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.cpp
  Continuum_mechanics/Umat/Finite/mooney_rivlin.cpp
  Continuum_mechanics/Umat/Finite/Neo_hookean_comp.cpp
  Continuum_mechanics/Umat/Finite/Neo_hookean_incomp.cpp
  Continuum_mechanics/Umat/Finite/saint_venant.cpp

  # Continuum Mechanics - UMAT Mechanical Combined
  Continuum_mechanics/Umat/Mechanical/Combined/Prony_Nfast_Plastic.cpp

  # Continuum Mechanics - UMAT Mechanical Damage
  Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.cpp
  Continuum_mechanics/Umat/Mechanical/Damage/damage_weibull.cpp

  # Continuum Mechanics - UMAT Mechanical Elasticity
  Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_isotropic.cpp
  Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.cpp
  Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.cpp

  # Continuum Mechanics - UMAT Mechanical Plasticity
  Continuum_mechanics/Umat/Mechanical/Plasticity/Ani_chaboche_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/DFA_chaboche_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/Generic_chaboche_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_chaboche_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_isoh_Nfast.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_chaboche_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.cpp
  Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_kin_iso_ccp.cpp

  # Continuum Mechanics - UMAT Mechanical SMA
  Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono.cpp
  Continuum_mechanics/Umat/Mechanical/SMA/SMA_mono_cubic.cpp
  Continuum_mechanics/Umat/Mechanical/SMA/unified_T.cpp
  Continuum_mechanics/Umat/Mechanical/SMA/aniso_T.cpp

  # Continuum Mechanics - UMAT Mechanical Viscoelasticity
  Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Prony_Nfast.cpp
  Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_fast.cpp
  Continuum_mechanics/Umat/Mechanical/Viscoelasticity/Zener_Nfast.cpp

  # Continuum Mechanics - UMAT Thermomechanical Elasticity
  Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.cpp
  Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.cpp
  Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_transverse_isotropic.cpp

  # Continuum Mechanics - UMAT Thermomechanical Plasticity
  Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic_ccp.cpp
  Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso_ccp.cpp

  # Continuum Mechanics - UMAT Thermomechanical SMA
  Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.cpp

  # Continuum Mechanics - UMAT Thermomechanical Viscoelasticity
  Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Prony_Nfast.cpp
  Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_fast.cpp
  Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_Nfast.cpp

  # Continuum Mechanics - UMAT Core
  Continuum_mechanics/Umat/fea_transfer.cpp
  Continuum_mechanics/Umat/umat_L_elastic.cpp
  Continuum_mechanics/Umat/umat_smart.cpp

  # Simulation - Geometry
  Simulation/Geometry/cylinder.cpp
  Simulation/Geometry/ellipsoid.cpp
  Simulation/Geometry/geometry.cpp
  Simulation/Geometry/layer.cpp

  # Simulation - Identification
  Simulation/Identification/constants.cpp
  Simulation/Identification/doe.cpp
  Simulation/Identification/generation.cpp
  Simulation/Identification/identification.cpp
  Simulation/Identification/individual.cpp
  Simulation/Identification/methods.cpp
  Simulation/Identification/opti_data.cpp
  Simulation/Identification/optimize.cpp
  Simulation/Identification/parameters.cpp
  Simulation/Identification/read.cpp
  Simulation/Identification/script.cpp

  # Simulation - Maths
  Simulation/Maths/lagrange.cpp
  Simulation/Maths/num_solve.cpp
  Simulation/Maths/random.cpp
  Simulation/Maths/rotation.cpp
  Simulation/Maths/rotation_class.cpp
  Simulation/Maths/solve.cpp
  Simulation/Maths/stats.cpp

  # Simulation - Phase
  Simulation/Phase/material_characteristics.cpp
  Simulation/Phase/output.cpp
  Simulation/Phase/phase_characteristics.cpp
  Simulation/Phase/read.cpp
  Simulation/Phase/state_variables.cpp
  Simulation/Phase/state_variables_M.cpp
  Simulation/Phase/state_variables_T.cpp
  Simulation/Phase/write.cpp

  # Simulation - Solver
  Simulation/Solver/block.cpp
  Simulation/Solver/output.cpp
  Simulation/Solver/read.cpp
  Simulation/Solver/solver.cpp
  Simulation/Solver/step.cpp
  Simulation/Solver/step_meca.cpp
  Simulation/Solver/step_thermomeca.cpp
)
