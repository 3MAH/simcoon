# Test files for simcoon
# Listed explicitly to ensure CMake detects changes

# Test source files (relative to CMAKE_SOURCE_DIR)
set(TEST_SRCS
  # Libraries - Continuum Mechanics
  test/Libraries/Continuum_mechanics/Tconstitutive.cpp
  test/Libraries/Continuum_mechanics/Tcontimech.cpp
  test/Libraries/Continuum_mechanics/Tcriteria.cpp
  test/Libraries/Continuum_mechanics/Tdamage.cpp
  test/Libraries/Continuum_mechanics/Tderivatives.cpp
  test/Libraries/Continuum_mechanics/Thyperelastic.cpp
  test/Libraries/Continuum_mechanics/Tkinematics.cpp
  test/Libraries/Continuum_mechanics/Tnatural_basis.cpp
  test/Libraries/Continuum_mechanics/Tobjective_rates.cpp
  test/Libraries/Continuum_mechanics/Trecovery_props.cpp
  test/Libraries/Continuum_mechanics/Tstress.cpp
  test/Libraries/Continuum_mechanics/Ttensor.cpp
  test/Libraries/Continuum_mechanics/Ttransfer.cpp

  # Libraries - Homogenization
  test/Libraries/Homogenization/Teshelby.cpp

  # Libraries - Maths
  test/Libraries/Maths/Tlagrange.cpp
  test/Libraries/Maths/Tnum_solve.cpp
  test/Libraries/Maths/Trotation.cpp
  test/Libraries/Maths/Trotation_class.cpp
  test/Libraries/Maths/Tstats.cpp

  # Libraries - Phase
  test/Libraries/Phase/Tphase_characteristics.cpp

  # Libraries - Umat
  test/Libraries/Umat/TAba2sim.cpp
  test/Libraries/Umat/Ttangent_assembly.cpp
  test/Libraries/Umat/Ttangent_EPICP.cpp
  test/Libraries/Umat/Ttangent_EPHIL.cpp
  test/Libraries/Umat/Ttangent_EPCHA.cpp
  test/Libraries/Umat/Tunified_TR.cpp
  test/Libraries/Umat/Modular/Tmodular_umat.cpp
  test/Libraries/Umat/Treference_umats.cpp
  test/Libraries/Umat/Ttutorial_umat.cpp

  # UMATs
  test/Umats/ELISO/TELISO.cpp
  test/Umats/ELIST/TELIST.cpp
  test/Umats/ELORT/TELORT.cpp
  test/Umats/EPCHA/TEPCHA.cpp
  test/Umats/EPCHG/TEPCHG.cpp
  test/Umats/EPHAC/TEPHAC.cpp
  test/Umats/EPICP/TEPICP.cpp
  test/Umats/EPKCP/TEPKCP.cpp
  test/Umats/HYPER/THYPER.cpp
  test/Umats/LOG_int/TLOG_int.cpp
  test/Umats/MIMTN/TMIMTN.cpp
  test/Umats/MIMTP/TMIMTP.cpp
  test/Umats/MIPLN/TMIPLN.cpp
  test/Umats/UMEXT/TUMEXT.cpp
)

# External test files (non-MSVC only)
set(TEST_EXTERN
  test_extern/Umats/UMABA/TUMABA.cpp
)

# Add test executables
foreach(testSrc ${TEST_SRCS})
  add_simcoon_test(${testSrc})
endforeach()

# Reference oracle kernels: the legacy UMAT implementations removed in the
# modular rationalization, kept verbatim under reference_kernels/ and compiled
# ONLY into Treference_umats (never into libsimcoon) to validate the modular
# adapters against independent physics.
file(GLOB REFERENCE_KERNEL_SRCS ${CMAKE_SOURCE_DIR}/test/Libraries/Umat/reference_kernels/*.cpp)
target_sources(Treference_umats PRIVATE ${REFERENCE_KERNEL_SRCS})
target_include_directories(Treference_umats PRIVATE ${CMAKE_SOURCE_DIR}/test/Libraries/Umat)

# Tutorial UMAT (examples/umat_tutorial): compiled into its validation test so
# the hands-on tutorial code is continuously checked against EPICP, the shared
# tangent assembly, and a finite-difference Jacobian.
target_sources(Ttutorial_umat PRIVATE ${CMAKE_SOURCE_DIR}/examples/umat_tutorial/umat_tutorial_J2.cpp)
target_include_directories(Ttutorial_umat PRIVATE ${CMAKE_SOURCE_DIR}/examples/umat_tutorial)

# Add external test executables (non-MSVC only)
if(NOT MSVC)
  foreach(testSrc ${TEST_EXTERN})
    add_simcoon_test(${testSrc})
  endforeach()
endif()
