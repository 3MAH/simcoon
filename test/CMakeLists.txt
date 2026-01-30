# Test files for simcoon v2.0
# C++ unit tests for core library functions
# UMAT solver tests use pytest: pytest python-setup/test/

set(TEST_SRCS
  # Continuum Mechanics
  test/Libraries/Continuum_mechanics/Tconstitutive.cpp
  test/Libraries/Continuum_mechanics/Tcontimech.cpp
  test/Libraries/Continuum_mechanics/Tcriteria.cpp
  test/Libraries/Continuum_mechanics/Thyperelastic.cpp
  test/Libraries/Continuum_mechanics/Tobjective_rates.cpp
  test/Libraries/Continuum_mechanics/Trecovery_props.cpp
  test/Libraries/Continuum_mechanics/Tstress.cpp
  test/Libraries/Continuum_mechanics/Ttransfer.cpp

  # Homogenization
  test/Libraries/Homogenization/Teshelby.cpp

  # Maths
  test/Libraries/Maths/Trotation.cpp
)

foreach(testSrc ${TEST_SRCS})
  add_simcoon_test(${testSrc})
endforeach()
