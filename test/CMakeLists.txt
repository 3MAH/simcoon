# Test files for simcoon
# Listed explicitly to ensure CMake detects changes

# Test source files (relative to CMAKE_SOURCE_DIR)
set(TEST_SRCS
  # Identification tests
  test/Identification/Tidentification.cpp

  # Libraries - Continuum Mechanics
  test/Libraries/Continuum_mechanics/Tconstitutive.cpp
  test/Libraries/Continuum_mechanics/Tcontimech.cpp
  test/Libraries/Continuum_mechanics/Tcriteria.cpp
  test/Libraries/Continuum_mechanics/Thyperelastic.cpp
  test/Libraries/Continuum_mechanics/Tobjective_rates.cpp
  test/Libraries/Continuum_mechanics/Trecovery_props.cpp
  test/Libraries/Continuum_mechanics/Tstress.cpp
  test/Libraries/Continuum_mechanics/Ttransfer.cpp

  # Libraries - Homogenization
  test/Libraries/Homogenization/Teshelby.cpp

  # Libraries - Maths
  test/Libraries/Maths/Trotation.cpp

  # Libraries - Phase
  test/Libraries/Phase/Tphase_characteristics.cpp

  # Libraries - Umat
  test/Libraries/Umat/TAba2sim.cpp
  test/Libraries/Umat/Modular/Tmodular_umat.cpp

  # UMATs
  test/Umats/ELISO/TELISO.cpp
  test/Umats/ELIST/TELIST.cpp
  test/Umats/ELORT/TELORT.cpp
  test/Umats/EPCHA/TEPCHA.cpp
  test/Umats/EPHAC/TEPHAC.cpp
  test/Umats/EPICP/TEPICP.cpp
  test/Umats/EPKCP/TEPKCP.cpp
  test/Umats/HYPER/THYPER.cpp
  test/Umats/LOG_int/TLOG_int.cpp
  test/Umats/MIMTN/TMIMTN.cpp
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

# Add external test executables (non-MSVC only)
if(NOT MSVC)
  foreach(testSrc ${TEST_EXTERN})
    add_simcoon_test(${testSrc})
  endforeach()
endif()
