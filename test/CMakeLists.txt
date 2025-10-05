# Test files for simcoon
# Listed explicitly to ensure CMake detects changes

# Test source files
set(TEST_SRCS
  # Identification tests
  Identification/Tidentification.cpp

  # Libraries - Continuum Mechanics
  Libraries/Continuum_mechanics/Tconstitutive.cpp
  Libraries/Continuum_mechanics/Tcontimech.cpp
  Libraries/Continuum_mechanics/Tcriteria.cpp
  Libraries/Continuum_mechanics/Thyperelastic.cpp
  Libraries/Continuum_mechanics/Tobjective_rates.cpp
  Libraries/Continuum_mechanics/Trecovery_props.cpp
  Libraries/Continuum_mechanics/Tstress.cpp
  Libraries/Continuum_mechanics/Ttransfer.cpp

  # Libraries - Homogenization
  Libraries/Homogenization/Teshelby.cpp

  # Libraries - Maths
  Libraries/Maths/Trotation.cpp

  # Libraries - Phase
  Libraries/Phase/Tphase_characteristics.cpp

  # Libraries - Umat
  Libraries/Umat/TAba2sim.cpp

  # UMATs
  Umats/ELISO/TELISO.cpp
  Umats/ELIST/TELIST.cpp
  Umats/ELORT/TELORT.cpp
  Umats/EPCHA/TEPCHA.cpp
  Umats/EPHAC/TEPHAC.cpp
  Umats/EPICP/TEPICP.cpp
  Umats/EPKCP/TEPKCP.cpp
  Umats/HYPER/THYPER.cpp
  Umats/LOG_int/TLOG_int.cpp
  Umats/MIMTN/TMIMTN.cpp
  Umats/MIPLN/TMIPLN.cpp
  Umats/UMEXT/TUMEXT.cpp
)

# External test files (non-MSVC only)
set(TEST_EXTERN
  ../test_extern/Umats/UMABA/TUMABA.cpp
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
