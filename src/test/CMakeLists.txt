SET(tests_src tests.f90 showTestRst.f90)
ADD_EXECUTABLE(tests ${tests_src})
TARGET_LINK_LIBRARIES(tests fosolvers)

SET(test_src test.f90 showTestRst.f90)
ADD_EXECUTABLE(test ${test_src})
TARGET_LINK_LIBRARIES(test fosolvers)
