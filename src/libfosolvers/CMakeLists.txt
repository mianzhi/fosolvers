SET(libfosolvers_src
  
)

SET(libfosolvers_lib
  
)

ADD_LIBRARY(fosolvers STATIC ${libfosolvers_src})
TARGET_LINK_LIBRARIES(fosolvers ${libfosolvers_lib})
