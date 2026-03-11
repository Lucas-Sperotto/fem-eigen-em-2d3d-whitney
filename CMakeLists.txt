cmake_minimum_required(VERSION 3.16)
project(tp3485_fem_eigen_em LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile-time project root used by output path helpers.
add_compile_definitions(TP3485_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(helm10_rect
  src/helm10/main_helm10_rect.cpp
  src/core/mesh2d_rect.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/helm10_scalar_system.cpp
)

add_executable(helm10_circle
  src/helm10/main_helm10_circle.cpp
  src/core/mesh2d_circle.cpp
  src/core/helm10_scalar_system.cpp
)

add_executable(helm10_coax
  src/helm10/main_helm10_coax.cpp
  src/core/mesh2d_coax.cpp
  src/core/helm10_scalar_system.cpp
)

add_executable(edge_rect
  src/edge/edge_dofs.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_assembly.cpp
  src/helmvec/main_edge_rect.cpp
  src/edge/mode_match_rect_edge.hpp
  src/core/mesh2d_rect.cpp
)

add_executable(edge_circle
  src/edge/edge_dofs.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_assembly.cpp
  src/helmvec/main_edge_circle.cpp
  src/edge/mode_match_circle_edge.hpp
  src/core/mesh2d_circle.cpp
)

add_executable(edge_coax
  src/edge/edge_dofs.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_assembly.cpp
  src/helmvec/main_edge_coax.cpp
  src/edge/mode_match_coax_edge.hpp
  src/core/mesh2d_coax.cpp
)

add_executable(mixed_rect
  src/helmvec1/main_mixed_rect.cpp
  src/helmvec1/helmvec1_mixed_system.cpp
  src/edge/edge_assembly.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_dofs.cpp
  src/core/helm10_scalar_system.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/mesh2d_rect.cpp
)

add_executable(mixed_circle
  src/helmvec1/main_mixed_circle.cpp
  src/helmvec1/helmvec1_mixed_system.cpp
  src/edge/edge_assembly.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_dofs.cpp
  src/core/helm10_scalar_system.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/mesh2d_circle.cpp
)

add_executable(mixed_coax
  src/helmvec1/main_mixed_coax.cpp
  src/helmvec1/helmvec1_mixed_system.cpp
  src/edge/edge_assembly.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_dofs.cpp
  src/core/helm10_scalar_system.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/mesh2d_coax.cpp
)

add_executable(helmvec2_rect
  src/helmvec2/main_helmvec2_rect.cpp
  src/helmvec2/helmvec2_coupled_system.cpp
  src/edge/edge_assembly.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_dofs.cpp
  src/core/helm10_scalar_system.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/mesh2d_rect.cpp
)

add_executable(helmvec3_rect
  src/helmvec3/main_helmvec3_rect.cpp
  src/helmvec2/helmvec2_coupled_system.cpp
  src/edge/edge_assembly.cpp
  src/edge/edge_basis.cpp
  src/edge/edge_dofs.cpp
  src/core/helm10_scalar_system.cpp
  src/core/fem_scalar_helm10.cpp
  src/core/mesh2d_rect.cpp
)

add_executable(fem3d0_rect
  src/fem3d0/main_fem3d0_rect.cpp
  src/core/mesh3d_rect.cpp
  src/core/mesh3d_cylinder.cpp
  src/core/mesh3d_sphere.cpp
  src/edge3d/edge3d_basis.cpp
  src/edge3d/edge3d_dofs.cpp
  src/edge3d/edge3d_assembly.cpp
)

add_executable(fem3d1_rect
  src/fem3d1/main_fem3d1_rect.cpp
  src/core/mesh3d_rect.cpp
  src/core/mesh3d_cylinder.cpp
  src/core/mesh3d_sphere.cpp
  src/edge3d/edge3d_basis.cpp
  src/edge3d/edge3d_dofs.cpp
  src/edge3d/edge3d_assembly.cpp
)

target_include_directories(helm10_rect PRIVATE src)
target_include_directories(helm10_circle PRIVATE src)
target_include_directories(helm10_coax PRIVATE src)
target_include_directories(edge_rect PRIVATE src)
target_include_directories(edge_circle PRIVATE src)
target_include_directories(edge_coax PRIVATE src)
target_include_directories(mixed_rect PRIVATE src)
target_include_directories(mixed_circle PRIVATE src)
target_include_directories(mixed_coax PRIVATE src)
target_include_directories(helmvec2_rect PRIVATE src)
target_include_directories(helmvec3_rect PRIVATE src)
target_include_directories(fem3d0_rect PRIVATE src)
target_include_directories(fem3d1_rect PRIVATE src)


# LAPACKE + BLAS/LAPACK
find_library(LAPACKE_LIB lapacke)
find_library(LAPACK_LIB lapack)
find_library(BLAS_LIB openblas)

if (NOT LAPACKE_LIB OR NOT LAPACK_LIB OR NOT BLAS_LIB)
  message(FATAL_ERROR "Nao achei LAPACKE/LAPACK/BLAS. Instale: liblapacke-dev libopenblas-dev")
endif()

target_link_libraries(helm10_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(helm10_circle PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(helm10_coax PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(edge_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(edge_circle PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(edge_coax PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(mixed_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(mixed_circle PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(mixed_coax PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(helmvec2_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(helmvec3_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(fem3d0_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
target_link_libraries(fem3d1_rect PRIVATE ${LAPACKE_LIB} ${LAPACK_LIB} ${BLAS_LIB})
