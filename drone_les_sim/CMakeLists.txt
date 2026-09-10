cmake_minimum_required(VERSION 3.16)
project(drone_les_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wpedantic)

# Recherche du package OpenMP (utilisé par SmagorinskyModel)
find_package(OpenMP REQUIRED)

# Recupération de tous les fichiers .cpp dans src/
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Exécutable
add_executable(sim_drone ${SOURCES})

# Inclusion des répertoires pour résoudre "physics/..."
target_include_directories(sim_drone PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/physics
)

# Liaisons des bibliothèques
target_link_libraries(sim_drone PRIVATE OpenMP::OpenMP_CXX)