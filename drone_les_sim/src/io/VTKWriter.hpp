// io/VTKWriter.hpp
#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "core/Mesh.hpp"
#include "core/Field2D.hpp"
#include <string>
namespace drone {

class VTKWriter {
public:
    // Déclaration exacte de la méthode write
    void write(const std::string& filename, const Mesh& mesh, 
               const Field2D& u, const Field2D& v, const Field2D& p);
};

}

#endif