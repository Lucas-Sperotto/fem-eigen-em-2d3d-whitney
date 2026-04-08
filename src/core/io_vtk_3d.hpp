/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/io_vtk_3d.hpp                                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios VTK para malhas tetraedricas 3D.                    */
/*****************************************************************************/

#pragma once

#include "mesh3d.hpp"

#include <fstream>
#include <string>
#include <vector>

/******************************************************************************/
/* FUNCAO: write_vtk_unstructured_tet_cell_vector                             */
/* DESCRICAO: Escreve malha tetraedrica com campo vetorial por tetraedro em   */
/* formato VTK legacy ASCII. Tambem grava a magnitude como escalar por celula.*/
/* ENTRADA: path: const std::string &; mesh: const Mesh3D &; cell_vx: const   */
/* std::vector<double> &; cell_vy: const std::vector<double> &; cell_vz:      */
/* const std::vector<double> &; cell_mag: const std::vector<double> &;        */
/* vector_name: const std::string &; scalar_name: const std::string &.        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void write_vtk_unstructured_tet_cell_vector(
    const std::string &path,
    const Mesh3D &mesh,
    const std::vector<double> &cell_vx,
    const std::vector<double> &cell_vy,
    const std::vector<double> &cell_vz,
    const std::vector<double> &cell_mag,
    const std::string &vector_name = "E",
    const std::string &scalar_name = "Emag")
{
  std::ofstream f(path);
  f << "# vtk DataFile Version 3.0\n";
  f << "tp3485 tet-cell-vector\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  f << "POINTS " << mesh.nodes.size() << " double\n";
  for (const auto &n : mesh.nodes)
    f << n.x << " " << n.y << " " << n.z << "\n";

  const int ncells = (int)mesh.tets.size();
  f << "CELLS " << ncells << " " << 5 * ncells << "\n";
  for (const auto &t : mesh.tets)
    f << "4 " << t.v[0] << " " << t.v[1] << " " << t.v[2] << " " << t.v[3] << "\n";

  f << "CELL_TYPES " << ncells << "\n";
  for (int i = 0; i < ncells; ++i)
    f << "10\n"; // VTK_TETRA

  f << "CELL_DATA " << ncells << "\n";
  f << "SCALARS " << scalar_name << " double 1\n";
  f << "LOOKUP_TABLE default\n";
  for (int i = 0; i < ncells; ++i)
    f << cell_mag[(size_t)i] << "\n";

  f << "VECTORS " << vector_name << " double\n";
  for (int i = 0; i < ncells; ++i)
    f << cell_vx[(size_t)i] << " " << cell_vy[(size_t)i] << " " << cell_vz[(size_t)i] << "\n";
}
