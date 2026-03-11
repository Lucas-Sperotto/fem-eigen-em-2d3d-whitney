/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/io_vtk_sv.hpp                                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios de malha e saida para suporte aos experimentos 2D e */
/* 3D.                                                                        */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh2d.hpp"
#include <vector>
#include <fstream>
#include <string>

/******************************************************************************/
/* FUNCAO: write_vtk_unstructured_tri_scalar_vector                           */
/* DESCRICAO: Escreve malha triangular com campo escalar e vetorial em VTK.   */
/* ENTRADA: path: const std::string &; mesh: const Mesh2D &; nodal_scalar:    */
/* const std::vector<double> &; nodal_vx: const std::vector<double> &;        */
/* nodal_vy: const std::vector<double> &; scalar_name: const std::string &;   */
/* vector_name: const std::string &.                                          */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void write_vtk_unstructured_tri_scalar_vector(
    const std::string &path,
    const Mesh2D &mesh,
    const std::vector<double> &nodal_scalar,
    const std::vector<double> &nodal_vx,
    const std::vector<double> &nodal_vy,
    const std::string &scalar_name = "phi",
    const std::string &vector_name = "Ft")
{
  std::ofstream f(path);
  f << "# vtk DataFile Version 3.0\n";
  f << "tp3485 scalar+vector\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  f << "POINTS " << mesh.nodes.size() << " double\n";
  for (auto &n : mesh.nodes)
    f << n.x << " " << n.y << " 0\n";

  const int ncells = (int)mesh.tris.size();
  f << "CELLS " << ncells << " " << 4 * ncells << "\n";
  for (auto &t : mesh.tris)
  {
    f << "3 " << t.v[0] << " " << t.v[1] << " " << t.v[2] << "\n";
  }

  f << "CELL_TYPES " << ncells << "\n";
  for (int i = 0; i < ncells; i++)
    f << "5\n"; // VTK_TRIANGLE

  f << "POINT_DATA " << mesh.nodes.size() << "\n";

  // scalar
  f << "SCALARS " << scalar_name << " double 1\n";
  f << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < nodal_scalar.size(); i++)
  {
    f << nodal_scalar[i] << "\n";
  }

  // vector
  f << "VECTORS " << vector_name << " double\n";
  for (size_t i = 0; i < mesh.nodes.size(); i++)
  {
    f << nodal_vx[i] << " " << nodal_vy[i] << " 0\n";
  }
}

/******************************************************************************/
/* FUNCAO: write_vtk_unstructured_tri_cell_vector                             */
/* DESCRICAO: Escreve campo vetorial por celula em malha triangular no formato*/
/* VTK.                                                                       */
/* ENTRADA: path: const std::string &; mesh: const Mesh2D &; cell_vx: const   */
/* std::vector<double> &; cell_vy: const std::vector<double> &; vector_name:  */
/* const std::string &.                                                       */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void write_vtk_unstructured_tri_cell_vector(
    const std::string &path,
    const Mesh2D &mesh,
    const std::vector<double> &cell_vx,
    const std::vector<double> &cell_vy,
    const std::string &vector_name = "Et")
{
  std::ofstream f(path);
  f << "# vtk DataFile Version 3.0\n";
  f << "tp3485 cell-vector\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  f << "POINTS " << mesh.nodes.size() << " double\n";
  for (auto &n : mesh.nodes)
    f << n.x << " " << n.y << " 0\n";

  const int ncells = (int)mesh.tris.size();
  f << "CELLS " << ncells << " " << 4 * ncells << "\n";
  for (auto &t : mesh.tris)
  {
    f << "3 " << t.v[0] << " " << t.v[1] << " " << t.v[2] << "\n";
  }

  f << "CELL_TYPES " << ncells << "\n";
  for (int i = 0; i < ncells; i++)
    f << "5\n"; // VTK_TRIANGLE

  f << "CELL_DATA " << ncells << "\n";
  f << "VECTORS " << vector_name << " double\n";
  for (int i = 0; i < ncells; ++i)
  {
    f << cell_vx[(size_t)i] << " " << cell_vy[(size_t)i] << " 0\n";
  }
}
