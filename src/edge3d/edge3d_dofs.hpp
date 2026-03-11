/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_dofs.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 3D de elementos de aresta tetraedricos (Whitney 1-form) e*/
/* montagem.                                                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, integrais*/
/* I1..I10.                                                                   */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh3d.hpp"
#include <array>
#include <vector>

struct Edge3D
{
  int n0 = -1;
  int n1 = -1;
  bool is_boundary = false;
};

enum class Edge3DBC
{
  PEC_TangentialZero,
  KeepBoundary
};

struct TetEdges
{
  std::array<int, 6> e;
  std::array<int, 6> sgn;
};

struct EdgeDofs3D
{
  std::vector<Edge3D> edges;
  std::vector<TetEdges> tet_edges;
  std::vector<int> edge_to_dof;
  int ndof = 0;
};

/******************************************************************************/
/* FUNCAO: build_edge_dofs_3d                                                 */
/* DESCRICAO: Construi conectividade global de arestas, orientacoes locais e  */
/* mapeamento aresta->DOF para elementos de Whitney 3D.                       */
/* ENTRADA: mesh: const Mesh3D &.                                             */
/* SAIDA: EdgeDofs3D.                                                         */
/******************************************************************************/
EdgeDofs3D build_edge_dofs_3d(const Mesh3D &mesh);
/******************************************************************************/
/* FUNCAO: build_edge_dofs_3d                                                 */
/* DESCRICAO: Variante com controle de contorno: elimina DOFs de arestas de   */
/* fronteira para PEC (Et tangencial nulo) ou preserva todos os DOFs.         */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC.                               */
/* SAIDA: EdgeDofs3D.                                                         */
/******************************************************************************/
EdgeDofs3D build_edge_dofs_3d(const Mesh3D &mesh, Edge3DBC bc);
