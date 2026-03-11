/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_basis.hpp                                       */
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

struct Vec3d
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct TetGeomEdge
{
  double V = 0.0;
  std::array<Vec3d, 4> X;
  std::array<Vec3d, 4> grad_lambda;
  // Ordem local das arestas: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).
  std::array<double, 6> L;
};

/******************************************************************************/
/* FUNCAO: tet_geom_edge                                                      */
/* DESCRICAO: Calcula dados geometricos locais do tetraedro (volume,          */
/* gradientes baricentricos e comprimentos de aresta) usados nas integrais    */
/* elementares da formulacao vetorial 3D (Secao 3.1).                         */
/* ENTRADA: mesh: const Mesh3D &; t: const Tet &.                             */
/* SAIDA: TetGeomEdge.                                                        */
/******************************************************************************/
TetGeomEdge tet_geom_edge(const Mesh3D &mesh, const Tet &t);
/******************************************************************************/
/* FUNCAO: whitney_W_local_3d                                                 */
/* DESCRICAO: Avalia a base vetorial de aresta W_ij = L_ij(la_i grad(la_j) -  */
/* la_j grad(la_i)) no tetraedro; termo usado no bloco de massa (Eq. 177).    */
/* ENTRADA: m: int; tg: const TetGeomEdge &; lambda: const std::array<double, */
/* 4> &.                                                                      */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d whitney_W_local_3d(int m, const TetGeomEdge &tg, const std::array<double, 4> &lambda);
/******************************************************************************/
/* FUNCAO: whitney_curl_local_3d                                              */
/* DESCRICAO: Avalia rotacional local da base de Whitney, usado no bloco      */
/* curl-curl da rigidez vetorial 3D (Eq. 176).                                */
/* ENTRADA: m: int; tg: const TetGeomEdge &.                                  */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d whitney_curl_local_3d(int m, const TetGeomEdge &tg);
