/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_basis.hpp                                           */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 2D de elementos de aresta (DOFs, base de Whitney e       */
/* montagem).                                                                 */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh2d.hpp"
#include "core/fem_scalar.hpp"
#include <array>

struct Vec2
{
  double x = 0, y = 0;
};

struct TriGeomEdge
{
  TriGeom g;               // geometria escalar do triangulo (A, b[], c[])
  std::array<Vec2, 3> X;   // coordenadas dos vertices
  std::array<double, 3> L; // comprimentos locais das arestas: (0-1), (1-2), (2-0)
};

/******************************************************************************/
/* FUNCAO: tri_geom_edge                                                      */
/* DESCRICAO: Pre-processa geometria local do triangulo (area, gradientes e comprimentos) para bases de Whitney. */
/* ENTRADA: m: const Mesh2D &; t: const Tri &.                                */
/* SAIDA: TriGeomEdge.                                                        */
/******************************************************************************/
TriGeomEdge tri_geom_edge(const Mesh2D &m, const Tri &t);

// Convencao local para as funcoes de base de Whitney:
//   m=0 -> aresta (0->1), m=1 -> aresta (1->2), m=2 -> aresta (2->0).
Vec2 whitney_W_local(int m, const TriGeomEdge &tg, const std::array<double, 3> &lambda);

// Em 2D, curl_z(W_m) e constante dentro de triangulo linear.
double whitney_curl_local(int m, const TriGeomEdge &tg);
