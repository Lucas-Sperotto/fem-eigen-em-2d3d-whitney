/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_basis.cpp                                           */
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

#include "edge_basis.hpp"
#include <cmath>

/******************************************************************************/
/* FUNCAO: grad_lambda                                                        */
/* DESCRICAO: Calcula os gradientes das funcoes baricentricas do triangulo, base da formulacao de aresta 2D. */
/* ENTRADA: i: int; g: const TriGeom &.                                       */
/* SAIDA: Vec2.                                                               */
/******************************************************************************/
static Vec2 grad_lambda(int i, const TriGeom &g)
{
  return {g.b[i] / (2.0 * g.A), g.c[i] / (2.0 * g.A)};
}

/******************************************************************************/
/* FUNCAO: dist                                                               */
/* DESCRICAO: Calcula distancia Euclidiana entre dois pontos 2D.              */
/* ENTRADA: a: Vec2; b: Vec2.                                                 */
/* SAIDA: double.                                                             */
/******************************************************************************/
static double dist(Vec2 a, Vec2 b)
{
  double dx = a.x - b.x, dy = a.y - b.y;
  return std::sqrt(dx * dx + dy * dy);
}

/******************************************************************************/
/* FUNCAO: tri_geom_edge                                                      */
/* DESCRICAO: Pre-processa geometria local do triangulo (area, gradientes e comprimentos) para bases de Whitney. */
/* ENTRADA: m: const Mesh2D &; t: const Tri &.                                */
/* SAIDA: TriGeomEdge.                                                        */
/******************************************************************************/
TriGeomEdge tri_geom_edge(const Mesh2D &m, const Tri &t)
{
  TriGeomEdge tg;
  tg.g = tri_geom(m, t);
  tg.X[0] = {m.nodes[t.v[0]].x, m.nodes[t.v[0]].y};
  tg.X[1] = {m.nodes[t.v[1]].x, m.nodes[t.v[1]].y};
  tg.X[2] = {m.nodes[t.v[2]].x, m.nodes[t.v[2]].y};

  tg.L[0] = dist(tg.X[0], tg.X[1]);
  tg.L[1] = dist(tg.X[1], tg.X[2]);
  tg.L[2] = dist(tg.X[2], tg.X[0]);
  return tg;
}

/******************************************************************************/
/* FUNCAO: whitney_W_local                                                    */
/* DESCRICAO: Avalia funcao de base de Whitney local no triangulo.            */
/* ENTRADA: m: int; tg: const TriGeomEdge &; lam: const std::array<double, 3> */
/* &.                                                                         */
/* SAIDA: Vec2.                                                               */
/******************************************************************************/
Vec2 whitney_W_local(int m, const TriGeomEdge &tg, const std::array<double, 3> &lam)
{
  // Mapeamento local das arestas para as 1-formas de Whitney.
  int i = (m == 0) ? 0 : (m == 1 ? 1 : 2);
  int j = (m == 0) ? 1 : (m == 1 ? 2 : 0);

  Vec2 gi = grad_lambda(i, tg.g);
  Vec2 gj = grad_lambda(j, tg.g);

  double Lij = tg.L[m];
  // Definicao da base de Whitney no triangulo:
  // W_ij = L_ij * (lambda_i * grad(lambda_j) - lambda_j * grad(lambda_i)).
  return {
      Lij * (lam[i] * gj.x - lam[j] * gi.x),
      Lij * (lam[i] * gj.y - lam[j] * gi.y)};
}

/******************************************************************************/
/* FUNCAO: cross2                                                             */
/* DESCRICAO: Calcula o pseudo-produto vetorial 2D (determinante) entre dois vetores planares. */
/* ENTRADA: a: Vec2; b: Vec2.                                                 */
/* SAIDA: double.                                                             */
/******************************************************************************/
static double cross2(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

/******************************************************************************/
/* FUNCAO: whitney_curl_local                                                 */
/* DESCRICAO: Avalia o rotacional escalar local da base de Whitney de aresta, */
/* usado no termo curl-curl da formulacao vetorial 2D.                        */
/* ENTRADA: m: int; tg: const TriGeomEdge &.                                  */
/* SAIDA: double.                                                             */
/******************************************************************************/
double whitney_curl_local(int m, const TriGeomEdge &tg)
{
  int i = (m == 0) ? 0 : (m == 1 ? 1 : 2);
  int j = (m == 0) ? 1 : (m == 1 ? 2 : 0);

  Vec2 gi = grad_lambda(i, tg.g);
  Vec2 gj = grad_lambda(j, tg.g);

  // Em elementos lineares:
  // curl_z(W_ij) = 2 * L_ij * (grad(lambda_i) x grad(lambda_j)).
  return 2.0 * tg.L[m] * cross2(gi, gj);
}
