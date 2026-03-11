/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_basis.cpp                                       */
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

#include "edge3d_basis.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

namespace
{
constexpr std::array<std::array<int, 2>, 6> kLocalEdges = {{
    {{0, 1}},
    {{0, 2}},
    {{0, 3}},
    {{1, 2}},
    {{1, 3}},
    {{2, 3}},
}};
// Ordem local adotada para as 6 arestas do tetraedro:
// (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).
// Esta ordem precisa ser consistente com a montagem local Sel/Tel.

/******************************************************************************/
/* FUNCAO: sub                                                                */
/* DESCRICAO: Subtrai dois vetores 3D componente a componente.                        */
/* ENTRADA: a: Vec3d; b: Vec3d.                                               */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d sub(Vec3d a, Vec3d b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

/******************************************************************************/
/* FUNCAO: dot                                                                */
/* DESCRICAO: Calcula produto escalar entre vetores 3D para compor normas e   */
/* integrais locais de massa/rigidez.                                         */
/* ENTRADA: a: Vec3d; b: Vec3d.                                               */
/* SAIDA: double.                                                             */
/******************************************************************************/
double dot(Vec3d a, Vec3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/******************************************************************************/
/* FUNCAO: cross                                                              */
/* DESCRICAO: Calcula produto vetorial 3D, usado no rotacional das bases de   */
/* Whitney e em operacoes geometricas no tetraedro.                           */
/* ENTRADA: a: Vec3d; b: Vec3d.                                               */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d cross(Vec3d a, Vec3d b)
{
  return {
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x};
}

/******************************************************************************/
/* FUNCAO: norm                                                               */
/* DESCRICAO: Retorna a norma Euclidiana de um vetor 3D.                              */
/* ENTRADA: a: Vec3d.                                                         */
/* SAIDA: double.                                                             */
/******************************************************************************/
double norm(Vec3d a) { return std::sqrt(dot(a, a)); }

/******************************************************************************/
/* FUNCAO: det3                                                               */
/* DESCRICAO: Calcula determinante 3x3, utilizado para Jacobiano e volume de tetraedros. */
/* ENTRADA: a11: double; a12: double; a13: double; a21: double; a22: double;  */
/* a23: double; a31: double; a32: double; a33: double.                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
double det3(
    double a11, double a12, double a13,
    double a21, double a22, double a23,
    double a31, double a32, double a33)
{
  return a11 * (a22 * a33 - a23 * a32) -
         a12 * (a21 * a33 - a23 * a31) +
         a13 * (a21 * a32 - a22 * a31);
}

/******************************************************************************/
/* FUNCAO: inverse4                                                           */
/* DESCRICAO: Inverte matriz 4x4 por Gauss-Jordan para recuperar coeficientes */
/* baricentricos usados nos gradientes de lambda_i.                           */
/* ENTRADA: A: const std::array<std::array<double, 4>, 4> &.                  */
/* SAIDA: std::array<std::array<double, 4>, 4>.                               */
/******************************************************************************/
std::array<std::array<double, 4>, 4> inverse4(const std::array<std::array<double, 4>, 4> &A)
{
  std::array<std::array<double, 8>, 4> aug{};
  for (int r = 0; r < 4; ++r)
  {
    for (int c = 0; c < 4; ++c)
      aug[r][c] = A[r][c];
    for (int c = 0; c < 4; ++c)
      aug[r][4 + c] = (r == c) ? 1.0 : 0.0;
  }

  for (int col = 0; col < 4; ++col)
  {
    int piv = col;
    double best = std::abs(aug[piv][col]);
    for (int r = col + 1; r < 4; ++r)
    {
      const double v = std::abs(aug[r][col]);
      if (v > best)
      {
        best = v;
        piv = r;
      }
    }
    if (best < 1e-15)
      throw std::runtime_error("Tetra degenerado: matriz singular ao calcular grad(lambda).");
    if (piv != col)
      std::swap(aug[piv], aug[col]);

    const double d = aug[col][col];
    for (int c = 0; c < 8; ++c)
      aug[col][c] /= d;

    for (int r = 0; r < 4; ++r)
    {
      if (r == col)
        continue;
      const double f = aug[r][col];
      if (std::abs(f) < 1e-18)
        continue;
      for (int c = 0; c < 8; ++c)
        aug[r][c] -= f * aug[col][c];
    }
  }

  std::array<std::array<double, 4>, 4> inv{};
  for (int r = 0; r < 4; ++r)
  {
    for (int c = 0; c < 4; ++c)
      inv[r][c] = aug[r][4 + c];
  }
  return inv;
}
} // namespace

/******************************************************************************/
/* FUNCAO: tet_geom_edge                                                      */
/* DESCRICAO: Extrai geometria local do tetraedro e constroi estruturas para   */
/* avaliacao das bases de Whitney; esta etapa alimenta diretamente os termos   */
/* de rigidez (Eq. 176) e massa (Eq. 177).                                    */
/* ENTRADA: mesh: const Mesh3D &; t: const Tet &.                             */
/* SAIDA: TetGeomEdge.                                                        */
/******************************************************************************/
TetGeomEdge tet_geom_edge(const Mesh3D &mesh, const Tet &t)
{
  TetGeomEdge tg;
  for (int i = 0; i < 4; ++i)
  {
    const Node3D &n = mesh.nodes[t.v[i]];
    tg.X[i] = {n.x, n.y, n.z};
  }

  Vec3d e1 = sub(tg.X[1], tg.X[0]);
  Vec3d e2 = sub(tg.X[2], tg.X[0]);
  Vec3d e3 = sub(tg.X[3], tg.X[0]);
  const double detJ = det3(
      e1.x, e2.x, e3.x,
      e1.y, e2.y, e3.y,
      e1.z, e2.z, e3.z);
  tg.V = std::abs(detJ) / 6.0;
  if (tg.V < 1e-18)
    throw std::runtime_error("Tetra degenerado: volume ~ 0.");

  std::array<std::array<double, 4>, 4> A = {{
      {{tg.X[0].x, tg.X[0].y, tg.X[0].z, 1.0}},
      {{tg.X[1].x, tg.X[1].y, tg.X[1].z, 1.0}},
      {{tg.X[2].x, tg.X[2].y, tg.X[2].z, 1.0}},
      {{tg.X[3].x, tg.X[3].y, tg.X[3].z, 1.0}},
  }};
  const auto Ainv = inverse4(A);

  // lambda_i = a_i x + b_i y + c_i z + d_i,
  // com grad(lambda_i) = (a_i,b_i,c_i) extraido da coluna i de A^{-1}.
  for (int i = 0; i < 4; ++i)
    tg.grad_lambda[i] = {Ainv[0][i], Ainv[1][i], Ainv[2][i]};

  for (int m = 0; m < 6; ++m)
  {
    const int i = kLocalEdges[m][0];
    const int j = kLocalEdges[m][1];
    tg.L[m] = norm(sub(tg.X[i], tg.X[j]));
  }

  return tg;
}

/******************************************************************************/
/* FUNCAO: whitney_W_local_3d                                                 */
/* DESCRICAO: Avalia base vetorial de Whitney no tetraedro local; contribuicao */
/* usada nas integrais de massa vetorial e nos coeficientes I5..I10.          */
/* ENTRADA: m: int; tg: const TetGeomEdge &; lambda: const std::array<double, */
/* 4> &.                                                                      */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d whitney_W_local_3d(int m, const TetGeomEdge &tg, const std::array<double, 4> &lambda)
{
  if (m < 0 || m >= 6)
    throw std::runtime_error("Indice de aresta local invalido em whitney_W_local_3d.");
  const int i = kLocalEdges[m][0];
  const int j = kLocalEdges[m][1];
  const Vec3d gi = tg.grad_lambda[i];
  const Vec3d gj = tg.grad_lambda[j];
  const double Lij = tg.L[m];
  // W_ij = L_ij * (lambda_i * grad(lambda_j) - lambda_j * grad(lambda_i)).
  return {
      Lij * (lambda[i] * gj.x - lambda[j] * gi.x),
      Lij * (lambda[i] * gj.y - lambda[j] * gi.y),
      Lij * (lambda[i] * gj.z - lambda[j] * gi.z)};
}

/******************************************************************************/
/* FUNCAO: whitney_curl_local_3d                                              */
/* DESCRICAO: Calcula rotacional local constante da base de Whitney no         */
/* tetraedro; termo usado no bloco curl-curl (Eq. 176, associado a I1..I4).   */
/* ENTRADA: m: int; tg: const TetGeomEdge &.                                  */
/* SAIDA: Vec3d.                                                              */
/******************************************************************************/
Vec3d whitney_curl_local_3d(int m, const TetGeomEdge &tg)
{
  if (m < 0 || m >= 6)
    throw std::runtime_error("Indice de aresta local invalido em whitney_curl_local_3d.");
  const int i = kLocalEdges[m][0];
  const int j = kLocalEdges[m][1];
  const Vec3d gi = tg.grad_lambda[i];
  const Vec3d gj = tg.grad_lambda[j];
  const Vec3d cij = cross(gi, gj);
  // curl( L_ij * (lambda_i grad(lambda_j) - lambda_j grad(lambda_i)) )
  // = 2 L_ij (grad(lambda_i) x grad(lambda_j)).
  return {2.0 * tg.L[m] * cij.x, 2.0 * tg.L[m] * cij.y, 2.0 * tg.L[m] * cij.z};
}
