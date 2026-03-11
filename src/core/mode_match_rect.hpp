/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mode_match_rect.hpp                                      */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo escalar 2D (malha, montagem e identificacao modal).      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh2d.hpp"
#include "dense.hpp"
#include "helm10_scalar_system.hpp"
#include <vector>
#include <cmath>
#include <utility>
#include <limits>
#include <algorithm>

struct RectModeID
{
  int m = 0, n = 0;
  double kc_ana = 0.0;
  double rho = 0.0; // |correlacao|
};

/******************************************************************************/
/* FUNCAO: kc_rect_analytic                                                   */
/* DESCRICAO: Retorna o valor analitico de k_c para o modo retangular (m,n), conforme separacao de variaveis. */
/* ENTRADA: a: double; b: double; m: int; n: int.                             */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double kc_rect_analytic(double a, double b, int m, int n)
{
  const double pi = 3.14159265358979323846;
  double kx = (m * pi) / a;
  double ky = (n * pi) / b;
  return std::sqrt(kx * kx + ky * ky);
}

/******************************************************************************/
/* FUNCAO: analytic_phi_rect_on_nodes                                         */
/* DESCRICAO: Avalia o potencial modal analitico retangular em todos os nos da*/
/* malha.                                                                     */
/* ENTRADA: mesh: const Mesh2D &; a: double; b: double; bc: ScalarBC; m: int; */
/* n: int; phi_nodes: std::vector<double> &.                                  */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void analytic_phi_rect_on_nodes(
    const Mesh2D &mesh, double a, double b, ScalarBC bc, int m, int n,
    std::vector<double> &phi_nodes)
{
  const double pi = 3.14159265358979323846;
  phi_nodes.assign(mesh.nodes.size(), 0.0);

  for (size_t i = 0; i < mesh.nodes.size(); i++)
  {
    double x = mesh.nodes[i].x;
    double y = mesh.nodes[i].y;

    if (bc == ScalarBC::TE_Neumann)
    {
      // cos cos (Neumann)
      phi_nodes[i] = std::cos(m * pi * x / a) * std::cos(n * pi * y / b);
    }
    else
    {
      // sin sin (Dirichlet)
      phi_nodes[i] = std::sin(m * pi * x / a) * std::sin(n * pi * y / b);
    }
  }
}

/******************************************************************************/
/* FUNCAO: restrict_nodes_to_dofs                                             */
/* DESCRICAO: Realiza transformacao de representacao de campos/modos para     */
/* analise e pos-processamento.                                               */
/* ENTRADA: phi_nodes: const std::vector<double> &; dof_map: const            */
/* std::vector<int> &; ndof: int; phi_dof: std::vector<double> &.             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void restrict_nodes_to_dofs(
    const std::vector<double> &phi_nodes,
    const std::vector<int> &dof_map,
    int ndof,
    std::vector<double> &phi_dof)
{
  phi_dof.assign((size_t)ndof, 0.0);
  for (size_t i = 0; i < dof_map.size(); i++)
  {
    int d = dof_map[i];
    if (d >= 0)
      phi_dof[d] = phi_nodes[i];
  }
}

/******************************************************************************/
/* FUNCAO: mass_inner                                                         */
/* DESCRICAO: Calcula produto interno ponderado pela matriz de massa, base    */
/* para correlacao e normalizacao modal.                                      */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &; y: const     */
/* std::vector<double> &.                                                     */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_inner(const DenseMat &T, const std::vector<double> &x, const std::vector<double> &y)
{
  // x^T T y
  const int n = T.n;
  double s = 0.0;
  for (int i = 0; i < n; i++)
  {
    double row = 0.0;
    const double xi = x[(size_t)i];
    for (int j = 0; j < n; j++)
    {
      row += T(i, j) * y[(size_t)j];
    }
    s += xi * row;
  }
  return s;
}

/******************************************************************************/
/* FUNCAO: mass_norm                                                          */
/* DESCRICAO: Calcula a norma energetica induzida pela matriz de massa.       */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &.              */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_norm(const DenseMat &T, const std::vector<double> &x)
{
  double v = mass_inner(T, x, x);
  return (v > 0) ? std::sqrt(v) : 0.0;
}

/******************************************************************************/
/* FUNCAO: mass_correlation_abs                                               */
/* DESCRICAO: Calcula a correlacao modal absoluta ponderada pela massa entre  */
/* dois vetores de modo.                                                      */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &; y: const     */
/* std::vector<double> &.                                                     */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_correlation_abs(const DenseMat &T, const std::vector<double> &x, const std::vector<double> &y)
{
  double xy = mass_inner(T, x, y);
  double nx = mass_norm(T, x);
  double ny = mass_norm(T, y);
  if (nx <= 0 || ny <= 0)
    return 0.0;
  return std::abs(xy / (nx * ny));
}

/******************************************************************************/
/* FUNCAO: extract_mode_from_Zcol                                             */
/* DESCRICAO: Extrai um autovetor (coluna) do formato column-major para vetor */
/* de trabalho.                                                               */
/* ENTRADA: Zcol: const std::vector<double> &; ndof: int; mode_idx: int; x:   */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void extract_mode_from_Zcol(
    const std::vector<double> &Zcol, int ndof, int mode_idx,
    std::vector<double> &x)
{
  // Zcol: column-major ndof x ndof
  x.assign((size_t)ndof, 0.0);
  auto idx_col = [&](int i, int j)
  { return (size_t)j * ndof + i; };
  for (int i = 0; i < ndof; i++)
  {
    x[(size_t)i] = Zcol[idx_col(i, mode_idx)];
  }
}

/******************************************************************************/
/* FUNCAO: candidates_rect                                                    */
/* DESCRICAO: Gera conjunto de pares (m,n) fisicamente validos conforme BC,   */
/* respeitando dominios TE (inclui zero) e TM (estritamente positivos).       */
/* ENTRADA: bc: ScalarBC; mmax: int; nmax: int.                               */
/* SAIDA: std::vector<std::pair<int, int>>.                                   */
/******************************************************************************/
inline std::vector<std::pair<int, int>> candidates_rect(ScalarBC bc, int mmax, int nmax)
{
  std::vector<std::pair<int, int>> cand;
  if (bc == ScalarBC::TE_Neumann)
  {
    for (int m = 0; m <= mmax; m++)
    {
      for (int n = 0; n <= nmax; n++)
      {
        if (m == 0 && n == 0)
          continue;
        cand.push_back({m, n});
      }
    }
  }
  else
  {
    for (int m = 1; m <= mmax; m++)
    {
      for (int n = 1; n <= nmax; n++)
      {
        cand.push_back({m, n});
      }
    }
  }
  return cand;
}

/******************************************************************************/
/* FUNCAO: match_rect_mode_by_mass_correlation                                */
/* DESCRICAO: Identifica o modo analitico retangular mais correlacionado ao   */
/* modo FEM.                                                                  */
/* ENTRADA: mesh: const Mesh2D &; a: double; b: double; sys: const            */
/* ScalarSystem &; Zcol: const std::vector<double> &; mode_idx: // autovetores*/
/* FEM int; bc: ScalarBC; mmax: int; nmax: int.                               */
/* SAIDA: RectModeID.                                                         */
/******************************************************************************/
inline RectModeID match_rect_mode_by_mass_correlation(
    const Mesh2D &mesh, double a, double b,
    const ScalarSystem &sys,
    const std::vector<double> &Zcol, // autovetores FEM
    int mode_idx,
    ScalarBC bc,
    int mmax = 8, int nmax = 8)
{
  std::vector<double> x;
  extract_mode_from_Zcol(Zcol, sys.ndof, mode_idx, x);

  auto cand = candidates_rect(bc, mmax, nmax);

  RectModeID best;
  best.rho = -1.0;

  std::vector<double> phi_nodes, phi_dof;
  for (auto [m, n] : cand)
  {
    analytic_phi_rect_on_nodes(mesh, a, b, bc, m, n, phi_nodes);
    restrict_nodes_to_dofs(phi_nodes, sys.dof_map, sys.ndof, phi_dof);

    double rho = mass_correlation_abs(sys.T, x, phi_dof);
    if (rho > best.rho)
    {
      best = {m, n, kc_rect_analytic(a, b, m, n), rho};
    }
  }
  return best;
}
