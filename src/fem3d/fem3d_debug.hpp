/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_debug.hpp                                         */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de depuracao para os casos vetoriais 3D. */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (178), Eq.      */
/* (181) e Eq. (182).                                                         */
/*****************************************************************************/

#pragma once

#include "core/mesh3d.hpp"
#include "edge3d/edge3d_basis.hpp"
#include "explicit/tet3d_edge_explicit.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace fem3d
{

/******************************************************************************/
/* FUNCAO: print_block_6x6                                                    */
/* DESCRICAO: Imprime uma matriz 6x6 em formato compacto para inspecao dos    */
/* blocos locais do tetraedro de aresta.                                      */
/* ENTRADA: name: const char *; M: const double[6][6].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_block_6x6(const char *name, const double M[6][6])
{
  std::cout << "  " << name << " =\n";
  for (int i = 0; i < 6; ++i)
  {
    std::cout << "   ";
    for (int j = 0; j < 6; ++j)
      std::cout << " " << M[i][j];
    std::cout << "\n";
  }
}

/******************************************************************************/
/* FUNCAO: print_first_tet_closed_form_debug                                  */
/* DESCRICAO: Imprime os blocos locais closed-form do primeiro tetraedro,     */
/* correspondentes as Eq. (181) e (182), para rastrear a contribuicao local   */
/* que depois alimenta a montagem global da Eq. (178).                        */
/* ENTRADA: mesh: const Mesh3D &; eps_r_tet: const std::vector<double> &;     */
/* mu_r_tet: const std::vector<double> &.                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_first_tet_closed_form_debug(
    const Mesh3D &mesh,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet)
{
  if (mesh.tets.empty() || eps_r_tet.empty() || mu_r_tet.empty())
    return;

  const Tet &t = mesh.tets.front();
  const TetGeomEdge tg = tet_geom_edge(mesh, t);
  double Sel[6][6] = {{0.0}};
  double Tel[6][6] = {{0.0}};
  explicit_tet3d::tet3d_edge_closed_form_eq_181_182(
      tg,
      1.0 / mu_r_tet.front(),
      eps_r_tet.front(),
      Sel,
      Tel);

  std::cout << "\n[debug] primeiro tetraedro: blocos locais closed-form 3D\n";
  std::cout << "  volume=" << tg.V
            << " eps_r=" << eps_r_tet.front()
            << " mu_r=" << mu_r_tet.front() << "\n";
  std::cout << "  artigo Eq. (181) e Eq. (182):\n";
  print_block_6x6("Sel", Sel);
  print_block_6x6("Tel", Tel);
}

/******************************************************************************/
/* FUNCAO: print_positive_k0_candidates_debug                                 */
/* DESCRICAO: Imprime uma amostra das primeiras raizes positivas k0 antes do  */
/* casamento com a tabela de referencia, para inspecao numerica do espectro.  */
/* ENTRADA: k0: const std::vector<double> &; max_print: int.                  */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_positive_k0_candidates_debug(
    const std::vector<double> &k0,
    int max_print = 20)
{
  std::cout << "\n[debug] primeiras raizes positivas k0 antes do matching:\n";
  const int n = std::min<int>((int)k0.size(), max_print);
  for (int i = 0; i < n; ++i)
    std::cout << "  " << (i + 1) << "  " << k0[(size_t)i] << "\n";
  if ((int)k0.size() > max_print)
    std::cout << "  ... (" << (k0.size() - (size_t)max_print) << " restantes)\n";
}

} // namespace fem3d
