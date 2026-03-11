/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d0/main_fem3d0_rect.cpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Driver FEM3D0 (montagem densa) para cavidades 3D.               */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, Tabelas  */
/* 12-15.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "core/lapack_eig.hpp"
#include "edge3d/edge3d_assembly.hpp"
#include "fem3d/fem3d_case_driver.hpp"
#include "fem3d/fem3d_compare.hpp"

#include <iostream>

namespace
{
/******************************************************************************/
/* FUNCAO: run_dense_case                                                     */
/* DESCRICAO: Executa um caso 3D com solver denso (FEM3D0), incluindo montagem, solucao e exportacao. */
/* Resolve S x = (k0^2) T x, com S/T oriundas da formulacao vetorial 3D       */
/* (Secao 3.1, integrais I1..I10).                                            */
/* ENTRADA: c: const fem3d::PreparedCase &.                                   */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void run_dense_case(const fem3d::PreparedCase &c)
{
  // Fluxo FEM3D0:
  // 1) monta S (curl-curl) e T (massa) em formato denso;
  // 2) resolve S x = lambda T x;
  // 3) converte lambda -> k0 = sqrt(lambda) e compara com as tabelas 12-15.
  const auto sys = build_helm3d_edge_system(
      c.mesh,
      Edge3DBC::PEC_TangentialZero,
      c.eps_r_tet,
      c.mu_r_tet);

  std::cout << "\n" << c.header << "\n";
  std::cout << "nodes=" << c.mesh.nodes.size()
            << ", tets=" << c.mesh.tets.size()
            << ", edges=" << sys.ed.edges.size()
            << ", dof=" << sys.ed.ndof << "\n";

  const auto eig = generalized_eigs_sym_vec(sys.S, sys.T);
  const int scan_limit = fem3d::scan_limit_for_table((int)c.rows.size());
  const auto k0 = fem3d::first_positive_k0(eig.w, scan_limit);
  const auto k0_match = fem3d::match_by_reference_with_degeneracy(c.rows, k0, scan_limit);
  fem3d::print_table_compare(
      "Comparacao (casamento agrupado para raizes analiticas degeneradas)",
      c.rows,
      k0_match);
}
} // namespace

/******************************************************************************/
/* FUNCAO: main                                                               */
/* DESCRICAO: Ponto de entrada do executavel: interpreta argumentos, prepara o*/
/* caso e executa o fluxo numerico principal para validacao da Secao 3.1.5.   */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: int.                                                                */
/******************************************************************************/
int main(int argc, char **argv)
{
  const fem3d::CliDefaults defaults{/*run_air=*/true, /*run_half=*/false, /*run_cyl=*/false, /*run_sphere=*/false};
  const auto opt = fem3d::parse_cli(argc, argv, defaults, "fem3d0_rect");
  if (!opt)
    return 0;

  fem3d::for_each_selected_case(
      *opt,
      "FEM3D0",
      [](const fem3d::PreparedCase &c)
      {
        run_dense_case(c);
      });

  return 0;
}
