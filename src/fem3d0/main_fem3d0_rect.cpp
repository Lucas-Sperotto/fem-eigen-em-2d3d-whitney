/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d0/main_fem3d0_rect.cpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Driver FEM3D0 (montagem densa) para cavidades 3D.               */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, Tabelas  */
/* 12-15.                                                                     */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa FEM3D0 do apendice  */
/* em FORTRAN. A montagem global correspondente a Eq. (178) ocorre em         */
/* src/edge3d/edge3d_assembly.cpp. Comentarios priorizam didatica,            */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/lapack_eig.hpp"
#include "core/timing_utils.hpp"
#include "edge3d/edge3d_assembly.hpp"
#include "fem3d/fem3d_case_driver.hpp"
#include "fem3d/fem3d_compare.hpp"
#include "fem3d/fem3d_debug.hpp"

#include <iostream>

namespace
{
/******************************************************************************/
/* FUNCAO: run_dense_case                                                     */
/* DESCRICAO: Executa um caso 3D com solver denso (FEM3D0), incluindo montagem, solucao e exportacao. */
/* Resolve S x = (k0^2) T x, com S/T oriundas da formulacao vetorial 3D       */
/* (Secao 3.1, integrais I1..I10).                                            */
/* ENTRADA: c: const fem3d::PreparedCase &; backend: ElementAssemblyBackend.  */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
timing::Breakdown run_dense_case(
    const fem3d::PreparedCase &c,
    ElementAssemblyBackend backend,
    bool debug_local_blocks,
    bool debug_candidates)
{
  timing::Breakdown perf;
  // Fluxo FEM3D0:
  // 1) monta S (curl-curl) e T (massa) em formato denso;
  // 2) resolve S x = lambda T x;
  // 3) converte lambda -> k0 = sqrt(lambda) e compara com as tabelas 12-15.
  timing::Stopwatch stage;
  const auto sys = tp3485::build_eq178_fem3d_system_dense(
      c.mesh,
      Edge3DBC::PEC_TangentialZero,
      c.eps_r_tet,
      c.mu_r_tet,
      backend);
  perf.assembly_ms += stage.elapsed_ms();

  if (debug_local_blocks)
    fem3d::print_first_tet_closed_form_debug(c.mesh, c.eps_r_tet, c.mu_r_tet);

  std::cout << "\n" << c.header << "\n";
  std::cout << "nodes=" << c.mesh.nodes.size()
            << ", tets=" << c.mesh.tets.size()
            << ", edges=" << sys.ed.edges.size()
            << ", dof=" << sys.ed.ndof
            << ", backend=" << element_assembly_backend_name(backend) << "\n";

  stage.reset();
  const auto eig = generalized_eigs_sym_vec(sys.S, sys.T);
  perf.solve_ms += stage.elapsed_ms();
  const int scan_limit = fem3d::scan_limit_for_table((int)c.rows.size());
  const auto k0 = fem3d::first_positive_k0(eig.w, scan_limit);
  if (debug_candidates)
    fem3d::print_positive_k0_candidates_debug(k0);
  const auto k0_match = fem3d::match_by_reference_with_degeneracy(c.rows, k0, scan_limit);
  fem3d::print_table_compare(
      "Comparacao (casamento agrupado para raizes analiticas degeneradas)",
      c.rows,
      k0_match);
  return perf;
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
  timing::Breakdown perf;
  timing::Stopwatch total_watch;
  const fem3d::CliDefaults defaults{/*run_air=*/true, /*run_half=*/false, /*run_cyl=*/false, /*run_sphere=*/false};
  const auto opt = fem3d::parse_cli(argc, argv, defaults, "fem3d0_rect");
  if (!opt)
    return 0;

  fem3d::for_each_selected_case(
      *opt,
      "FEM3D0",
      [backend = opt->backend, debug_local_blocks = opt->debug_local_blocks, debug_candidates = opt->debug_candidates, &perf](const fem3d::PreparedCase &c)
      {
        const auto case_perf = run_dense_case(c, backend, debug_local_blocks, debug_candidates);
        perf.assembly_ms += case_perf.assembly_ms;
        perf.solve_ms += case_perf.solve_ms;
      });

  perf.total_ms = total_watch.elapsed_ms();
  timing::print_breakdown("fem3d0_rect", perf);
  return 0;
}
