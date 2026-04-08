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
#include "core/execution_log.hpp"
#include "core/lapack_eig.hpp"
#include "core/run_timing_fem3d_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "edge3d/edge3d_assembly.hpp"
#include "fem3d/fem3d_case_driver.hpp"
#include "fem3d/fem3d_case_output.hpp"
#include "fem3d/fem3d_compare.hpp"
#include "fem3d/fem3d_debug.hpp"
#include "fem3d/fem3d_field_output.hpp"
#include "fem3d0/fem3d0_driver_entry.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
inline std::vector<double> k0_values_from_candidates(
    const std::vector<fem3d::PositiveK0Candidate> &candidates)
{
  std::vector<double> out;
  out.reserve(candidates.size());
  for (const auto &cand : candidates)
    out.push_back(cand.k0);
  return out;
}

inline std::vector<double> matched_k0_from_indices(
    const std::vector<double> &lambda,
    const std::vector<int> &matched_indices,
    double tol = 1e-10)
{
  std::vector<double> out(matched_indices.size(), std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < matched_indices.size(); ++i)
  {
    const int eig_index = matched_indices[i];
    if (eig_index < 0 || eig_index >= (int)lambda.size())
      continue;
    const double l = lambda[(size_t)eig_index];
    if (l <= tol)
      continue;
    out[i] = std::sqrt(l);
  }
  return out;
}

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
  const auto dirs = fem3d_output::ensure_case_dirs("fem3d0", c.case_name);
  execution_log::ExecutionLogScope case_log((dirs.root / "run.log").string());
  if (!case_log.active())
  {
    std::cerr << "[warn] nao foi possivel abrir run.log do caso "
              << c.case_name << ": " << case_log.error_message() << "\n";
  }

  timing::Breakdown perf;
  timing::Stopwatch total_watch;
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
  std::cout << "output_root=" << dirs.root.string() << "\n";
  std::cout << "csv_dir=" << dirs.csv.string() << "\n";
  std::cout << "linop_dir=" << dirs.linop.string() << "\n";
  std::cout << "vtk_dir=" << dirs.vtk.string() << "\n";
  std::cout << "nodes=" << c.mesh.nodes.size()
            << ", tets=" << c.mesh.tets.size()
            << ", edges=" << sys.ed.edges.size()
            << ", dof=" << sys.ed.ndof
            << ", backend=" << element_assembly_backend_name(backend) << "\n";

  stage.reset();
  const auto eig = generalized_eigs_sym_vec(sys.S, sys.T);
  perf.solve_ms += stage.elapsed_ms();

  stage.reset();
  const int scan_limit = fem3d::scan_limit_for_table((int)c.rows.size());
  const auto k0_candidates = fem3d::first_positive_k0_candidates(eig.w, scan_limit);
  if (debug_candidates)
    fem3d::print_positive_k0_candidates_debug(k0_values_from_candidates(k0_candidates));
  const auto matched_indices =
      fem3d::match_indices_by_reference_with_degeneracy(c.rows, k0_candidates, scan_limit);
  const auto k0_match = matched_k0_from_indices(eig.w, matched_indices);
  fem3d::print_table_compare(
      "Comparacao (casamento agrupado para raizes analiticas degeneradas)",
      c.rows,
      k0_match);

  const std::string prefix = "fem3d0_" + c.case_name;
  auto mode_rows = fem3d_output::build_mode_records(c.rows, k0_match);
  for (size_t i = 0; i < mode_rows.size(); ++i)
  {
    const int eig_index = (i < matched_indices.size()) ? matched_indices[i] : -1;
    mode_rows[i].matched_eig_index = eig_index;
    if (eig_index < 0)
    {
      mode_rows[i].field_status = "missing";
      continue;
    }

    const auto artifacts = fem3d_field::export_mode(
        dirs,
        c.mesh,
        sys.ed,
        eig.Zcol,
        prefix,
        eig_index,
        static_cast<int>(i) + 1,
        mode_rows[i].mode_label);
    mode_rows[i].field_status = artifacts.field_status;
    mode_rows[i].fields_csv_file = artifacts.fields_csv_file;
    mode_rows[i].vtk_file = artifacts.vtk_file;
  }

  const auto modes_path = dirs.csv / (prefix + "_modes.csv");
  if (!fem3d_output::write_modes_csv(modes_path, mode_rows))
    throw std::runtime_error("Falha ao escrever CSV de modos: " + modes_path.string());

  if (!spectral_csv::write_symmetric_problem_exports(dirs.linop, prefix, sys.S, sys.T, eig))
    throw std::runtime_error("Falha ao exportar artefatos espectrais do caso " + c.case_name);

  perf.post_ms += stage.elapsed_ms();
  perf.total_ms = total_watch.elapsed_ms();

  run_timing_fem3d_csv::Record timing_record;
  timing_record.solver = "fem3d0";
  timing_record.case_name = c.case_name;
  timing_record.storage = "dense";
  timing_record.backend = element_assembly_backend_name(backend);
  timing_record.debug_local_blocks = debug_local_blocks ? 1 : 0;
  timing_record.debug_candidates = debug_candidates ? 1 : 0;
  timing_record.nx = c.grid.nx;
  timing_record.ny = c.grid.ny;
  timing_record.nz = c.grid.nz;
  timing_record.mesh_nodes = static_cast<int>(c.mesh.nodes.size());
  timing_record.mesh_tets = static_cast<int>(c.mesh.tets.size());
  timing_record.mesh_edges = static_cast<int>(sys.ed.edges.size());
  timing_record.ndof = sys.ed.ndof;
  timing_record.assembly_ms = perf.assembly_ms;
  timing_record.solve_ms = perf.solve_ms;
  timing_record.post_ms = perf.post_ms;
  timing_record.total_ms = perf.total_ms;

  const auto run_timing_path = dirs.root / "run_timing.csv";
  if (!run_timing_fem3d_csv::write_csv(run_timing_path.string(), timing_record))
    throw std::runtime_error("Falha ao escrever run_timing.csv do caso " + c.case_name);

  std::cout << "saved_modes_csv=" << modes_path.string() << "\n";
  std::cout << "saved_run_timing=" << run_timing_path.string() << "\n";
  timing::print_breakdown("fem3d0_" + c.case_name, perf);
  return perf;
}
} // namespace

/******************************************************************************/
/* FUNCAO: run_fem3d0_with_options                                            */
/* DESCRICAO: Executa o fluxo principal do FEM3D0 a partir de opcoes de CLI   */
/* ja validadas, permitindo reutilizar o mesmo nucleo em drivers dedicados    */
/* por caso e no executavel legado multicase.                                 */
/* ENTRADA: opt: const fem3d::CliOptions &; banner_name: const char *.        */
/* SAIDA: int.                                                                */
/******************************************************************************/
int run_fem3d0_with_options(const fem3d::CliOptions &opt, const char *banner_name)
{
  timing::Breakdown perf;
  timing::Stopwatch total_watch;

  const auto solver_dirs = fem3d_output::ensure_solver_dirs("fem3d0");
  execution_log::ExecutionLogScope solver_log((solver_dirs.root / "run.log").string());
  if (!solver_log.active())
  {
    std::cerr << "[warn] nao foi possivel abrir run.log global do FEM3D0: "
              << solver_log.error_message() << "\n";
  }

  fem3d::for_each_selected_case(
      opt,
      "FEM3D0",
      [backend = opt.backend, debug_local_blocks = opt.debug_local_blocks, debug_candidates = opt.debug_candidates, &perf](const fem3d::PreparedCase &c)
      {
        const auto case_perf = run_dense_case(c, backend, debug_local_blocks, debug_candidates);
        perf.assembly_ms += case_perf.assembly_ms;
        perf.solve_ms += case_perf.solve_ms;
        perf.post_ms += case_perf.post_ms;
      });

  perf.total_ms = total_watch.elapsed_ms();
  timing::print_breakdown(banner_name, perf);
  return 0;
}

/******************************************************************************/
/* FUNCAO: run_fem3d0_case_driver                                             */
/* DESCRICAO: Entrada compartilhada para executaveis dedicados a um unico     */
/* caso 3D, com CLI restrita a malha/backend/depuracao.                       */
/* ENTRADA: argc: int; argv: char **; fixed_case: fem3d::CaseId; bin_name.    */
/* SAIDA: int.                                                                */
/******************************************************************************/
int run_fem3d0_case_driver(int argc, char **argv, fem3d::CaseId fixed_case, const char *bin_name)
{
  const auto opt = fem3d::parse_single_case_cli(argc, argv, fixed_case, bin_name);
  if (!opt)
    return 0;
  return run_fem3d0_with_options(*opt, bin_name);
}

/******************************************************************************/
/* FUNCAO: main                                                               */
/* DESCRICAO: Ponto de entrada do executavel: interpreta argumentos, prepara o*/
/* caso e executa o fluxo numerico principal para validacao da Secao 3.1.5.   */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: int.                                                                */
/******************************************************************************/
#ifndef FEM3D0_NO_LEGACY_MAIN
int main(int argc, char **argv)
{
  const fem3d::CliDefaults defaults{/*run_air=*/true, /*run_half=*/false, /*run_cyl=*/false, /*run_sphere=*/false};
  const auto opt = fem3d::parse_cli(argc, argv, defaults, "fem3d0_rect");
  if (!opt)
    return 0;
  return run_fem3d0_with_options(*opt, "fem3d0_rect");
}
#endif
