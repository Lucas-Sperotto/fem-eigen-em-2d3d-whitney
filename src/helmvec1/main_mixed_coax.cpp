/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/main_mixed_coax.cpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema misto vetorial+escalar para kc, separando blocos        */
/* transverso/longitudinal.                                                   */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.2, Eq.    */
/* (92).                                                                      */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC1 do         */
/* apendice em FORTRAN. A montagem global correspondente a Eq. (92) ocorre em */
/* src/helmvec1/helmvec1_mixed_system.cpp. Comentarios priorizam didatica,    */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/error_metrics.hpp"
#include "core/mode_match_coax.hpp"
#include "core/execution_log.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_coax.hpp"
#include "core/output_paths.hpp"
#include "core/run_timing_mixed_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "mixed_coax_edge_match_wrapper.hpp"
#include "helmvec1_mixed_system.hpp"
#include "mixed_case_output.hpp"
#include "mixed_cli_options.hpp"
#include "mixed_debug.hpp"
#include "mixed_field_output.hpp"
#include "mixed_mode_match.hpp"
#include "mixed_mode_utils.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

namespace
{
constexpr int kMixedModeExportLimit = 20;

struct CoaxMatchedRow
{
    int positive_rank = 0;
    int m = 0;
    int p = 0;
    double kc_ana = 0.0;
    double kc_fem = 0.0;
    double err_percent = 0.0;
    double rho_abs = 0.0;
};

void print_coax_matched_table(
    const std::string &title,
    const std::vector<CoaxMatchedRow> &rows,
    double lambda_filter,
    int max_rows = 8)
{
    std::cout << "\n" << title << "\n";
    std::cout << "lambda_filter(kc^2) = " << lambda_filter << "\n";
    std::cout << " i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";
    int shown = 0;
    for (const CoaxMatchedRow &row : rows)
    {
        if (shown >= max_rows)
            break;
        std::cout << " " << row.positive_rank
                  << "  (" << row.m << "," << row.p << ")  "
                  << row.kc_ana << "  "
                  << row.kc_fem << "  "
                  << row.err_percent << "  "
                  << row.rho_abs << "\n";
        ++shown;
    }
}

void append_coax_mode_records(
    const char *formulation,
    const char *dominant_block,
    const char *component_label,
    const char *family,
    const std::vector<BlockEnergyMode> &modes,
    const char *case_name,
    const helmvec1_output::CaseDirs &dirs,
    const Mesh2D &mesh,
    const MixedSystem92 &sys,
    const std::vector<double> &mixed_zcol,
    bool use_edge_block,
    bool is_te_family,
    double r1,
    double r2,
    std::vector<helmvec1_output::ModeCsvRecord> &records,
    std::vector<CoaxMatchedRow> &matched_rows)
{
    for (size_t idx = 0; idx < modes.size(); ++idx)
    {
        if ((int)idx >= kMixedModeExportLimit)
            break;
        const BlockEnergyMode &mode = modes[idx];
        CoaxMatchedRow table_row;
        helmvec1_output::ModeCsvRecord rec;
        rec.formulation = formulation;
        rec.dominant_block = dominant_block;
        rec.component_label = component_label;
        rec.family = family;
        rec.positive_rank = static_cast<int>(idx) + 1;
        rec.eig_index = mode.eig_index;
        rec.r1_m = helmvec1_output::format_float(r1);
        rec.r2_m = helmvec1_output::format_float(r2);
        rec.kc_fem = helmvec1_output::format_float(mode.kc);
        rec.kc_r1_fem = helmvec1_output::format_float(mode.kc * r1);
        rec.match_space = use_edge_block ? "edge_Tt" : "scalar_Tz";
        rec.match_method = "mass_correlation";
        rec.edge_energy = mode.edge_energy;
        rec.scalar_energy = mode.scalar_energy;
        rec.dominant_energy_ratio = mode.dominant_energy_ratio;
        rec.mode_status = std::string(dominant_block) + "_dominant";

        if (use_edge_block)
        {
            const std::vector<double> block_zcol =
                helmvec1_match::extract_edge_block_mode_as_column_major(
                    mixed_zcol, sys.nt, sys.nz, mode.eig_index);
            const CoaxEdgeModeID id = is_te_family
                                          ? match_coax_edge_mode_by_mass_correlation_TE(
                                                mesh, r1, r2, sys.edge.T, sys.edge.ed, block_zcol, 0, 6, 6)
                                          : match_coax_edge_mode_by_mass_correlation_TM(
                                                mesh, r1, r2, sys.edge.T, sys.edge.ed, block_zcol, 0, 6, 6);
            rec.m = std::to_string(id.m);
            rec.p = std::to_string(id.p);
            rec.kc_ana = helmvec1_output::format_float(id.kc_ana);
            rec.kc_r1_ana = helmvec1_output::format_float(id.kc_ana * r1);
            const double err = error_metrics::absolute_relative_error_percent(id.kc_ana, mode.kc);
            rec.error_percent = helmvec1_output::format_float(err);
            rec.rho_abs = id.rho;
            rec.mode_label = std::string(family) + "_m" + rec.m + "_p" + rec.p;
            table_row = {rec.positive_rank, id.m, id.p, id.kc_ana, mode.kc, err, id.rho};

            const auto artifacts = helmvec1_field::export_edge_mode(
                case_name, dirs, mesh, sys, mixed_zcol, mode.eig_index, rec);
            rec.field_data_kind = artifacts.field_data_kind;
            rec.field_status = artifacts.field_status;
            rec.fields_csv_file = artifacts.fields_csv_file;
            rec.vtk_file = artifacts.vtk_file;
        }
        else
        {
            const std::vector<double> block_zcol =
                helmvec1_match::extract_scalar_block_mode_as_column_major(
                    mixed_zcol, sys.nt, sys.nz, mode.eig_index);
            const CoaxModeID id = match_coax_mode_by_mass_correlation(
                mesh,
                r1,
                r2,
                sys.scal,
                block_zcol,
                0,
                is_te_family ? ScalarBC::TE_Neumann : ScalarBC::TM_Dirichlet,
                6,
                6);
            rec.m = std::to_string(id.m);
            rec.p = std::to_string(id.p);
            rec.kc_ana = helmvec1_output::format_float(id.kc_ana);
            rec.kc_r1_ana = helmvec1_output::format_float(id.kc_ana * r1);
            const double err = error_metrics::absolute_relative_error_percent(id.kc_ana, mode.kc);
            rec.error_percent = helmvec1_output::format_float(err);
            rec.rho_abs = id.rho;
            rec.mode_label = std::string(family) + "_m" + rec.m + "_p" + rec.p;
            table_row = {rec.positive_rank, id.m, id.p, id.kc_ana, mode.kc, err, id.rho};

            const auto artifacts = helmvec1_field::export_scalar_mode(
                case_name, dirs, mesh, sys, mixed_zcol, mode.eig_index, rec);
            rec.field_data_kind = artifacts.field_data_kind;
            rec.field_status = artifacts.field_status;
            rec.fields_csv_file = artifacts.fields_csv_file;
            rec.vtk_file = artifacts.vtk_file;
        }

        records.push_back(rec);
        matched_rows.push_back(table_row);
    }
}

} // namespace

/******************************************************************************/
/* FUNCAO: main                                                               */
/* DESCRICAO: Ponto de entrada do executavel: interpreta argumentos, prepara o*/
/* caso e executa o fluxo numerico principal.                                 */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: int.                                                                */
/******************************************************************************/
int main(int argc, char **argv)
{
    timing::Breakdown perf;
    timing::Stopwatch total_watch;
    timing::Breakdown perf_e;
    timing::Breakdown perf_h;
    // Caso coaxial de referencia para a secao 2.2.2 (blocos escalar+vetorial homogeneos).
    const double r1 = 1.0;
    const double r2 = 4.0;
    int nr = 10, nt = 48;
    helmvec1::MixedCliOptions cli;
    try
    {
        cli = helmvec1::parse_mixed_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./mixed_coax [nr nt] [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (cli.positionals.size() >= 2)
    {
        nr = std::atoi(cli.positionals[0].c_str());
        nt = std::atoi(cli.positionals[1].c_str());
    }

    const auto dirs = helmvec1_output::ensure_case_dirs("coax");
    execution_log::ExecutionLogScope execution_log(
        output_paths::file_in(dirs.root, "run.log"));
    if (execution_log.active())
        std::cout << "Log file: " << execution_log.file_path() << "\n";
    else
        std::cerr << "Aviso: nao foi possivel criar log em "
                  << execution_log.file_path() << ": "
                  << execution_log.error_message() << "\n";
    std::cout << "Output dir: \"" << dirs.root.string() << "\"\n";
    std::cout << "CSV dir: \"" << dirs.csv.string() << "\"\n";
    std::cout << "VTK dir: \"" << dirs.vtk.string() << "\"\n";
    std::cout << "LinOp dir: \"" << dirs.linop.string() << "\"\n";

    Mesh2D mesh = make_coax_mesh(r1, r2, nr, nt);
    std::cout << "Mixed coax: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " r1=" << r1 << " r2=" << r2
              << " nr=" << nr << " nt=" << nt
              << " backend=" << element_assembly_backend_name(cli.backend) << "\n";

    std::vector<double> eps(mesh.tris.size(), 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);
    const double lambda_filter = 1e-10;
    std::vector<helmvec1_output::ModeCsvRecord> mode_records;

    if (cli.debug_local_blocks)
        helmvec1_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    timing::Stopwatch stage;
    auto sys_e = tp3485::build_eq92_helmvec1_system_E(mesh, eps, mu, cli.backend);
    perf_e.assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += perf_e.assembly_ms;
    stage.reset();
    auto res_e = generalized_eigs_sym_vec(sys_e.S, sys_e.T);
    perf_e.solve_ms = stage.elapsed_ms();
    perf.solve_ms += perf_e.solve_ms;
    if (!spectral_csv::write_symmetric_problem_exports(dirs.linop, "mixed_coax_E", sys_e.S, sys_e.T, res_e) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_E_St_crs.csv", sys_e.St) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_E_Tt_crs.csv", sys_e.Tt) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_E_Sz_crs.csv", sys_e.Sz) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_E_Tz_crs.csv", sys_e.Tz))
    {
        std::cerr << "Erro ao escrever artefatos espectrais da formulacao E em " << dirs.linop << "\n";
        return 4;
    }
    stage.reset();
    std::vector<BlockEnergyMode> modes_te_edge_e;
    std::vector<BlockEnergyMode> modes_tm_scalar_e;
    collect_modes_by_block_energy(
        res_e,
        sys_e.nt,
        sys_e.nz,
        lambda_filter,
        modes_te_edge_e,
        modes_tm_scalar_e);
    const std::vector<double> k_te_edge_e = extract_kc_values(modes_te_edge_e);
    const std::vector<double> k_tm_scalar_e = extract_kc_values(modes_tm_scalar_e);
    if (cli.debug_candidates)
    {
        helmvec1_debug::print_block_candidates_debug(
            "[E] candidatos TE (bloco de aresta):",
            k_te_edge_e,
            "[E] candidatos TM (bloco escalar):",
            k_tm_scalar_e);
    }

    std::cout << "\n[E-formulation | coax]\n";
    std::vector<CoaxMatchedRow> matched_te_edge_e;
    std::vector<CoaxMatchedRow> matched_tm_scalar_e;
    append_coax_mode_records(
        "E", "edge", "Et", "TE",
        modes_te_edge_e, "coax", dirs, mesh, sys_e, res_e.Zcol, true, true, r1, r2,
        mode_records, matched_te_edge_e);
    append_coax_mode_records(
        "E", "scalar", "Ez", "TM",
        modes_tm_scalar_e, "coax", dirs, mesh, sys_e, res_e.Zcol, false, false, r1, r2,
        mode_records, matched_tm_scalar_e);
    print_coax_matched_table("TE (edge block)", matched_te_edge_e, lambda_filter);
    print_coax_matched_table("TM (scalar block)", matched_tm_scalar_e, lambda_filter);
    perf_e.post_ms = stage.elapsed_ms();
    perf.post_ms += perf_e.post_ms;

    stage.reset();
    auto sys_h = tp3485::build_eq92_helmvec1_system_H(mesh, eps, mu, cli.backend);
    perf_h.assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += perf_h.assembly_ms;
    stage.reset();
    auto res_h = generalized_eigs_sym_vec(sys_h.S, sys_h.T);
    perf_h.solve_ms = stage.elapsed_ms();
    perf.solve_ms += perf_h.solve_ms;
    if (!spectral_csv::write_symmetric_problem_exports(dirs.linop, "mixed_coax_H", sys_h.S, sys_h.T, res_h) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_H_St_crs.csv", sys_h.St) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_H_Tt_crs.csv", sys_h.Tt) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_H_Sz_crs.csv", sys_h.Sz) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_coax_H_Tz_crs.csv", sys_h.Tz))
    {
        std::cerr << "Erro ao escrever artefatos espectrais da formulacao H em " << dirs.linop << "\n";
        return 4;
    }
    stage.reset();
    std::vector<BlockEnergyMode> modes_tm_edge_h;
    std::vector<BlockEnergyMode> modes_te_scalar_h;
    collect_modes_by_block_energy(
        res_h,
        sys_h.nt,
        sys_h.nz,
        lambda_filter,
        modes_tm_edge_h,
        modes_te_scalar_h);
    const std::vector<double> k_tm_edge_h = extract_kc_values(modes_tm_edge_h);
    const std::vector<double> k_te_scalar_h = extract_kc_values(modes_te_scalar_h);
    if (cli.debug_candidates)
    {
        helmvec1_debug::print_block_candidates_debug(
            "[H] candidatos TM (bloco de aresta):",
            k_tm_edge_h,
            "[H] candidatos TE (bloco escalar):",
            k_te_scalar_h);
    }

    std::cout << "\n[H-formulation | coax]\n";
    std::vector<CoaxMatchedRow> matched_tm_edge_h;
    std::vector<CoaxMatchedRow> matched_te_scalar_h;
    append_coax_mode_records(
        "H", "edge", "Ht", "TM",
        modes_tm_edge_h, "coax", dirs, mesh, sys_h, res_h.Zcol, true, false, r1, r2,
        mode_records, matched_tm_edge_h);
    append_coax_mode_records(
        "H", "scalar", "Hz", "TE",
        modes_te_scalar_h, "coax", dirs, mesh, sys_h, res_h.Zcol, false, true, r1, r2,
        mode_records, matched_te_scalar_h);
    print_coax_matched_table("TM (edge block)", matched_tm_edge_h, lambda_filter);
    print_coax_matched_table("TE (scalar block)", matched_te_scalar_h, lambda_filter);
    perf_h.post_ms = stage.elapsed_ms();
    perf.post_ms += perf_h.post_ms;

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("mixed_coax", perf);

    const std::string modes_csv_path =
        output_paths::file_in(dirs.csv, "mixed_coax_modes.csv");
    if (helmvec1_output::write_modes_csv(dirs.csv / "mixed_coax_modes.csv", mode_records))
        std::cout << "Saved: " << modes_csv_path << "\n";
    else
        std::cerr << "Aviso: falha ao escrever " << modes_csv_path << "\n";

    run_timing_mixed_csv::Record timing_record;
    timing_record.case_label = "mixed_coax";
    timing_record.geometry = "coax";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.debug_local_blocks = cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cli.debug_candidates ? 1 : 0;
    timing_record.r1_m = helmvec1_output::format_float(r1);
    timing_record.r2_m = helmvec1_output::format_float(r2);
    timing_record.nr = std::to_string(nr);
    timing_record.nt = std::to_string(nt);
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.e_assembly_ms = perf_e.assembly_ms;
    timing_record.e_solve_ms = perf_e.solve_ms;
    timing_record.e_post_ms = perf_e.post_ms;
    timing_record.h_assembly_ms = perf_h.assembly_ms;
    timing_record.h_solve_ms = perf_h.solve_ms;
    timing_record.h_post_ms = perf_h.post_ms;
    timing_record.assembly_ms_total = perf.assembly_ms;
    timing_record.solve_ms_total = perf.solve_ms;
    timing_record.post_ms_total = perf.post_ms;
    timing_record.total_ms = perf.total_ms;
    const std::string timing_csv_path = output_paths::file_in(dirs.root, "run_timing.csv");
    if (run_timing_mixed_csv::write_csv(timing_csv_path, timing_record))
        std::cout << "Saved: " << timing_csv_path << "\n";
    else
        std::cerr << "Aviso: falha ao escrever " << timing_csv_path << "\n";

    return 0;
}
