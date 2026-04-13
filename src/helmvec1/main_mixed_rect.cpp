/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/main_mixed_rect.cpp                                  */
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
#include "core/mode_match_rect.hpp"
#include "core/execution_log.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/output_paths.hpp"
#include "core/run_timing_mixed_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "mixed_rect_edge_match_wrapper.hpp"
#include "helmvec1_mixed_system.hpp"
#include "mixed_case_output.hpp"
#include "mixed_cli_options.hpp"
#include "mixed_debug.hpp"
#include "mixed_field_output.hpp"
#include "mixed_mode_match.hpp"
#include "mixed_mode_utils.hpp"
#include "mixed_rect_reference.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iostream>
#include <vector>

namespace
{
constexpr int kMixedModeExportLimit = 20;

struct RectMatchedRow
{
    int positive_rank = 0;
    int m = 0;
    int n = 0;
    double kc_ana = 0.0;
    double kc_fem = 0.0;
    double err_percent = 0.0;
    double rho_abs = 0.0;
};

void print_rect_matched_table(
    const std::string &title,
    const std::vector<RectMatchedRow> &rows,
    double lambda_filter,
    int max_rows = 10)
{
    std::cout << "\n" << title << "\n";
    std::cout << "lambda_filter(kc^2) = " << lambda_filter << "\n";
    std::cout << " i  (m,n)   kc_ana      kc_fem      err(%)    |rho|\n";
    int shown = 0;
    for (const RectMatchedRow &row : rows)
    {
        if (shown >= max_rows)
            break;
        std::cout << " " << row.positive_rank
                  << "  (" << row.m << "," << row.n << ")  "
                  << row.kc_ana << "  "
                  << row.kc_fem << "  "
                  << row.err_percent << "  "
                  << row.rho_abs << "\n";
        ++shown;
    }
}

void append_rect_mode_records(
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
    double a,
    double b,
    std::vector<helmvec1_output::ModeCsvRecord> &records,
    std::vector<RectMatchedRow> &matched_rows)
{
    for (size_t idx = 0; idx < modes.size(); ++idx)
    {
        if ((int)idx >= kMixedModeExportLimit)
            break;
        const BlockEnergyMode &mode = modes[idx];
        RectMatchedRow table_row;
        helmvec1_output::ModeCsvRecord rec;
        rec.formulation = formulation;
        rec.dominant_block = dominant_block;
        rec.component_label = component_label;
        rec.family = family;
        rec.positive_rank = static_cast<int>(idx) + 1;
        rec.eig_index = mode.eig_index;
        rec.ar_m = helmvec1_output::format_float(a);
        rec.b_m = helmvec1_output::format_float(b);
        rec.kc_fem = helmvec1_output::format_float(mode.kc);
        rec.kc_ar_fem = helmvec1_output::format_float(mode.kc * a);
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
            const RectEdgeModeID id = is_te_family
                                          ? match_rect_edge_mode_by_mass_correlation_te(
                                                mesh, a, b, sys.edge.T, sys.edge.ed, block_zcol, 0, 8, 8)
                                          : match_rect_edge_mode_by_mass_correlation_tm(
                                                mesh, a, b, sys.edge.T, sys.edge.ed, block_zcol, 0, 8, 8);
            rec.m = std::to_string(id.m);
            rec.n = std::to_string(id.n);
            rec.kc_ana = helmvec1_output::format_float(id.kc_ana);
            rec.kc_ar_ana = helmvec1_output::format_float(id.kc_ana * a);
            const double err = error_metrics::absolute_relative_error_percent(id.kc_ana, mode.kc);
            rec.error_percent = helmvec1_output::format_float(err);
            rec.rho_abs = id.rho;
            rec.mode_label = std::string(family) + "_m" + rec.m + "_n" + rec.n;

            table_row = {rec.positive_rank, id.m, id.n, id.kc_ana, mode.kc, err, id.rho};

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
            const RectModeID id = match_rect_mode_by_mass_correlation(
                mesh,
                a,
                b,
                sys.scal,
                block_zcol,
                0,
                is_te_family ? ScalarBC::TE_Neumann : ScalarBC::TM_Dirichlet,
                8,
                8);
            rec.m = std::to_string(id.m);
            rec.n = std::to_string(id.n);
            rec.kc_ana = helmvec1_output::format_float(id.kc_ana);
            rec.kc_ar_ana = helmvec1_output::format_float(id.kc_ana * a);
            const double err = error_metrics::absolute_relative_error_percent(id.kc_ana, mode.kc);
            rec.error_percent = helmvec1_output::format_float(err);
            rec.rho_abs = id.rho;
            rec.mode_label = std::string(family) + "_m" + rec.m + "_n" + rec.n;

            table_row = {rec.positive_rank, id.m, id.n, id.kc_ana, mode.kc, err, id.rho};

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
    // Malha padrao (2*nx*ny triangulos). Para comparacao em formato de tabela,
    // usamos uma discretizacao uniforme moderadamente refinada e expomos (nx, ny).
    int nx = 10, ny = 20;
    helmvec1::MixedCliOptions cli;
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso: ./mixed_rect [nx ny] [--backend closed-form|gauss|efgmi]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--nx NX] [--ny NY]"
           << " (nao misture com os posicionais principais)\n";
        os << "Compatibilidade: os posicionais principais continuam aceitos, mas estao deprecated;"
           << " prefira os aliases nomeados acima.\n";
    };
    if (helmvec1::mixed_cli_requests_help(argc, argv))
    {
        print_usage(std::cout);
        return 0;
    }
    try
    {
        cli = helmvec1::parse_mixed_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const bool has_named_rect_args =
        cli.nx_was_provided ||
        cli.ny_was_provided;

    if (cli.nr_was_provided || cli.nt_was_provided)
    {
        std::cerr << "Erro: mixed_rect nao aceita --nr/--nt; use --nx/--ny.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_rect_args && !cli.positionals.empty())
    {
        std::cerr << "Erro: nao misture aliases nomeados principais com argumentos posicionais principais em mixed_rect; escolha apenas um estilo de chamada.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_rect_args)
    {
        if (cli.nx_was_provided)
            nx = cli.nx;
        if (cli.ny_was_provided)
            ny = cli.ny;
    }
    else
    {
        if (!cli.positionals.empty() && cli.positionals.size() != 2)
        {
            std::cerr << "Erro: mixed_rect exige exatamente [nx ny] quando usa argumentos posicionais.\n";
            print_usage(std::cerr);
            return 2;
        }
        try
        {
            if (cli.positionals.size() == 2)
            {
                nx = helmvec1::parse_positive_mixed_cli_int(cli.positionals[0], "nx");
                ny = helmvec1::parse_positive_mixed_cli_int(cli.positionals[1], "ny");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Erro ao interpretar argumentos posicionais de mixed_rect: "
                      << e.what() << "\n";
            print_usage(std::cerr);
            return 2;
        }
    }

    const double a = 1.0;
    const double b = 0.5;
    const double lambda_filter = 1e-10; // filtro em lambda = kc^2

    const auto dirs = helmvec1_output::ensure_case_dirs("rect");
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
    if (!has_named_rect_args && !cli.positionals.empty())
    {
        std::cerr << "Aviso: os argumentos posicionais principais de mixed_rect continuam aceitos por compatibilidade, mas estao deprecated; prefira --nx/--ny.\n";
    }

    Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    std::vector<double> eps(mesh.tris.size(), 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);

    if (cli.debug_local_blocks)
        helmvec1_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    const int nprint = 10;
    // Formulacao E: x=[Et; Ez], Sx=(kc^2)Tx com blocos da Eq. (92).
    std::cout << "backend=" << element_assembly_backend_name(cli.backend) << "\n";
    std::cout << "Mixed rect: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " a=" << a << " b=" << b
              << " nx=" << nx << " ny=" << ny << "\n";

    std::vector<helmvec1_output::ModeCsvRecord> mode_records;

    timing::Stopwatch stage;
    auto sys_e = tp3485::build_eq92_helmvec1_system_E(mesh, eps, mu, cli.backend);
    perf_e.assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += perf_e.assembly_ms;
    stage.reset();
    auto res_e = generalized_eigs_sym_vec(sys_e.S, sys_e.T);
    perf_e.solve_ms = stage.elapsed_ms();
    perf.solve_ms += perf_e.solve_ms;
    if (!spectral_csv::write_symmetric_problem_exports(dirs.linop, "mixed_rect_E", sys_e.S, sys_e.T, res_e) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_E_St_crs.csv", sys_e.St) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_E_Tt_crs.csv", sys_e.Tt) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_E_Sz_crs.csv", sys_e.Sz) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_E_Tz_crs.csv", sys_e.Tz))
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

    std::vector<RectMatchedRow> matched_te_edge_e;
    std::vector<RectMatchedRow> matched_tm_scalar_e;
    append_rect_mode_records(
        "E", "edge", "Et", "TE",
        modes_te_edge_e, "rect", dirs, mesh, sys_e, res_e.Zcol, true, true, a, b,
        mode_records, matched_te_edge_e);
    append_rect_mode_records(
        "E", "scalar", "Ez", "TM",
        modes_tm_scalar_e, "rect", dirs, mesh, sys_e, res_e.Zcol, false, false, a, b,
        mode_records, matched_tm_scalar_e);
    print_rect_matched_table("[E-formulation] TE cutoffs (edge block)", matched_te_edge_e, lambda_filter);
    print_rect_matched_table("[E-formulation] TM cutoffs (scalar block)", matched_tm_scalar_e, lambda_filter);
    perf_e.post_ms = stage.elapsed_ms();
    perf.post_ms += perf_e.post_ms;

    // Formulacao H: operador constitutivo dual por troca eps/mu.
    stage.reset();
    auto sys_h = tp3485::build_eq92_helmvec1_system_H(mesh, eps, mu, cli.backend);
    perf_h.assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += perf_h.assembly_ms;
    stage.reset();
    auto res_h = generalized_eigs_sym_vec(sys_h.S, sys_h.T);
    perf_h.solve_ms = stage.elapsed_ms();
    perf.solve_ms += perf_h.solve_ms;
    if (!spectral_csv::write_symmetric_problem_exports(dirs.linop, "mixed_rect_H", sys_h.S, sys_h.T, res_h) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_H_St_crs.csv", sys_h.St) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_H_Tt_crs.csv", sys_h.Tt) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_H_Sz_crs.csv", sys_h.Sz) ||
        !spectral_csv::write_dense_crs_csv(dirs.linop / "mixed_rect_H_Tz_crs.csv", sys_h.Tz))
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

    std::vector<RectMatchedRow> matched_tm_edge_h;
    std::vector<RectMatchedRow> matched_te_scalar_h;
    append_rect_mode_records(
        "H", "edge", "Ht", "TM",
        modes_tm_edge_h, "rect", dirs, mesh, sys_h, res_h.Zcol, true, false, a, b,
        mode_records, matched_tm_edge_h);
    append_rect_mode_records(
        "H", "scalar", "Hz", "TE",
        modes_te_scalar_h, "rect", dirs, mesh, sys_h, res_h.Zcol, false, true, a, b,
        mode_records, matched_te_scalar_h);
    // Na formulacao H, o bloco de aresta corresponde a familia TM.
    print_rect_matched_table("[H-formulation] TM cutoffs (edge block)", matched_tm_edge_h, lambda_filter);
    // Na formulacao H, o bloco escalar corresponde a familia TE.
    print_rect_matched_table("[H-formulation] TE cutoffs (scalar block)", matched_te_scalar_h, lambda_filter);
    perf_h.post_ms = stage.elapsed_ms();
    perf.post_ms += perf_h.post_ms;

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("mixed_rect", perf);

    const std::string modes_csv_path =
        output_paths::file_in(dirs.csv, "mixed_rect_modes.csv");
    if (helmvec1_output::write_modes_csv(dirs.csv / "mixed_rect_modes.csv", mode_records))
        std::cout << "Saved: " << modes_csv_path << "\n";
    else
        std::cerr << "Aviso: falha ao escrever " << modes_csv_path << "\n";

    run_timing_mixed_csv::Record timing_record;
    timing_record.case_label = "mixed_rect";
    timing_record.geometry = "rect";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.debug_local_blocks = cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cli.debug_candidates ? 1 : 0;
    timing_record.ar_m = helmvec1_output::format_float(a);
    timing_record.b_m = helmvec1_output::format_float(b);
    timing_record.nx = std::to_string(nx);
    timing_record.ny = std::to_string(ny);
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
