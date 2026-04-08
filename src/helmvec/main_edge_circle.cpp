/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/main_edge_circle.cpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Executavel da formulacao vetorial transversal com elementos de  */
/* aresta.                                                                    */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC do apendice */
/* em FORTRAN. A montagem global correspondente a Eq. (65) ocorre em          */
/* src/edge/edge_assembly.cpp. Comentarios priorizam didatica,                */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/error_metrics.hpp"
#include "core/execution_log.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_circle.hpp"
#include "core/output_paths.hpp"
#include "core/run_timing_edge_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "edge/edge_assembly.hpp"
#include "edge/mode_match_circle_edge.hpp"
#include "helmvec/edge_case_output.hpp"
#include "helmvec/edge_cli_options.hpp"
#include "helmvec/edge_debug.hpp"
#include "helmvec/edge_mode_post.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct CircleModeCsvRecord
{
    std::string family;
    std::string transverse_label;
    std::string mode_label;
    int positive_rank = 0;
    int eig_index = 0;
    int m = 0;
    int p = 0;
    double r_m = 0.0;
    double kc_fem = 0.0;
    double kc_ana = 0.0;
    double kc_r_fem = 0.0;
    double kc_r_ana = 0.0;
    double error_percent = 0.0;
    double rho_abs = 0.0;
    std::string field_status;
    std::string fields_csv_file;
    std::string vtk_file;
};

/******************************************************************************/
/* FUNCAO: run_case_circle                                                    */
/* DESCRICAO: Executa o caso 2.2.1 circular completo: montagem, autovalores, pos-processamento e relatorio. */
/* ENTRADA: tag: const char *; mesh: const Mesh2D &; r: double; bc: EdgeBC;   */
/* is_te: bool; out_dir: const std::filesystem::path &; legacy_vtk_name: const*/
/* char *; export_modes: int; backend: ElementAssemblyBackend;                */
/* debug_candidates: bool.                                                    */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static timing::Breakdown run_case_circle(
    const char *tag,
    const Mesh2D &mesh,
    double r,
    EdgeBC bc,
    bool is_te,
    const helmvec_output::CaseDirs &dirs,
    const char *legacy_vtk_name,
    int export_modes,
    ElementAssemblyBackend backend,
    bool debug_candidates,
    std::vector<CircleModeCsvRecord> &mode_records)
{
    timing::Breakdown perf;
    std::cout << "\n[" << tag << "]\n";

    timing::Stopwatch stage;
    const auto sys = tp3485::build_eq65_helmvec_system(mesh, bc, 1.0, 1.0, backend);
    perf.assembly_ms += stage.elapsed_ms();
    std::cout << "edges=" << sys.ed.edges.size()
              << " edge_dofs=" << sys.ed.ndof << "\n";

    stage.reset();
    const auto res = generalized_eigs_sym_vec(sys.S, sys.T);
    perf.solve_ms += stage.elapsed_ms();

    const std::string branch_prefix =
        std::string("edge_circle_") + (is_te ? "te" : "tm");
    if (!spectral_csv::write_symmetric_problem_exports(
            dirs.linop,
            branch_prefix,
            sys.S,
            sys.T,
            res))
    {
        throw std::runtime_error("Erro ao escrever artefatos espectrais de " + branch_prefix);
    }

    if (debug_candidates)
        helmvec_debug::print_positive_kc_candidates_debug(res.w, 1e-9, 20);

    std::cout << "first kc:\n";
    helmvec_post::print_positive_kc(res.w, 12);

    std::cout << "\nTabela (" << tag << "): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";

    int shown = 0;
    for (int i = 0; i < (int)res.w.size() && shown < 8; ++i)
    {
        if (res.w[(size_t)i] < 1e-9)
            continue;

        const double kc_fem = std::sqrt(res.w[(size_t)i]);
        const CircEdgeModeID id = is_te
                                      ? match_circle_edge_mode_by_mass_correlation_TE(mesh, r, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_circle_edge_mode_by_mass_correlation_TM(mesh, r, sys.T, sys.ed, res.Zcol, i, 8, 8);

        const double err = error_metrics::absolute_relative_error_percent(id.kc_ana, kc_fem);
        std::cout << std::setw(1) << (shown + 1)
                  << "  (" << id.m << "," << id.p << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho
                  << "\n";
        ++shown;
    }

    stage.reset();
    int exported = 0;
    for (int i = 0; i < (int)res.w.size() && exported < export_modes; ++i)
    {
        if (res.w[(size_t)i] < 1e-9)
            continue;

        ++exported;
        const CircEdgeModeID id = is_te
                                      ? match_circle_edge_mode_by_mass_correlation_TE(mesh, r, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_circle_edge_mode_by_mass_correlation_TM(mesh, r, sys.T, sys.ed, res.Zcol, i, 8, 8);

        std::vector<double> cell_vx;
        std::vector<double> cell_vy;
        helmvec_post::reconstruct_cell_field_from_edge_mode(mesh, sys, res.Zcol, i, cell_vx, cell_vy);

        std::ostringstream name;
        name << "edge_circle_" << (is_te ? "te" : "tm")
             << "_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported
             << "_" << (is_te ? "Et" : "Ht") << ".vtk";

        std::ostringstream fields_name;
        fields_name << "edge_circle_fields_" << (is_te ? "TE" : "TM")
                    << "_m" << id.m
                    << "_p" << id.p
                    << "_rank" << std::setw(2) << std::setfill('0') << exported
                    << ".csv";

        const std::string vtk_filename = name.str();
        const std::string fields_csv_filename = fields_name.str();
        if (!helmvec_output::write_mode_fields_csv(
                dirs.csv / fields_csv_filename,
                mesh,
                cell_vx,
                cell_vy,
                is_te ? "E" : "H"))
        {
            throw std::runtime_error("Erro ao escrever CSV de campos: " + fields_csv_filename);
        }

        write_vtk_unstructured_tri_cell_vector(
            output_paths::file_in(dirs.vtk, vtk_filename),
            mesh,
            cell_vx,
            cell_vy,
            is_te ? "Et" : "Ht");
        std::cout << "Saved: " << fields_csv_filename << "\n";
        std::cout << "Saved: " << vtk_filename << " (CELL_DATA vectors)\n";

        if (exported == 1)
        {
            write_vtk_unstructured_tri_cell_vector(
                output_paths::file_in(dirs.vtk, legacy_vtk_name),
                mesh,
                cell_vx,
                cell_vy,
                is_te ? "Et" : "Ht");
            std::cout << "Saved: " << legacy_vtk_name << " (CELL_DATA vectors)\n";
        }

        CircleModeCsvRecord rec;
        rec.family = is_te ? "TE" : "TM";
        rec.transverse_label = is_te ? "Et" : "Ht";
        rec.mode_label = rec.family + "_m" + std::to_string(id.m) + "_p" + std::to_string(id.p);
        rec.positive_rank = exported;
        rec.eig_index = i;
        rec.m = id.m;
        rec.p = id.p;
        rec.r_m = r;
        rec.kc_fem = std::sqrt(res.w[(size_t)i]);
        rec.kc_ana = id.kc_ana;
        rec.kc_r_fem = rec.kc_fem * r;
        rec.kc_r_ana = rec.kc_ana * r;
        rec.error_percent = error_metrics::absolute_relative_error_percent(rec.kc_ana, rec.kc_fem);
        rec.rho_abs = id.rho;
        rec.field_status = "cell_centroid_unit_peak_normalized";
        rec.fields_csv_file = fields_csv_filename;
        rec.vtk_file = vtk_filename;
        mode_records.push_back(rec);
    }

    perf.post_ms += stage.elapsed_ms();
    return perf;
}

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
    const double r = 1.0;
    int nr = 8;
    int nt = 48;
    int export_modes = 20;
    helmvec::EdgeCliOptions cli;

    try
    {
        cli = helmvec::parse_edge_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./edge_circle [nr nt [nmodos]] [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (cli.positionals.size() >= 2)
    {
        nr = std::atoi(cli.positionals[0].c_str());
        nt = std::atoi(cli.positionals[1].c_str());
    }
    if (cli.positionals.size() >= 3)
    {
        export_modes = std::max(0, std::atoi(cli.positionals[2].c_str()));
    }

    const auto dirs = helmvec_output::ensure_case_dirs("circle");
    execution_log::ExecutionLogScope execution_log(
        output_paths::file_in(dirs.root, "run.log"));
    if (execution_log.active())
        std::cout << "Log file: " << execution_log.file_path() << "\n";
    else
        std::cerr << "Aviso: nao foi possivel criar log em "
                  << execution_log.file_path() << ": "
                  << execution_log.error_message() << "\n";

    std::cout << "Output dir: " << dirs.root << "\n";
    std::cout << "CSV dir: " << dirs.csv << "\n";
    std::cout << "VTK dir: " << dirs.vtk << "\n";
    std::cout << "LinOp dir: " << dirs.linop << "\n";
    std::cout << "Backend de aresta: " << element_assembly_backend_name(cli.backend) << "\n";

    const auto mesh = make_circle_mesh(r, nr, nt);
    std::cout << "Edge circle: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size() << "\n";
    std::cout << "R=" << r << " nr=" << nr << " nt=" << nt << "\n";

    if (cli.debug_local_blocks)
        helmvec_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    std::vector<CircleModeCsvRecord> mode_records;
    const auto te_perf = run_case_circle(
        "TE (Edge, PEC: Et=0 on boundary)",
        mesh,
        r,
        EdgeBC::TE_PEC_TangentialZero,
        true,
        dirs,
        "edge_circle_Et.vtk",
        export_modes,
        cli.backend,
        cli.debug_candidates,
        mode_records);
    perf.assembly_ms += te_perf.assembly_ms;
    perf.solve_ms += te_perf.solve_ms;
    perf.post_ms += te_perf.post_ms;

    const auto tm_perf = run_case_circle(
        "TM (Edge, PEC: Hn=0, keep boundary edges)",
        mesh,
        r,
        EdgeBC::TM_PEC_NormalZero,
        false,
        dirs,
        "edge_circle_Ht.vtk",
        export_modes,
        cli.backend,
        cli.debug_candidates,
        mode_records);
    perf.assembly_ms += tm_perf.assembly_ms;
    perf.solve_ms += tm_perf.solve_ms;
    perf.post_ms += tm_perf.post_ms;

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("edge_circle", perf);

    const std::string modes_csv_path = output_paths::file_in(dirs.csv, "edge_circle_modes.csv");
    std::ofstream modes_csv(modes_csv_path);
    if (!modes_csv)
    {
        std::cerr << "Erro ao abrir CSV de modos: " << modes_csv_path << "\n";
        return 1;
    }
    modes_csv << std::setprecision(16);
    modes_csv << "family,transverse_label,mode_label,positive_rank,eig_index,m,p,r_m,"
                 "kc_fem,kc_ana,kc_r_fem,kc_r_ana,error_percent,rho_abs,"
                 "field_status,fields_csv_file,vtk_file\n";
    for (const auto &rec : mode_records)
    {
        modes_csv << rec.family << ","
                  << rec.transverse_label << ","
                  << rec.mode_label << ","
                  << rec.positive_rank << ","
                  << rec.eig_index << ","
                  << rec.m << ","
                  << rec.p << ","
                  << rec.r_m << ","
                  << rec.kc_fem << ","
                  << rec.kc_ana << ","
                  << rec.kc_r_fem << ","
                  << rec.kc_r_ana << ","
                  << rec.error_percent << ","
                  << rec.rho_abs << ","
                  << rec.field_status << ","
                  << rec.fields_csv_file << ","
                  << rec.vtk_file << "\n";
    }
    modes_csv.close();
    std::cout << "Saved: " << modes_csv_path << "\n";

    run_timing_edge_csv::Record timing_record;
    timing_record.case_label = "edge_circle";
    timing_record.geometry = "circle";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.debug_local_blocks = cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cli.debug_candidates ? 1 : 0;
    timing_record.r_m = std::to_string(r);
    timing_record.nr = std::to_string(nr);
    timing_record.nt = std::to_string(nt);
    timing_record.nmodos = export_modes;
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.te_assembly_ms = te_perf.assembly_ms;
    timing_record.te_solve_ms = te_perf.solve_ms;
    timing_record.te_post_ms = te_perf.post_ms;
    timing_record.tm_assembly_ms = tm_perf.assembly_ms;
    timing_record.tm_solve_ms = tm_perf.solve_ms;
    timing_record.tm_post_ms = tm_perf.post_ms;
    timing_record.assembly_ms_total = perf.assembly_ms;
    timing_record.solve_ms_total = perf.solve_ms;
    timing_record.post_ms_total = perf.post_ms;
    timing_record.total_ms = perf.total_ms;
    const std::string timing_csv_path = output_paths::file_in(dirs.root, "run_timing.csv");
    if (run_timing_edge_csv::write_csv(timing_csv_path, timing_record))
        std::cout << "Saved: " << timing_csv_path << "\n";
    else
        std::cerr << "Aviso: falha ao escrever " << timing_csv_path << "\n";

    return 0;
}
