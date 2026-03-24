/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/main_edge_coax.cpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
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

#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_coax.hpp"
#include "core/output_paths.hpp"
#include "core/timing_utils.hpp"
#include "edge/edge_assembly.hpp"
#include "edge/mode_match_coax_edge.hpp"
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

/******************************************************************************/
/* FUNCAO: run_case                                                           */
/* DESCRICAO: Executa o fluxo principal do caso (montagem, solve, comparacao analitica e escrita de saidas). */
/* ENTRADA: title: const char *; mesh: const Mesh2D &; r1: double; r2: double;*/
/* is_te: bool; bc: EdgeBC; out_dir: const std::filesystem::path &;           */
/* legacy_vtk_name: const char *; export_modes: int; backend:                 */
/* ElementAssemblyBackend; debug_candidates: bool.                            */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static timing::Breakdown run_case(
    const char *title,
    const Mesh2D &mesh,
    double r1,
    double r2,
    bool is_te,
    EdgeBC bc,
    const std::filesystem::path &out_dir,
    const char *legacy_vtk_name,
    int export_modes,
    ElementAssemblyBackend backend,
    bool debug_candidates)
{
    timing::Breakdown perf;
    std::cout << "\n[" << title << "]\n";

    timing::Stopwatch stage;
    const auto sys = build_helm10_edge_system(mesh, bc, 1.0, 1.0, backend);
    perf.assembly_ms += stage.elapsed_ms();
    std::cout << "edges=" << sys.ed.edges.size()
              << " edge_dofs=" << sys.ed.ndof << "\n";

    stage.reset();
    const auto res = generalized_eigs_sym_vec(sys.S, sys.T);
    perf.solve_ms += stage.elapsed_ms();

    if (debug_candidates)
        helmvec_debug::print_positive_kc_candidates_debug(res.w, 1e-9, 20);

    std::cout << "first kc:\n";
    helmvec_post::print_positive_kc(res.w, 12);

    std::cout << "\nTabela (" << title << "): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";

    int shown = 0;
    for (int i = 0; i < (int)res.w.size() && shown < 8; ++i)
    {
        if (res.w[(size_t)i] < 1e-9)
            continue;

        const double kc_fem = std::sqrt(res.w[(size_t)i]);
        const CoaxEdgeModeID id = is_te
                                      ? match_coax_edge_mode_by_mass_correlation_TE(mesh, r1, r2, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_coax_edge_mode_by_mass_correlation_TM(mesh, r1, r2, sys.T, sys.ed, res.Zcol, i, 8, 8);

        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;
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
        const CoaxEdgeModeID id = is_te
                                      ? match_coax_edge_mode_by_mass_correlation_TE(mesh, r1, r2, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_coax_edge_mode_by_mass_correlation_TM(mesh, r1, r2, sys.T, sys.ed, res.Zcol, i, 8, 8);

        std::vector<double> cell_vx;
        std::vector<double> cell_vy;
        helmvec_post::reconstruct_cell_field_from_edge_mode(mesh, sys, res.Zcol, i, cell_vx, cell_vy);

        std::ostringstream name;
        name << "edge_coax_" << (is_te ? "te" : "tm")
             << "_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported
             << "_" << (is_te ? "Et" : "Ht") << ".vtk";

        write_vtk_unstructured_tri_cell_vector(
            output_paths::file_in(out_dir, name.str()),
            mesh,
            cell_vx,
            cell_vy,
            is_te ? "Et" : "Ht");
        std::cout << "Saved: " << name.str() << " (CELL_DATA vectors)\n";

        if (exported == 1)
        {
            write_vtk_unstructured_tri_cell_vector(
                output_paths::file_in(out_dir, legacy_vtk_name),
                mesh,
                cell_vx,
                cell_vy,
                is_te ? "Et" : "Ht");
            std::cout << "Saved: " << legacy_vtk_name << " (CELL_DATA vectors)\n";
        }
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
    const double r1 = 0.3;
    const double r2 = 1.0;
    int nr = 12;
    int nt = 64;
    int export_modes = 8;
    helmvec::EdgeCliOptions cli;

    try
    {
        cli = helmvec::parse_edge_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./edge_coax [nr nt [nmodos]] [--backend gauss|closed-form]"
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

    const auto out_dir = output_paths::ensure_case_dir("2d/2.2.1_edge/coax");
    std::cout << "Output dir: " << out_dir << "\n";
    std::cout << "Backend de aresta: " << element_assembly_backend_name(cli.backend) << "\n";

    const auto mesh = make_coax_mesh(r1, r2, nr, nt);
    std::cout << "Edge coax: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size() << "\n";
    std::cout << "r1=" << r1 << " r2=" << r2
              << " nr=" << nr << " nt=" << nt << "\n";

    if (cli.debug_local_blocks)
        helmvec_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    const auto te_perf = run_case(
        "TE (Edge, PEC: Et=0 on boundary)",
        mesh,
        r1,
        r2,
        true,
        EdgeBC::TE_PEC_TangentialZero,
        out_dir,
        "edge_coax_Et.vtk",
        export_modes,
        cli.backend,
        cli.debug_candidates);
    perf.assembly_ms += te_perf.assembly_ms;
    perf.solve_ms += te_perf.solve_ms;
    perf.post_ms += te_perf.post_ms;

    const auto tm_perf = run_case(
        "TM (Edge, PEC: Hn=0, keep boundary edges)",
        mesh,
        r1,
        r2,
        false,
        EdgeBC::TM_PEC_NormalZero,
        out_dir,
        "edge_coax_Ht.vtk",
        export_modes,
        cli.backend,
        cli.debug_candidates);
    perf.assembly_ms += tm_perf.assembly_ms;
    perf.solve_ms += tm_perf.solve_ms;
    perf.post_ms += tm_perf.post_ms;

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("edge_coax", perf);

    return 0;
}
