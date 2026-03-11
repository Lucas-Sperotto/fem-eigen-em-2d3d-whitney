/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/main_edge_rect.cpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Executavel da formulacao vetorial transversal com elementos de  */
/* aresta.                                                                    */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/output_paths.hpp"
#include "edge/edge_assembly.hpp"
#include "edge/mode_match_rect_edge.hpp"
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
/* ENTRADA: tag: const char *; mesh: const Mesh2D &; a: double; b: double; bc:*/
/* EdgeBC; is_te: bool; out_dir: const std::filesystem::path &;               */
/* legacy_vtk_name: const char *; export_modes: int.                          */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static void run_case(
    const char *tag,
    const Mesh2D &mesh,
    double a,
    double b,
    EdgeBC bc,
    bool is_te,
    const std::filesystem::path &out_dir,
    const char *legacy_vtk_name,
    int export_modes)
{
    std::cout << "\n[" << tag << "]\n";
    const auto sys = build_helm10_edge_system(mesh, bc, 1.0, 1.0);
    std::cout << "edges=" << sys.ed.edges.size()
              << " edge_dofs=" << sys.ed.ndof << "\n";

    const auto res = generalized_eigs_sym_vec(sys.S, sys.T);

    std::cout << "first kc:\n";
    helmvec_post::print_positive_kc(res.w, 12);

    std::cout << "\nTabela (" << tag << "): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,n)   kc_ana      kc_fem      err(%)    |rho|\n";

    int shown = 0;
    for (int i = 0; i < (int)res.w.size() && shown < 8; ++i)
    {
        if (res.w[(size_t)i] < 1e-9)
            continue;

        const double kc_fem = std::sqrt(res.w[(size_t)i]);
        const RectEdgeModeID id = is_te
                                      ? match_rect_edge_mode_by_mass_correlation_te(mesh, a, b, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_rect_edge_mode_by_mass_correlation_tm(mesh, a, b, sys.T, sys.ed, res.Zcol, i, 8, 8);

        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;
        std::cout << std::setw(1) << (shown + 1) << "  (" << id.m << "," << id.n << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho << "\n";
        ++shown;
    }

    int exported = 0;
    for (int i = 0; i < (int)res.w.size() && exported < export_modes; ++i)
    {
        if (res.w[(size_t)i] < 1e-9)
            continue;

        ++exported;
        const RectEdgeModeID id = is_te
                                      ? match_rect_edge_mode_by_mass_correlation_te(mesh, a, b, sys.T, sys.ed, res.Zcol, i, 8, 8)
                                      : match_rect_edge_mode_by_mass_correlation_tm(mesh, a, b, sys.T, sys.ed, res.Zcol, i, 8, 8);

        std::vector<double> cell_vx;
        std::vector<double> cell_vy;
        helmvec_post::reconstruct_cell_field_from_edge_mode(mesh, sys, res.Zcol, i, cell_vx, cell_vy);

        std::ostringstream name;
        name << "edge_rect_" << (is_te ? "te" : "tm")
             << "_m" << id.m
             << "_n" << id.n
             << "_rank" << std::setw(2) << std::setfill('0') << exported
             << "_" << (is_te ? "Et" : "Ht") << ".vtk";

        write_vtk_unstructured_tri_cell_vector(
            output_paths::file_in(out_dir, name.str()),
            mesh,
            cell_vx,
            cell_vy,
            is_te ? "Et" : "Ht");
        std::cout << "Saved: " << name.str() << "\n";

        if (exported == 1)
        {
            write_vtk_unstructured_tri_cell_vector(
                output_paths::file_in(out_dir, legacy_vtk_name),
                mesh,
                cell_vx,
                cell_vy,
                is_te ? "Et" : "Ht");
            std::cout << "Saved: " << legacy_vtk_name << "\n";
        }
    }
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
    const double a = 1.0;
    const double b = 0.5;

    int nx = 14;
    int ny = 14;
    int export_modes = 8;
    if (argc >= 3)
    {
        nx = std::atoi(argv[1]);
        ny = std::atoi(argv[2]);
    }
    if (argc >= 4)
    {
        export_modes = std::max(1, std::atoi(argv[3]));
    }

    const auto out_dir = output_paths::ensure_case_dir("2d/2.2.1_edge/rect");
    std::cout << "Output dir: " << out_dir << "\n";

    const Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    std::cout << "Edge rect: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size() << "\n";

    run_case(
        "TE (Edge, PEC: Et=0 on boundary)",
        mesh,
        a,
        b,
        EdgeBC::TE_PEC_TangentialZero,
        true,
        out_dir,
        "edge_rect_Et.vtk",
        export_modes);

    run_case(
        "TM (Edge, PEC: Hn=0, keep boundary edges)",
        mesh,
        a,
        b,
        EdgeBC::TM_PEC_NormalZero,
        false,
        out_dir,
        "edge_rect_Ht.vtk",
        export_modes);

    return 0;
}
