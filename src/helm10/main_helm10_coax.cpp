/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_helm10_coax.cpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Executavel da formulacao escalar 2D para kc (guia               */
/* retangular/circular/coaxial).                                              */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1, Tabelas  */
/* 1-3.                                                                       */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELM10 do apendice  */
/* em FORTRAN. A montagem global correspondente a Eq. (43) ocorre em          */
/* src/core/helm10_scalar_system.cpp. Comentarios priorizam didatica,         */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "core/helm10_scalar_system.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_coax.hpp"
#include "core/mode_match_coax.hpp"
#include "core/output_paths.hpp"
#include "core/timing_utils.hpp"
#include "helm10/scalar_cli_options.hpp"
#include "helm10/scalar_debug.hpp"
#include "helm10/scalar_mode_post.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
    const double r1 = 1.0;
    const double r2 = 4.0;

    int nr = 10;
    int nt = 48;
    int export_modes = 8;
    helm10::ScalarCliOptions cli;
    try
    {
        cli = helm10::parse_scalar_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./helm10_coax [nr nt [nmodos]] [--backend gauss|closed-form]"
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

    const auto out_dir = output_paths::ensure_case_dir("2d/2.1_scalar/coax");
    std::cout << "Output dir: " << out_dir << "\n";
    std::cout << "Backend escalar: " << element_assembly_backend_name(cli.backend) << "\n";

    const Mesh2D mesh = make_coax_mesh(r1, r2, nr, nt);
    std::cout << "Coax mesh: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " r1=" << r1 << " r2=" << r2
              << " nr=" << nr << " nt=" << nt << "\n\n";

    if (cli.debug_local_blocks)
        helm10_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    // Modos TE da formulacao escalar com Neumann.
    timing::Stopwatch stage;
    const auto sys_te = build_helm10_scalar_system(mesh, ScalarBC::TE_Neumann, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    const auto te = generalized_eigs_sym_vec(sys_te.S, sys_te.T);
    perf.solve_ms += stage.elapsed_ms();
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(te.w, 1e-9);
    std::cout << "[TE scalar (Neumann) - coax]\n";
    helm10_post::print_positive_kc(te.w, 12, true);

    std::cout << "\nTabela 3 (TE coax): FEM vs Analitico (det J'/Y')\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";
    int shown = 0;
    for (int i = 0; i < (int)te.w.size() && shown < 8; ++i)
    {
        if (te.w[(size_t)i] < 1e-9)
            continue;

        const double kc_fem = std::sqrt(te.w[(size_t)i]);
        const auto id = match_coax_mode_by_mass_correlation(
            mesh,
            r1,
            r2,
            sys_te,
            te.Zcol,
            i,
            ScalarBC::TE_Neumann,
            6,
            6);
        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;

        std::cout << (shown + 1) << "  (" << id.m << "," << id.p << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho << "\n";
        ++shown;
    }

    // Modos TM da formulacao escalar com Dirichlet.
    stage.reset();
    const auto sys_tm = build_helm10_scalar_system(mesh, ScalarBC::TM_Dirichlet, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    const auto tm = generalized_eigs_sym_vec(sys_tm.S, sys_tm.T);
    perf.solve_ms += stage.elapsed_ms();
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(tm.w, 0.0);
    std::cout << "\n[TM scalar (Dirichlet) - coax]\n";
    helm10_post::print_positive_kc(tm.w, 12, false);

    std::cout << "\nTabela 3 (TM coax): FEM vs Analitico (det J/Y)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";
    shown = 0;
    for (int i = 0; i < (int)tm.w.size() && shown < 8; ++i)
    {
        if (tm.w[(size_t)i] < 0.0)
            continue;

        const double kc_fem = std::sqrt(tm.w[(size_t)i]);
        const auto id = match_coax_mode_by_mass_correlation(
            mesh,
            r1,
            r2,
            sys_tm,
            tm.Zcol,
            i,
            ScalarBC::TM_Dirichlet,
            6,
            6);
        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;

        std::cout << (shown + 1) << "  (" << id.m << "," << id.p << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho << "\n";
        ++shown;
    }

    auto write_mode = [&](const auto &sys,
                          const auto &res,
                          int eig_idx,
                          const std::string &vtk_name) {
        const auto phi = helm10_post::extract_mode_nodal_from_Z(mesh, sys, res.Zcol, eig_idx);
        std::vector<double> fx;
        std::vector<double> fy;
        helm10_post::compute_smoothed_transverse_from_scalar(mesh, phi, fx, fy);
        write_vtk_unstructured_tri_scalar_vector(
            output_paths::file_in(out_dir, vtk_name),
            mesh,
            phi,
            fx,
            fy,
            "phi",
            "Ft");
        std::cout << "Saved: " << vtk_name << "\n";
    };

    stage.reset();
    int exported_te = 0;
    for (int i = 0; i < (int)te.w.size() && exported_te < export_modes; ++i)
    {
        if (te.w[(size_t)i] < 1e-9)
            continue;

        ++exported_te;
        const auto id = match_coax_mode_by_mass_correlation(
            mesh,
            r1,
            r2,
            sys_te,
            te.Zcol,
            i,
            ScalarBC::TE_Neumann,
            8,
            8);

        std::ostringstream name;
        name << "te_coax_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_te
             << "_sv.vtk";
        write_mode(sys_te, te, i, name.str());
        if (exported_te == 1)
            write_mode(sys_te, te, i, "te_coax_sv.vtk");
    }

    int exported_tm = 0;
    for (int i = 0; i < (int)tm.w.size() && exported_tm < export_modes; ++i)
    {
        if (tm.w[(size_t)i] < 0.0)
            continue;

        ++exported_tm;
        const auto id = match_coax_mode_by_mass_correlation(
            mesh,
            r1,
            r2,
            sys_tm,
            tm.Zcol,
            i,
            ScalarBC::TM_Dirichlet,
            8,
            8);

        std::ostringstream name;
        name << "tm_coax_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_tm
             << "_sv.vtk";
        write_mode(sys_tm, tm, i, name.str());
        if (exported_tm == 1)
            write_mode(sys_tm, tm, i, "tm_coax_sv.vtk");
    }

    perf.post_ms += stage.elapsed_ms();
    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("helm10_coax", perf);

    return 0;
}
