/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_helm10_circle.cpp                                 */
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
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "core/helm10_scalar_system.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_circle.hpp"
#include "core/mode_match_circle.hpp"
#include "core/output_paths.hpp"
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
    const double r = 1.0;
    int nr = 10;
    int nt = 48;
    int export_modes = 8;
    if (argc >= 3)
    {
        nr = std::atoi(argv[1]);
        nt = std::atoi(argv[2]);
    }
    if (argc >= 4)
    {
        export_modes = std::max(1, std::atoi(argv[3]));
    }

    const auto out_dir = output_paths::ensure_case_dir("2d/2.1_scalar/circle");
    std::cout << "Output dir: " << out_dir << "\n";

    const Mesh2D mesh = make_circle_mesh(r, nr, nt);
    std::cout << "Circle mesh: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " R=" << r << " nr=" << nr << " nt=" << nt << "\n\n";

    // Modos TE vindos do problema escalar com Neumann.
    std::cout << "[TE scalar (Neumann) - circle]\n";
    const auto sys_te = build_helm10_scalar_system(mesh, ScalarBC::TE_Neumann);
    const auto te = generalized_eigs_sym_vec(sys_te.S, sys_te.T);
    helm10_post::print_positive_kc(te.w, 10, true);

    std::cout << "\nTabela 2 (TE): FEM vs Analitico (Jm' roots)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";
    int shown = 0;
    for (int i = 0; i < (int)te.w.size() && shown < 8; ++i)
    {
        if (te.w[(size_t)i] < 1e-9)
            continue;

        const double kc_fem = std::sqrt(te.w[(size_t)i]);
        const auto id = match_circle_mode_by_mass_correlation(
            mesh,
            r,
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

    // Modos TM vindos do problema escalar com Dirichlet.
    std::cout << "\n[TM scalar (Dirichlet) - circle]\n";
    const auto sys_tm = build_helm10_scalar_system(mesh, ScalarBC::TM_Dirichlet);
    const auto tm = generalized_eigs_sym_vec(sys_tm.S, sys_tm.T);
    helm10_post::print_positive_kc(tm.w, 10, false);

    std::cout << "\nTabela 2 (TM): FEM vs Analitico (Jm roots)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      err(%)    |rho|\n";
    shown = 0;
    for (int i = 0; i < (int)tm.w.size() && shown < 8; ++i)
    {
        if (tm.w[(size_t)i] < 0.0)
            continue;

        const double kc_fem = std::sqrt(tm.w[(size_t)i]);
        const auto id = match_circle_mode_by_mass_correlation(
            mesh,
            r,
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

    int exported_te = 0;
    for (int i = 0; i < (int)te.w.size() && exported_te < export_modes; ++i)
    {
        if (te.w[(size_t)i] < 1e-9)
            continue;

        ++exported_te;
        const auto id = match_circle_mode_by_mass_correlation(
            mesh,
            r,
            sys_te,
            te.Zcol,
            i,
            ScalarBC::TE_Neumann,
            8,
            8);

        std::ostringstream name;
        name << "te_circle_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_te
             << "_sv.vtk";
        write_mode(sys_te, te, i, name.str());
        if (exported_te == 1)
            write_mode(sys_te, te, i, "te_circle_sv.vtk");
    }

    int exported_tm = 0;
    for (int i = 0; i < (int)tm.w.size() && exported_tm < export_modes; ++i)
    {
        if (tm.w[(size_t)i] < 0.0)
            continue;

        ++exported_tm;
        const auto id = match_circle_mode_by_mass_correlation(
            mesh,
            r,
            sys_tm,
            tm.Zcol,
            i,
            ScalarBC::TM_Dirichlet,
            8,
            8);

        std::ostringstream name;
        name << "tm_circle_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_tm
             << "_sv.vtk";
        write_mode(sys_tm, tm, i, name.str());
        if (exported_tm == 1)
            write_mode(sys_tm, tm, i, "tm_circle_sv.vtk");
    }

    return 0;
}
