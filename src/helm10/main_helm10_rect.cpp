/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_helm10_rect.cpp                                   */
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
#include "core/mesh2d_rect.hpp"
#include "core/mode_match_rect.hpp"
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

    const auto out_dir = output_paths::ensure_case_dir("2d/2.1_scalar/rect");
    std::cout << "Output dir: " << out_dir << "\n";

    const Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    std::cout << "Rect mesh: nodes=" << mesh.nodes.size() << " tris=" << mesh.tris.size() << "\n";
    std::cout << "a=" << a << " b=" << b << " nx=" << nx << " ny=" << ny << "\n\n";

    // TE (Neumann) block in the scalar formulation.
    std::cout << "[TE scalar (Neumann)]\n";
    const auto sys_te = build_helm10_scalar_system(mesh, ScalarBC::TE_Neumann);
    const auto te_res = generalized_eigs_sym_vec(sys_te.S, sys_te.T);
    helm10_post::print_positive_kc(te_res.w, 12, true);

    std::cout << "\nTabela 1 (TE): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,n)   kc_ana      kc_fem      err(%)    |rho|\n";
    int shown = 0;
    for (int i = 0; i < (int)te_res.w.size() && shown < 8; ++i)
    {
        if (te_res.w[(size_t)i] < 1e-9)
            continue; // remove constant mode

        const double kc_fem = std::sqrt(te_res.w[(size_t)i]);
        const auto id = match_rect_mode_by_mass_correlation(
            mesh,
            a,
            b,
            sys_te,
            te_res.Zcol,
            i,
            ScalarBC::TE_Neumann,
            8,
            8);
        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;

        std::cout << std::setw(1) << (shown + 1) << "  ("
                  << id.m << "," << id.n << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho
                  << "\n";
        ++shown;
    }

    // TM (Dirichlet) block in the scalar formulation.
    std::cout << "\n[TM scalar (Dirichlet)]\n";
    const auto sys_tm = build_helm10_scalar_system(mesh, ScalarBC::TM_Dirichlet);
    const auto tm_res = generalized_eigs_sym_vec(sys_tm.S, sys_tm.T);
    helm10_post::print_positive_kc(tm_res.w, 12, false);

    std::cout << "\nTabela 1 (TM): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,n)   kc_ana      kc_fem      err(%)    |rho|\n";
    shown = 0;
    for (int i = 0; i < (int)tm_res.w.size() && shown < 8; ++i)
    {
        if (tm_res.w[(size_t)i] < 0.0)
            continue;

        const double kc_fem = std::sqrt(tm_res.w[(size_t)i]);
        const auto id = match_rect_mode_by_mass_correlation(
            mesh,
            a,
            b,
            sys_tm,
            tm_res.Zcol,
            i,
            ScalarBC::TM_Dirichlet,
            8,
            8);
        const double err = 100.0 * std::abs(kc_fem - id.kc_ana) / id.kc_ana;

        std::cout << std::setw(1) << (shown + 1) << "  ("
                  << id.m << "," << id.n << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << id.kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho
                  << "\n";
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
        std::cout << "Saved: " << vtk_name << " (phi + Ft)\n";
    };

    int exported_te = 0;
    for (int i = 0; i < (int)te_res.w.size() && exported_te < export_modes; ++i)
    {
        if (te_res.w[(size_t)i] < 1e-9)
            continue;

        ++exported_te;
        const auto id = match_rect_mode_by_mass_correlation(
            mesh,
            a,
            b,
            sys_te,
            te_res.Zcol,
            i,
            ScalarBC::TE_Neumann,
            8,
            8);

        std::ostringstream name;
        name << "te_rect_m" << id.m
             << "_n" << id.n
             << "_rank" << std::setw(2) << std::setfill('0') << exported_te
             << "_sv.vtk";
        write_mode(sys_te, te_res, i, name.str());

        if (exported_te == 1)
        {
            write_mode(sys_te, te_res, i, "te10_rect_sv.vtk");
        }
    }

    int exported_tm = 0;
    for (int i = 0; i < (int)tm_res.w.size() && exported_tm < export_modes; ++i)
    {
        if (tm_res.w[(size_t)i] < 0.0)
            continue;

        ++exported_tm;
        const auto id = match_rect_mode_by_mass_correlation(
            mesh,
            a,
            b,
            sys_tm,
            tm_res.Zcol,
            i,
            ScalarBC::TM_Dirichlet,
            8,
            8);

        std::ostringstream name;
        name << "tm_rect_m" << id.m
             << "_n" << id.n
             << "_rank" << std::setw(2) << std::setfill('0') << exported_tm
             << "_sv.vtk";
        write_mode(sys_tm, tm_res, i, name.str());

        if (exported_tm == 1)
        {
            write_mode(sys_tm, tm_res, i, "tm11_rect_sv.vtk");
        }
    }
    return 0;
}
