/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_helm10_circle.cpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
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

#include "article/tp3485_systems.hpp"
#include "core/error_metrics.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_circle.hpp"
#include "core/mode_match_circle.hpp"
#include "core/execution_log.hpp"
#include "core/output_paths.hpp"
#include "core/run_timing_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "helm10/field_reconstruction.hpp"
#include "helm10/scalar_cli_options.hpp"
#include "helm10/scalar_debug.hpp"
#include "helm10/scalar_mode_post.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
double material_adjusted_kc(double kc_geometry, double eps_r, double mu_r)
{
    return kc_geometry / std::sqrt(eps_r * mu_r);
}

double max_exported_kc_from_spectrum(
    const std::vector<double> &lambdas,
    int export_modes,
    double lambda_tol)
{
    double max_kc = 0.0;
    int exported = 0;
    for (double lambda : lambdas)
    {
        if (lambda < lambda_tol)
            continue;

        max_kc = std::max(max_kc, std::sqrt(lambda));
        ++exported;
        if (exported >= export_modes)
            break;
    }
    return max_kc;
}

struct SelectedMediumInfo
{
    helm10::field_reconstruction::HomogeneousMedium medium;
    std::string frequency_policy;
    double reference_tm_kc = 0.0;
    double reference_tm_ztm_actual = 0.0;
};

double tm_ztm_from_medium(
    const helm10::field_reconstruction::HomogeneousMedium &medium,
    double kc)
{
    const double beta_squared = medium.k * medium.k - kc * kc;
    if (beta_squared <= 0.0)
        return 0.0;
    return std::sqrt(beta_squared) / (medium.omega_rad_s * medium.epsilon);
}

SelectedMediumInfo select_reconstruction_medium(
    const helm10::ScalarCliOptions &cli,
    double max_exported_kc_all,
    double max_exported_tm_kc)
{
    using helm10::field_reconstruction::make_homogeneous_medium_for_tm_reference_ztm;
    using helm10::field_reconstruction::make_homogeneous_medium_from_k;
    using helm10::field_reconstruction::make_homogeneous_medium_above_kc;
    using helm10::field_reconstruction::make_homogeneous_medium_from_relative;

    SelectedMediumInfo info;
    if (cli.frequency_was_provided)
    {
        info.frequency_policy = "explicita_do_usuario";
        info.medium = make_homogeneous_medium_from_relative(
            cli.frequency_hz,
            cli.eps_r,
            cli.mu_r);
        info.reference_tm_kc = max_exported_tm_kc;
        if (max_exported_tm_kc > 0.0)
            info.reference_tm_ztm_actual = tm_ztm_from_medium(info.medium, max_exported_tm_kc);
        return info;
    }

    if (max_exported_tm_kc > 0.0)
    {
        info.reference_tm_kc = max_exported_tm_kc;
        try
        {
            const auto tm_reference_medium =
                make_homogeneous_medium_for_tm_reference_ztm(
                    max_exported_tm_kc,
                    1.0,
                    cli.eps_r,
                    cli.mu_r);
            const double propagation_floor =
                std::max(1.0e-9, max_exported_kc_all) * 1.000001;

            if (tm_reference_medium.k >= propagation_floor)
            {
                info.medium = tm_reference_medium;
                info.frequency_policy = "automatica_tm_referencia_ztm1";
            }
            else
            {
                info.medium = make_homogeneous_medium_from_k(
                    propagation_floor,
                    cli.eps_r,
                    cli.mu_r);
                info.frequency_policy =
                    "automatica_tm_referencia_ztm1_ajustada_para_propagar_todos";
            }

            info.reference_tm_ztm_actual =
                tm_ztm_from_medium(info.medium, max_exported_tm_kc);
            return info;
        }
        catch (const std::exception &)
        {
        }
    }

    info.medium = make_homogeneous_medium_above_kc(
        max_exported_kc_all,
        cli.eps_r,
        cli.mu_r);
    info.frequency_policy = "automatica_acima_dos_modos_exportados";
    return info;
}

std::string circle_mode_label(const std::string &family, int m, int p)
{
    std::ostringstream oss;
    oss << family << "_m" << m << "_p" << p;
    return oss.str();
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
    const double r = 1.0;
    int nr = 10;
    int nt = 48;
    int export_modes = 20;
    helm10::ScalarCliOptions cli;
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso: ./helm10_circle [nr nt [nmodos]] [--backend closed-form|gauss]"
           << " [--freq-hz F] [--eps-r E] [--mu-r M]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--nr NR] [--nt NT] [--nmodos M]"
           << " (nao misture com os posicionais principais)\n";
    };
    if (helm10::scalar_cli_requests_help(argc, argv))
    {
        print_usage(std::cout);
        return 0;
    }
    try
    {
        cli = helm10::parse_scalar_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const bool has_named_polar_args =
        cli.nr_was_provided ||
        cli.nt_was_provided ||
        cli.nmodos_was_provided;

    if (cli.ar_m_was_provided || cli.nx_was_provided || cli.ny_was_provided)
    {
        std::cerr << "Erro: helm10_circle nao aceita --ar-m/--nx/--ny; use --nr/--nt.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_polar_args && !cli.positionals.empty())
    {
        std::cerr << "Erro: nao misture aliases nomeados principais com os argumentos posicionais de helm10_circle.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (!has_named_polar_args && !cli.positionals.empty() && cli.positionals.size() < 2)
    {
        std::cerr << "Erro: use ./helm10_circle [nr nt [nmodos]]"
                  << " [--backend closed-form|gauss]"
                  << " [--freq-hz F] [--eps-r E] [--mu-r M]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        print_usage(std::cerr);
        return 2;
    }

    try
    {
        if (has_named_polar_args)
        {
            if (cli.nr_was_provided)
                nr = cli.nr;
            if (cli.nt_was_provided)
                nt = cli.nt;
            if (cli.nmodos_was_provided)
                export_modes = cli.nmodos;
        }
        else
        {
            if (cli.positionals.size() >= 2)
            {
                nr = helm10::parse_positive_cli_int(cli.positionals[0], "nr");
                nt = helm10::parse_positive_cli_int(cli.positionals[1], "nt");
            }
            if (cli.positionals.size() >= 3)
            {
                export_modes = helm10::parse_nonnegative_cli_int(cli.positionals[2], "nmodos");
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos posicionais de helm10_circle: "
                  << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const auto case_dir = output_paths::ensure_case_dir("helm10/circle");
    const auto vtk_dir = output_paths::ensure_case_dir("helm10/circle/vtk");
    const auto csv_dir = output_paths::ensure_case_dir("helm10/circle/csv");
    const auto linop_dir = output_paths::ensure_case_dir("helm10/circle/linop");
    const std::string timing_csv_path = output_paths::file_in(case_dir, "run_timing.csv");
    execution_log::ExecutionLogScope execution_log(
        output_paths::file_in(case_dir, "run.log"));
    if (execution_log.active())
        std::cout << "Log file: " << execution_log.file_path() << "\n";
    else
        std::cerr << "Aviso: nao foi possivel abrir log de execucao em "
                  << execution_log.file_path() << ": "
                  << execution_log.error_message() << "\n";
    std::cout << "Output dir: " << case_dir << "\n";
    std::cout << "VTK dir: " << vtk_dir << "\n";
    std::cout << "CSV dir: " << csv_dir << "\n";
    std::cout << "LinOp dir: " << linop_dir << "\n";
    std::cout << "Backend escalar: " << element_assembly_backend_name(cli.backend) << "\n";

    const Mesh2D mesh = make_circle_mesh(r, nr, nt);
    const std::vector<double> eps_r_tri(mesh.tris.size(), cli.eps_r);
    const std::vector<double> mu_r_tri(mesh.tris.size(), cli.mu_r);
    std::cout << "Circle mesh: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " R=" << r << " nr=" << nr << " nt=" << nt << "\n\n";

    const std::string modes_csv_path = output_paths::file_in(csv_dir, "helm10_circle_modes.csv");
    std::ofstream modes_csv(modes_csv_path);
    if (!modes_csv)
    {
        std::cerr << "Erro ao abrir CSV de modos do caso circular.\n";
        return 3;
    }
    modes_csv << std::setprecision(16);
    modes_csv << "family,longitudinal_label,mode_label,positive_rank,eig_index,m,p,r_m,"
              << "freq_hz,eps_r,mu_r,k_medium,kc_fem,kc_ana,kc_r_fem,kc_r_ana,"
              << "beta,ztm,below_cutoff,error_percent,rho_abs,field_status,fields_csv_file\n";

    if (cli.debug_local_blocks)
        helm10_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    // Modos TE vindos do problema escalar com Neumann.
    std::cout << "[TE scalar (Neumann) - circle]\n";
    timing::Stopwatch stage;
    double te_assembly_ms = 0.0;
    double te_solve_ms = 0.0;
    double tm_assembly_ms = 0.0;
    double tm_solve_ms = 0.0;
    const auto sys_te = tp3485::build_eq43_helm10_system(
        mesh,
        ScalarBC::TE_Neumann,
        eps_r_tri,
        mu_r_tri,
        cli.backend);
    te_assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += te_assembly_ms;
    stage.reset();
    const auto te = generalized_eigs_sym_vec(sys_te.S, sys_te.T);
    te_solve_ms = stage.elapsed_ms();
    perf.solve_ms += te_solve_ms;
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(te.w, 1e-9);
    helm10_post::print_positive_kc(te.w, 10, true);

    std::cout << "\nTabela 2 (TE): FEM vs Analitico (Jm' roots)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      kcR_fem     err(%)    |rho|\n";
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
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);

        std::cout << (shown + 1) << "  (" << id.m << "," << id.p << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(10) << kc_fem * r << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho << "\n";
        ++shown;
    }

    // Modos TM vindos do problema escalar com Dirichlet.
    std::cout << "\n[TM scalar (Dirichlet) - circle]\n";
    stage.reset();
    const auto sys_tm = tp3485::build_eq43_helm10_system(
        mesh,
        ScalarBC::TM_Dirichlet,
        eps_r_tri,
        mu_r_tri,
        cli.backend);
    tm_assembly_ms = stage.elapsed_ms();
    perf.assembly_ms += tm_assembly_ms;
    stage.reset();
    const auto tm = generalized_eigs_sym_vec(sys_tm.S, sys_tm.T);
    tm_solve_ms = stage.elapsed_ms();
    perf.solve_ms += tm_solve_ms;
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(tm.w, 0.0);
    helm10_post::print_positive_kc(tm.w, 10, false);

    if (!spectral_csv::write_symmetric_problem_exports(
            linop_dir,
            "helm10_circle_te",
            sys_te.S,
            sys_te.T,
            te) ||
        !spectral_csv::write_symmetric_problem_exports(
            linop_dir,
            "helm10_circle_tm",
            sys_tm.S,
            sys_tm.T,
            tm))
    {
        std::cerr << "Erro ao escrever artefatos espectrais em " << linop_dir << "\n";
        return 4;
    }
    std::cout << "Saved: " << (linop_dir / "helm10_circle_te_S_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_te_T_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_te_eigenvalues.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_te_eigenvectors.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_tm_S_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_tm_T_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_tm_eigenvalues.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_circle_tm_eigenvectors.csv") << "\n";

    const double max_exported_kc_te =
        max_exported_kc_from_spectrum(te.w, export_modes, 1e-9);
    const double max_exported_kc_tm =
        max_exported_kc_from_spectrum(tm.w, export_modes, 0.0);
    const double max_exported_kc = std::max(max_exported_kc_te, max_exported_kc_tm);
    const auto medium_info = select_reconstruction_medium(
        cli,
        max_exported_kc,
        max_exported_kc_tm);
    const auto &medium = medium_info.medium;
    std::cout << "Reconstrucao de campos: freq=" << medium.frequency_hz
              << " Hz, eps_r=" << medium.eps_r
              << ", mu_r=" << medium.mu_r
              << ", k=" << medium.k << " rad/m"
              << ", politica=" << medium_info.frequency_policy << "\n";
    if (medium_info.reference_tm_kc > 0.0)
    {
        std::cout << "Referencia automatica TM: kc_ref=" << medium_info.reference_tm_kc
                  << ", Ztm_ref_real=" << medium_info.reference_tm_ztm_actual << "\n";
    }
    std::cout << "Aviso didatico: os campos TM exportados usam gradiente+sinal como saida principal, "
                 "sem multiplicacao por Ztm; beta e Ztm ficam apenas no modes.csv.\n";

    std::cout << "\nTabela 2 (TM): FEM vs Analitico (Jm roots)\n";
    std::cout << "i  (m,p)   kc_ana      kc_fem      kcR_fem     err(%)    |rho|\n";
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
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);

        std::cout << (shown + 1) << "  (" << id.m << "," << id.p << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(10) << kc_fem * r << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho << "\n";
        ++shown;
    }

    auto write_mode = [&](const auto &sys,
                          const auto &res,
                          ScalarBC bc,
                          int eig_idx,
                          const std::string &vtk_name) {
        const auto psi = helm10_post::extract_mode_nodal_from_Z(mesh, sys, res.Zcol, eig_idx);
        const auto field_result = helm10::field_reconstruction::reconstruct_transverse_fields(
            mesh,
            psi,
            (bc == ScalarBC::TE_Neumann)
                ? helm10::field_reconstruction::LongitudinalScalarKind::TE_Hz
                : helm10::field_reconstruction::LongitudinalScalarKind::TM_Ez,
            std::sqrt(res.w[(size_t)eig_idx]),
            medium,
            true);
        write_vtk_unstructured_tri_scalar_vector(
            output_paths::file_in(vtk_dir, vtk_name),
            mesh,
            field_result.psi,
            field_result.ex,
            field_result.ey,
            "psi",
            "Et");
        std::cout << "Saved: " << vtk_name << " (psi + Et)\n";
    };

    stage.reset();
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
        const double kc_fem = std::sqrt(te.w[(size_t)i]);
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const auto psi_raw = helm10_post::extract_mode_nodal_from_Z(mesh, sys_te, te.Zcol, i, false);
        const auto field_result = helm10::field_reconstruction::reconstruct_transverse_fields(
            mesh,
            psi_raw,
            helm10::field_reconstruction::LongitudinalScalarKind::TE_Hz,
            kc_fem,
            medium,
            false);
        const std::string label = circle_mode_label("TE", id.m, id.p);
        std::ostringstream fields_name;
        fields_name << "helm10_circle_fields_" << label
                    << "_rank" << std::setw(2) << std::setfill('0') << exported_te
                    << ".csv";
        const std::string fields_csv_filename = fields_name.str();
        std::ofstream mode_fields_csv(output_paths::file_in(csv_dir, fields_csv_filename));
        if (!mode_fields_csv)
        {
            throw std::runtime_error("Erro ao abrir CSV de campos por modo: " + fields_csv_filename);
        }
        mode_fields_csv << std::setprecision(16);
        mode_fields_csv << "node_id,x_m,y_m,psi,dpsi_dx,dpsi_dy,Ex,Ey\n";
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);

        if (field_result.below_cutoff)
        {
            std::cout << "Aviso: TE " << label << " abaixo do corte. "
                      << field_result.status_message << "\n";
        }

        modes_csv << "TE,"
                  << field_result.longitudinal_label << ","
                  << label << ","
                  << exported_te << ","
                  << (i + 1) << ","
                  << id.m << ","
                  << id.p << ","
                  << r << ","
                  << medium.frequency_hz << ","
                  << medium.eps_r << ","
                  << medium.mu_r << ","
                  << field_result.k << ","
                  << kc_fem << ","
                  << kc_ana << ","
                  << kc_fem * r << ","
                  << kc_ana * r << ","
                  << field_result.beta << ","
                  << field_result.ztm << ","
                  << (field_result.below_cutoff ? 1 : 0) << ","
                  << err << ","
                  << id.rho << ","
                  << field_result.status_label << ","
                  << fields_csv_filename << "\n";

        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            const auto &node = mesh.nodes[node_id];
            mode_fields_csv << node_id << ","
                            << node.x << ","
                            << node.y << ","
                            << field_result.psi[node_id] << ","
                            << field_result.dpsi_dx[node_id] << ","
                            << field_result.dpsi_dy[node_id] << ","
                            << field_result.ex[node_id] << ","
                            << field_result.ey[node_id] << "\n";
        }
        mode_fields_csv.close();

        std::ostringstream name;
        name << "te_circle_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_te
             << "_sv.vtk";
        write_mode(sys_te, te, ScalarBC::TE_Neumann, i, name.str());
        if (exported_te == 1)
            write_mode(sys_te, te, ScalarBC::TE_Neumann, i, "te_circle_sv.vtk");
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
        const double kc_fem = std::sqrt(tm.w[(size_t)i]);
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const auto psi_raw = helm10_post::extract_mode_nodal_from_Z(mesh, sys_tm, tm.Zcol, i, false);
        const auto field_result = helm10::field_reconstruction::reconstruct_transverse_fields(
            mesh,
            psi_raw,
            helm10::field_reconstruction::LongitudinalScalarKind::TM_Ez,
            kc_fem,
            medium,
            false);
        const std::string label = circle_mode_label("TM", id.m, id.p);
        std::ostringstream fields_name;
        fields_name << "helm10_circle_fields_" << label
                    << "_rank" << std::setw(2) << std::setfill('0') << exported_tm
                    << ".csv";
        const std::string fields_csv_filename = fields_name.str();
        std::ofstream mode_fields_csv(output_paths::file_in(csv_dir, fields_csv_filename));
        if (!mode_fields_csv)
        {
            throw std::runtime_error("Erro ao abrir CSV de campos por modo: " + fields_csv_filename);
        }
        mode_fields_csv << std::setprecision(16);
        mode_fields_csv << "node_id,x_m,y_m,psi,dpsi_dx,dpsi_dy,Ex,Ey\n";
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);

        if (field_result.below_cutoff)
        {
            std::cout << "Aviso: TM " << label << " abaixo do corte. "
                      << field_result.status_message << "\n";
        }

        modes_csv << "TM,"
                  << field_result.longitudinal_label << ","
                  << label << ","
                  << exported_tm << ","
                  << (i + 1) << ","
                  << id.m << ","
                  << id.p << ","
                  << r << ","
                  << medium.frequency_hz << ","
                  << medium.eps_r << ","
                  << medium.mu_r << ","
                  << field_result.k << ","
                  << kc_fem << ","
                  << kc_ana << ","
                  << kc_fem * r << ","
                  << kc_ana * r << ","
                  << field_result.beta << ","
                  << field_result.ztm << ","
                  << (field_result.below_cutoff ? 1 : 0) << ","
                  << err << ","
                  << id.rho << ","
                  << field_result.status_label << ","
                  << fields_csv_filename << "\n";

        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            const auto &node = mesh.nodes[node_id];
            mode_fields_csv << node_id << ","
                            << node.x << ","
                            << node.y << ","
                            << field_result.psi[node_id] << ","
                            << field_result.dpsi_dx[node_id] << ","
                            << field_result.dpsi_dy[node_id] << ","
                            << field_result.ex[node_id] << ","
                            << field_result.ey[node_id] << "\n";
        }
        mode_fields_csv.close();

        std::ostringstream name;
        name << "tm_circle_m" << id.m
             << "_p" << id.p
             << "_rank" << std::setw(2) << std::setfill('0') << exported_tm
             << "_sv.vtk";
        write_mode(sys_tm, tm, ScalarBC::TM_Dirichlet, i, name.str());
        if (exported_tm == 1)
            write_mode(sys_tm, tm, ScalarBC::TM_Dirichlet, i, "tm_circle_sv.vtk");
    }

    perf.post_ms += stage.elapsed_ms();
    std::cout << "Saved: " << modes_csv_path << "\n";
    modes_csv.close();
    perf.total_ms = total_watch.elapsed_ms();
    run_timing_csv::Record timing_record;
    timing_record.case_label = "helm10_circle";
    timing_record.geometry = "circle";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.debug_local_blocks = cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cli.debug_candidates ? 1 : 0;
    timing_record.frequency_was_provided = cli.frequency_was_provided ? 1 : 0;
    timing_record.frequency_input_hz = cli.frequency_was_provided ? cli.frequency_hz : 0.0;
    timing_record.frequency_policy = medium_info.frequency_policy;
    timing_record.frequency_used_hz = medium.frequency_hz;
    timing_record.eps_r = cli.eps_r;
    timing_record.mu_r = cli.mu_r;
    timing_record.k_medium = medium.k;
    timing_record.reference_tm_kc = medium_info.reference_tm_kc;
    timing_record.reference_tm_ztm_actual = medium_info.reference_tm_ztm_actual;
    timing_record.r_m = std::to_string(r);
    timing_record.nr = std::to_string(nr);
    timing_record.nt = std::to_string(nt);
    timing_record.nmodos = export_modes;
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.te_assembly_ms = te_assembly_ms;
    timing_record.te_solve_ms = te_solve_ms;
    timing_record.tm_assembly_ms = tm_assembly_ms;
    timing_record.tm_solve_ms = tm_solve_ms;
    timing_record.assembly_ms_total = perf.assembly_ms;
    timing_record.solve_ms_total = perf.solve_ms;
    timing_record.post_ms = perf.post_ms;
    timing_record.total_ms = perf.total_ms;
    if (run_timing_csv::write_csv(timing_csv_path, timing_record))
        std::cout << "Saved: " << timing_csv_path << "\n";
    else
        std::cerr << "Aviso: nao foi possivel salvar tempos em "
                  << timing_csv_path << "\n";
    timing::print_breakdown("helm10_circle", perf);

    return 0;
}
