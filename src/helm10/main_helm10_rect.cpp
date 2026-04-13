/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_helm10_rect.cpp                                   */
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
#include "core/mesh2d_rect.hpp"
#include "core/mode_match_rect.hpp"
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
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct RectRunConfig
{
    double ar_m = 1.0;
    int nx = 10;
    int ny = 20;
    int mode_limit = 10;
    bool used_positional_cli = false;
    bool used_legacy_cli = false;
};

struct ModeExportRecord
{
    std::string family;
    std::string longitudinal_label;
    std::string field_status;
    double medium_k = 0.0;
    double beta = 0.0;
    double tm_ztm = 0.0;
    int below_cutoff = 0;
    int positive_rank = 0;
    int eig_index = 0;
    int m = 0;
    int n = 0;
    double kc_fem = 0.0;
    double kc_ana = 0.0;
    double kc_ar_fem = 0.0;
    double kc_ar_ana = 0.0;
    double err_percent = 0.0;
    double rho_abs = 0.0;
};

/******************************************************************************/
/* FUNCAO: token_is_integer_like                                              */
/* DESCRICAO: Detecta argumento numerico inteiro simples para manter a leitura*/
/* do CLI legado sem criar ambiguidade com flags nomeadas.                    */
/* ENTRADA: token: const std::string &.                                       */
/* SAIDA: bool.                                                               */
/******************************************************************************/
bool token_is_integer_like(const std::string &token)
{
    if (token.empty())
        return false;

    size_t start = 0;
    if (token[0] == '+' || token[0] == '-')
        start = 1;
    if (start >= token.size())
        return false;

    for (size_t i = start; i < token.size(); ++i)
    {
        if (!std::isdigit(static_cast<unsigned char>(token[i])))
            return false;
    }
    return true;
}

/******************************************************************************/
/* FUNCAO: parse_positive_double                                              */
/* DESCRICAO: Converte argumento textual para real positivo, com mensagem     */
/* clara em caso de erro.                                                     */
/* ENTRADA: token: const std::string &; name: const std::string &.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
double parse_positive_double(const std::string &token, const std::string &name)
{
    const double value = std::stod(token);
    if (value <= 0.0)
        throw std::runtime_error(name + " deve ser > 0");
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_positive_int                                                 */
/* DESCRICAO: Converte argumento textual para inteiro positivo.               */
/* ENTRADA: token: const std::string &; name: const std::string &.            */
/* SAIDA: int.                                                                */
/******************************************************************************/
int parse_positive_int(const std::string &token, const std::string &name)
{
    const int value = std::stoi(token);
    if (value <= 0)
        throw std::runtime_error(name + " deve ser > 0");
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_nonnegative_int                                              */
/* DESCRICAO: Converte argumento textual para inteiro nao negativo, util para */
/* controlar quantidade de modos exportados.                                  */
/* ENTRADA: token: const std::string &; name: const std::string &.            */
/* SAIDA: int.                                                                */
/******************************************************************************/
int parse_nonnegative_int(const std::string &token, const std::string &name)
{
    const int value = std::stoi(token);
    if (value < 0)
        throw std::runtime_error(name + " deve ser >= 0");
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_rect_run_config                                              */
/* DESCRICAO: Interpreta argumentos posicionais do caso retangular, aceitando */
/* a forma nova com `ar` e a forma legada `[nx ny [nmodos]]`.                 */
/* ENTRADA: cli: const helm10::ScalarCliOptions &.                            */
/* SAIDA: RectRunConfig.                                                      */
/******************************************************************************/
RectRunConfig parse_rect_run_config(const helm10::ScalarCliOptions &cli)
{
    RectRunConfig cfg;
    const auto &p = cli.positionals;

    if (p.empty())
        return cfg;

    cfg.used_positional_cli = true;

    if (p.size() > 4)
    {
        throw std::runtime_error(
            "argumentos posicionais em excesso; use [ar_m [nx [ny [nmodos]]]] "
            "ou o legado [nx ny [nmodos]]");
    }

    const bool legacy_cli =
        p.size() >= 2 &&
        p.size() <= 3 &&
        token_is_integer_like(p[0]) &&
        token_is_integer_like(p[1]);

    if (legacy_cli)
    {
        cfg.nx = parse_positive_int(p[0], "nx");
        cfg.ny = parse_positive_int(p[1], "ny");
        if (p.size() == 3)
            cfg.mode_limit = parse_nonnegative_int(p[2], "nmodos");
        cfg.used_legacy_cli = true;
        return cfg;
    }

    cfg.ar_m = parse_positive_double(p[0], "ar_m");
    if (p.size() >= 2)
        cfg.nx = parse_positive_int(p[1], "nx");
    if (p.size() >= 3)
        cfg.ny = parse_positive_int(p[2], "ny");
    else if (p.size() == 2)
        cfg.ny = cfg.nx;
    if (p.size() == 4)
        cfg.mode_limit = parse_nonnegative_int(p[3], "nmodos");

    return cfg;
}

/******************************************************************************/
/* FUNCAO: mode_label                                                         */
/* DESCRICAO: Gera rotulo textual simples do modo para CSVs e mensagens.      */
/* ENTRADA: family: const std::string &; m: int; n: int.                      */
/* SAIDA: std::string.                                                        */
/******************************************************************************/
std::string mode_label(const std::string &family, int m, int n)
{
    std::ostringstream oss;
    oss << family << "_m" << m << "_n" << n;
    return oss.str();
}

/******************************************************************************/
/* FUNCAO: material_adjusted_kc                                               */
/* DESCRICAO: Ajusta o valor analitico puramente geometrico de kc para o meio */
/* homogeneo usado na montagem FEM com eps_r e mu_r uniformes.                */
/* ENTRADA: kc_geometry: double; eps_r: double; mu_r: double.                 */
/* SAIDA: double.                                                             */
/******************************************************************************/
double material_adjusted_kc(
    double kc_geometry,
    double eps_r,
    double mu_r)
{
    return kc_geometry / std::sqrt(eps_r * mu_r);
}

/******************************************************************************/
/* FUNCAO: max_exported_kc_from_spectrum                                      */
/* DESCRICAO: Varre o espectro ordenado e retorna o maior kc dentre os modos  */
/* positivos que serao realmente exportados. Isso permite escolher uma        */
/* frequencia automatica que deixe todos esses modos acima do corte.          */
/* ENTRADA: lambdas: const std::vector<double> &; export_modes: int;          */
/* lambda_tol: double.                                                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
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

/******************************************************************************/
/* FUNCAO: select_reconstruction_medium                                       */
/* DESCRICAO: Escolhe o meio homogeneo usado na reconstrucao. Se o usuario    */
/* informar `--freq-hz`, respeitamos o valor. Caso contrario, tentamos usar o */
/* maior kc TM exportado como modo de referencia, impondo `Ztm = 1`. Se isso  */
/* nao bastar para propagar todos os modos pedidos, elevamos a frequencia ate */
/* cobrir o maior kc exportado.                                               */
/* ENTRADA: cli: const helm10::ScalarCliOptions &; max_exported_kc_all: double;*/
/* max_exported_tm_kc: double.                                                */
/* SAIDA: SelectedMediumInfo.                                                 */
/******************************************************************************/
SelectedMediumInfo select_reconstruction_medium(
    const helm10::ScalarCliOptions &cli,
    double max_exported_kc_all,
    double max_exported_tm_kc)
{
    using helm10::field_reconstruction::make_homogeneous_medium_above_kc;
    using helm10::field_reconstruction::make_homogeneous_medium_for_tm_reference_ztm;
    using helm10::field_reconstruction::make_homogeneous_medium_from_k;
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
            // Se a politica Ztm=1 nao for viavel para o meio escolhido,
            // caimos para a regra generica de frequencia acima do maior kc.
        }
    }

    info.medium = make_homogeneous_medium_above_kc(
        max_exported_kc_all,
        cli.eps_r,
        cli.mu_r);
    info.frequency_policy = "automatica_acima_dos_modos_exportados";
    return info;
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
    helm10::ScalarCliOptions cli;
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso preferencial: ./helm10_rect [ar_m [nx [ny [nmodos]]]]"
           << " [--backend closed-form|gauss|efgmi]"
           << " [--freq-hz F] [--eps-r E] [--mu-r M]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Uso legado ainda aceito: ./helm10_rect [nx ny [nmodos]]"
           << " [--backend closed-form|gauss|efgmi]"
           << " [--freq-hz F] [--eps-r E] [--mu-r M]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--ar-m AR] [--nx NX] [--ny NY] [--nmodos M]"
           << " (nao misture com os posicionais principais)\n";
        os << "Compatibilidade: os posicionais principais continuam aceitos, mas estao deprecated;"
           << " prefira os aliases nomeados acima.\n";
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

    RectRunConfig run;
    const bool has_named_rect_args =
        cli.ar_m_was_provided ||
        cli.nx_was_provided ||
        cli.ny_was_provided ||
        cli.nmodos_was_provided;

    if (cli.nr_was_provided || cli.nt_was_provided)
    {
        std::cerr << "Erro: helm10_rect nao aceita --nr/--nt; use --ar-m/--nx/--ny.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_rect_args && !cli.positionals.empty())
    {
        std::cerr << "Erro: nao misture aliases nomeados principais com argumentos posicionais principais em helm10_rect; escolha apenas um estilo de chamada.\n";
        print_usage(std::cerr);
        return 2;
    }

    try
    {
        if (has_named_rect_args)
        {
            if (cli.ar_m_was_provided)
                run.ar_m = cli.ar_m;
            if (cli.nx_was_provided)
                run.nx = cli.nx;
            if (cli.ny_was_provided)
                run.ny = cli.ny;
            if (cli.nmodos_was_provided)
                run.mode_limit = cli.nmodos;
        }
        else
        {
            run = parse_rect_run_config(cli);
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos posicionais de helm10_rect: "
                  << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const double a = run.ar_m;
    const double b = 0.5 * run.ar_m;
    const int nx = run.nx;
    const int ny = run.ny;
    const int export_modes = run.mode_limit;
    const int match_span = std::max(8, export_modes);
    const auto case_dir = output_paths::ensure_case_dir("helm10/rect");
    const auto vtk_dir = output_paths::ensure_case_dir("helm10/rect/vtk");
    const auto csv_dir = output_paths::ensure_case_dir("helm10/rect/csv");
    const auto linop_dir = output_paths::ensure_case_dir("helm10/rect/linop");
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
    if (run.used_legacy_cli)
    {
        std::cerr << "Aviso: a assinatura posicional legada [nx ny [nmodos]] de helm10_rect continua aceita por compatibilidade, mas esta deprecated; prefira --ar-m/--nx/--ny/--nmodos.\n";
    }
    else if (run.used_positional_cli)
    {
        std::cerr << "Aviso: os argumentos posicionais principais de helm10_rect continuam aceitos por compatibilidade, mas estao deprecated; prefira --ar-m/--nx/--ny/--nmodos.\n";
    }

    const Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    const std::vector<double> eps_r_tri(mesh.tris.size(), cli.eps_r);
    const std::vector<double> mu_r_tri(mesh.tris.size(), cli.mu_r);
    std::cout << "Rect mesh: nodes=" << mesh.nodes.size() << " tris=" << mesh.tris.size() << "\n";
    std::cout << "ar=" << a << " m, b=" << b << " m, nx=" << nx << " ny=" << ny
              << " nmodos=" << export_modes << "\n\n";

    const std::string modes_csv_path = output_paths::file_in(csv_dir, "helm10_rect_modes.csv");
    std::ofstream modes_csv(modes_csv_path);
    if (!modes_csv)
    {
        std::cerr << "Erro ao abrir CSV de modos: " << modes_csv_path << "\n";
        return 3;
    }

    modes_csv << std::setprecision(16);
    modes_csv << "family,longitudinal_label,mode_label,positive_rank,eig_index,m,n,ar_m,b_m,"
              << "freq_hz,eps_r,mu_r,k_medium,kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,"
              << "beta,ztm,below_cutoff,error_percent,rho_abs,field_status,fields_csv_file\n";

    if (cli.debug_local_blocks)
        helm10_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    // TE (Neumann) block in the scalar formulation.
    std::cout << "[TE scalar (Neumann)]\n";
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
    const auto te_res = generalized_eigs_sym_vec(sys_te.S, sys_te.T);
    te_solve_ms = stage.elapsed_ms();
    perf.solve_ms += te_solve_ms;
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(te_res.w, 1e-9);
    helm10_post::print_positive_kc(te_res.w, 12, true);

    std::cout << "\nTabela 1 (TE): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,n)   kc_ana      kc_fem      kcAr_fem    err(%)    |rho|\n";
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
            match_span,
            match_span);
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);
        const double kc_ar_fem = kc_fem * a;
        std::cout << std::setw(1) << (shown + 1) << "  ("
                  << id.m << "," << id.n << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(10) << kc_ar_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho
                  << "\n";
        ++shown;
    }

    // TM (Dirichlet) block in the scalar formulation.
    std::cout << "\n[TM scalar (Dirichlet)]\n";
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
    const auto tm_res = generalized_eigs_sym_vec(sys_tm.S, sys_tm.T);
    tm_solve_ms = stage.elapsed_ms();
    perf.solve_ms += tm_solve_ms;
    if (cli.debug_candidates)
        helm10_debug::print_positive_kc_candidates_debug(tm_res.w, 0.0);
    helm10_post::print_positive_kc(tm_res.w, 12, false);

    if (!spectral_csv::write_symmetric_problem_exports(
            linop_dir,
            "helm10_rect_te",
            sys_te.S,
            sys_te.T,
            te_res) ||
        !spectral_csv::write_symmetric_problem_exports(
            linop_dir,
            "helm10_rect_tm",
            sys_tm.S,
            sys_tm.T,
            tm_res))
    {
        std::cerr << "Erro ao escrever artefatos espectrais em " << linop_dir << "\n";
        return 4;
    }
    std::cout << "Saved: " << (linop_dir / "helm10_rect_te_S_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_te_T_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_te_eigenvalues.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_te_eigenvectors.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_tm_S_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_tm_T_crs.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_tm_eigenvalues.csv") << "\n";
    std::cout << "Saved: " << (linop_dir / "helm10_rect_tm_eigenvectors.csv") << "\n";

    const double max_exported_kc_te =
        max_exported_kc_from_spectrum(te_res.w, export_modes, 1e-9);
    const double max_exported_kc_tm =
        max_exported_kc_from_spectrum(tm_res.w, export_modes, 0.0);
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

    std::cout << "\nTabela 1 (TM): FEM vs Analitico (match por correlacao com T)\n";
    std::cout << "i  (m,n)   kc_ana      kc_fem      kcAr_fem    err(%)    |rho|\n";
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
            match_span,
            match_span);
        const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
        const double err = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);
        const double kc_ar_fem = kc_fem * a;
        std::cout << std::setw(1) << (shown + 1) << "  ("
                  << id.m << "," << id.n << ")  "
                  << std::setw(9) << std::fixed << std::setprecision(6) << kc_ana << "  "
                  << std::setw(9) << kc_fem << "  "
                  << std::setw(10) << kc_ar_fem << "  "
                  << std::setw(7) << std::setprecision(3) << err << "  "
                  << std::setw(6) << std::setprecision(4) << id.rho
                  << "\n";
        ++shown;
    }

    auto write_mode_vtk = [&](const auto &sys,
                              const auto &res,
                              ScalarBC bc,
                              int eig_idx,
                              const std::string &vtk_name) {
        const auto psi = helm10_post::extract_mode_nodal_from_Z(
            mesh, sys, res.Zcol, eig_idx, true);
        const double kc_fem = std::sqrt(res.w[(size_t)eig_idx]);
        const auto field_result = helm10::field_reconstruction::reconstruct_transverse_fields(
            mesh,
            psi,
            (bc == ScalarBC::TE_Neumann)
                ? helm10::field_reconstruction::LongitudinalScalarKind::TE_Hz
                : helm10::field_reconstruction::LongitudinalScalarKind::TM_Ez,
            kc_fem,
            medium,
            cli.backend,
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

    auto export_branch = [&](const std::string &family,
                             ScalarBC bc,
                             const auto &sys,
                             const auto &res,
                             double lambda_tol,
                             const char *legacy_vtk_name) {
        const std::string family_lower =
            (bc == ScalarBC::TE_Neumann) ? "te" : "tm";
        int exported = 0;
        for (int i = 0; i < (int)res.w.size() && exported < export_modes; ++i)
        {
            if (res.w[(size_t)i] < lambda_tol)
                continue;

            ++exported;
            const double kc_fem = std::sqrt(res.w[(size_t)i]);
            const auto id = match_rect_mode_by_mass_correlation(
                mesh,
                a,
                b,
                sys,
                res.Zcol,
                i,
                bc,
                match_span,
                match_span);
            const double kc_ana = material_adjusted_kc(id.kc_ana, cli.eps_r, cli.mu_r);
            const auto psi_raw = helm10_post::extract_mode_nodal_from_Z(
                mesh, sys, res.Zcol, i, false);
            const auto raw_field = helm10::field_reconstruction::reconstruct_transverse_fields(
                mesh,
                psi_raw,
                (bc == ScalarBC::TE_Neumann)
                    ? helm10::field_reconstruction::LongitudinalScalarKind::TE_Hz
                    : helm10::field_reconstruction::LongitudinalScalarKind::TM_Ez,
                kc_fem,
                medium,
                cli.backend,
                false);

            ModeExportRecord rec;
            rec.family = family;
            rec.longitudinal_label = raw_field.longitudinal_label;
            rec.field_status = raw_field.status_label;
            rec.medium_k = raw_field.k;
            rec.beta = raw_field.beta;
            rec.tm_ztm = raw_field.ztm;
            rec.below_cutoff = raw_field.below_cutoff ? 1 : 0;
            rec.positive_rank = exported;
            rec.eig_index = i + 1;
            rec.m = id.m;
            rec.n = id.n;
            rec.kc_fem = kc_fem;
            rec.kc_ana = kc_ana;
            rec.kc_ar_fem = kc_fem * a;
            rec.kc_ar_ana = kc_ana * a;
            rec.err_percent = error_metrics::absolute_relative_error_percent(kc_ana, kc_fem);
            rec.rho_abs = id.rho;

            if (raw_field.below_cutoff)
            {
                std::cout << "Aviso: " << family
                          << " " << mode_label(family, rec.m, rec.n)
                          << " abaixo do corte. "
                          << raw_field.status_message << "\n";
            }

            const std::string label = mode_label(rec.family, rec.m, rec.n);
            std::ostringstream fields_name;
            fields_name << "helm10_rect_fields_" << label
                        << "_rank" << std::setw(2) << std::setfill('0') << rec.positive_rank
                        << ".csv";
            const std::string fields_csv_filename = fields_name.str();
            std::ofstream mode_fields_csv(output_paths::file_in(csv_dir, fields_csv_filename));
            if (!mode_fields_csv)
            {
                throw std::runtime_error("Erro ao abrir CSV de campos por modo: " + fields_csv_filename);
            }
            mode_fields_csv << std::setprecision(16);
            mode_fields_csv << "node_id,x_m,y_m,psi,dpsi_dx,dpsi_dy,Ex,Ey\n";
            modes_csv << rec.family << ","
                      << rec.longitudinal_label << ","
                      << label << ","
                      << rec.positive_rank << ","
                      << rec.eig_index << ","
                      << rec.m << ","
                      << rec.n << ","
                      << a << ","
                      << b << ","
                      << medium.frequency_hz << ","
                      << medium.eps_r << ","
                      << medium.mu_r << ","
                      << rec.medium_k << ","
                      << rec.kc_fem << ","
                      << rec.kc_ana << ","
                      << rec.kc_ar_fem << ","
                      << rec.kc_ar_ana << ","
                      << rec.beta << ","
                      << rec.tm_ztm << ","
                      << rec.below_cutoff << ","
                      << rec.err_percent << ","
                      << rec.rho_abs << ","
                      << rec.field_status << ","
                      << fields_csv_filename << "\n";

            for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
            {
                const auto &node = mesh.nodes[node_id];
                mode_fields_csv << node_id << ","
                                << node.x << ","
                                << node.y << ","
                                << raw_field.psi[node_id] << ","
                                << raw_field.dpsi_dx[node_id] << ","
                                << raw_field.dpsi_dy[node_id] << ","
                                << raw_field.ex[node_id] << ","
                                << raw_field.ey[node_id] << "\n";
            }
            mode_fields_csv.close();

            std::ostringstream name;
            name << family_lower << "_rect_m" << rec.m
                 << "_n" << rec.n
                 << "_rank" << std::setw(2) << std::setfill('0') << rec.positive_rank
                 << "_sv.vtk";
            write_mode_vtk(sys, res, bc, i, name.str());

            if (rec.positive_rank == 1)
            {
                write_mode_vtk(sys, res, bc, i, legacy_vtk_name);
            }
        }
    };

    stage.reset();
    export_branch("TE", ScalarBC::TE_Neumann, sys_te, te_res, 1e-9, "te10_rect_sv.vtk");
    export_branch("TM", ScalarBC::TM_Dirichlet, sys_tm, tm_res, 0.0, "tm11_rect_sv.vtk");
    perf.post_ms += stage.elapsed_ms();
    std::cout << "Saved: " << modes_csv_path << "\n";
    modes_csv.close();
    perf.total_ms = total_watch.elapsed_ms();
    run_timing_csv::Record timing_record;
    timing_record.case_label = "helm10_rect";
    timing_record.geometry = "rect";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.legacy_cli_used = run.used_legacy_cli ? 1 : 0;
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
    timing_record.ar_m = std::to_string(a);
    timing_record.b_m = std::to_string(b);
    timing_record.nx = std::to_string(nx);
    timing_record.ny = std::to_string(ny);
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
    timing::print_breakdown("helm10_rect", perf);

    return 0;
}
