/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec3/main_helmvec3_rect.cpp                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Implementacao compartilhada dos drivers publicos do HELMVEC3.   */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.4, Eq.    */
/* (126)-(127), Fig. 12-13, Tabelas 9-10.                                     */
/*****************************************************************************/
/* Observacao: Este arquivo concentra o miolo numerico comum das entradas     */
/* publicas `helmvec3_fig12_rect` e `helmvec3_fig13_rect`. A montagem global  */
/* correspondente a Eq. (136) ocorre em src/helmvec2/helmvec2_coupled_system. */
/* Comentarios priorizam didatica, rastreabilidade e validacao.               */
/*****************************************************************************/

#include "helmvec3/helmvec3_driver_entry.hpp"
#include "article/tp3485_systems.hpp"
#include "core/error_metrics.hpp"
#include "core/execution_log.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/run_timing_dispersion_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "explicit/tri2d_coupled_explicit.hpp"
#include "helmvec2/helmvec23_shared.hpp"
#include "helmvec2/coupled_cli_options.hpp"
#include "helmvec2/helmvec2_coupled_system.hpp"
#include "helmvec3_case_output.hpp"
#include "helmvec3_field_output.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct Table10Block
{
    double d_over_a = 0.0;
    std::vector<double> analytic_beta_over_k0;
    std::vector<double> helmvec3_beta_over_k0;
};

struct RatioCandidate
{
    int candidate_rank = 0;
    int eig_index = -1;
    double beta_ratio = std::numeric_limits<double>::quiet_NaN();
    double ez_ratio = 0.0;
};

struct SelectedRatioPoint
{
    double beta_ratio = std::numeric_limits<double>::quiet_NaN();
    int selected_candidate_rank = 0;
    int selected_eig_index = -1;
    double ez_ratio = 0.0;
    std::string match_status = "missing_candidate";
    helmvec3_field::SpatialArtifacts spatial;
};

struct BetaPointSolve
{
    CoupledBetaSystem sys;
    GenEigGeneralResult eig;
    std::vector<RatioCandidate> candidates;
};

/******************************************************************************/
/* FUNCAO: print_block_3x3                                                    */
/* DESCRICAO: Imprime uma matriz 3x3 em formato compacto para depuracao       */
/* didatica dos blocos locais da formulacao acoplada.                         */
/* ENTRADA: name: const char *; M: const double[3][3].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void print_block_3x3(const char *name, const double M[3][3])
{
    std::cout << "  " << name << " =\n";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << "   ";
        for (int j = 0; j < 3; ++j)
            std::cout << " " << M[i][j];
        std::cout << "\n";
    }
}

std::string label_number(double value)
{
    std::string text = std::to_string(value);
    for (char &ch : text)
    {
        if (ch == '.')
            ch = '_';
    }
    while (!text.empty() && text.back() == '0')
        text.pop_back();
    if (!text.empty() && text.back() == '_')
        text.pop_back();
    if (text.empty())
        text = "0";
    return text;
}

/******************************************************************************/
/* FUNCAO: print_first_triangle_closed_form_debug                             */
/* DESCRICAO: Imprime os blocos locais closed-form do primeiro triangulo para */
/* a formulacao de 2.2.4, mostrando tanto as Eq. (137)-(142) quanto a forma  */
/* rearranjada usada no sistema global Eq. (136).                             */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps: const std::vector<double>  */
/* &; mu: const std::vector<double> &.                                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void print_first_triangle_closed_form_debug(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps,
    const std::vector<double> &mu)
{
    if (mesh.tris.empty())
        return;

    const Tri &t = mesh.tris.front();
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    const auto blk142 = explicit_tri2d::tri2d_beta_closed_form_eq_137_142(
        tg,
        k0,
        eps.front(),
        mu.front());
    const auto blk136 = explicit_tri2d::tri2d_beta_rearranged_closed_form_eq_136(
        tg,
        k0,
        eps.front(),
        mu.front());

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form de 2.2.4\n";
    std::cout << "  k0_amostra=" << k0 << "\n";
    std::cout << "  artigo Eq. (137)-(142):\n";
    print_block_3x3("Sel_tt", blk142.Sel_tt);
    print_block_3x3("Tel_tz", blk142.Tel_tz);
    print_block_3x3("Tel_zt", blk142.Tel_zt);
    print_block_3x3("Sel_zz", blk142.Sel_zz);
    print_block_3x3("Tel_tt", blk142.Tel_tt);
    print_block_3x3("Tel_zz", blk142.Tel_zz);

    std::cout << "  rearranjo usado no codigo para Eq. (136):\n";
    print_block_3x3("P_tt", blk136.P_tt);
    print_block_3x3("P_zz", blk136.P_zz);
    print_block_3x3("Q_tt", blk136.Q_tt);
    print_block_3x3("Q_tz", blk136.Q_tz);
    print_block_3x3("Q_zt", blk136.Q_zt);
    print_block_3x3("Q_zz", blk136.Q_zz);
}

/******************************************************************************/
/* FUNCAO: beta_ratio_candidates_from_k0                                      */
/* DESCRICAO: Resolve o sistema de beta para um k0 fixo e filtra candidatos   */
/* fisicos de beta/k0 para construir ramos de dispersao da Secao 2.2.4        */
/* (Eq. 126-127).                                                             */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; k0: double; backend: ElementAssemblyBackend.        */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
BetaPointSolve solve_beta_point(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    double k0,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr,
    const std::filesystem::path *linop_dir = nullptr,
    const std::string *artifact_prefix = nullptr)
{
    timing::Stopwatch stage;
    BetaPointSolve result;
    result.sys = tp3485::build_eq136_helmvec3_system_E(mesh, k0, eps, mu, backend);
    if (perf != nullptr)
        perf->assembly_ms += stage.elapsed_ms();
    stage.reset();
    result.eig = generalized_eigs_real_vec(result.sys.P, result.sys.Q);
    if (perf != nullptr)
        perf->solve_ms += stage.elapsed_ms();
    if (linop_dir != nullptr && artifact_prefix != nullptr)
    {
        const std::string &prefix = *artifact_prefix;
        if (!spectral_csv::write_general_problem_exports_named(*linop_dir, prefix, "P", "Q", result.sys.P, result.sys.Q, result.eig) ||
            !spectral_csv::write_dense_crs_csv((*linop_dir) / (prefix + "_P_tt_crs.csv"), result.sys.P_tt) ||
            !spectral_csv::write_dense_crs_csv((*linop_dir) / (prefix + "_P_zz_crs.csv"), result.sys.P_zz) ||
            !spectral_csv::write_dense_crs_csv((*linop_dir) / (prefix + "_Q_tt_crs.csv"), result.sys.Q_tt) ||
            !spectral_csv::write_rect_like_crs_csv(
                (*linop_dir) / (prefix + "_Q_tz_crs.csv"),
                result.sys.Q_tz.nr,
                result.sys.Q_tz.nc,
                [&](int i, int j)
                { return result.sys.Q_tz(i, j); }) ||
            !spectral_csv::write_rect_like_crs_csv(
                (*linop_dir) / (prefix + "_Q_zt_crs.csv"),
                result.sys.Q_zt.nr,
                result.sys.Q_zt.nc,
                [&](int i, int j)
                { return result.sys.Q_zt(i, j); }) ||
            !spectral_csv::write_dense_crs_csv((*linop_dir) / (prefix + "_Q_zz_crs.csv"), result.sys.Q_zz))
        {
            throw std::runtime_error("Falha ao escrever artefatos espectrais de " + prefix);
        }
    }

    double eps_max = 1.0;
    for (double e : eps)
        eps_max = std::max(eps_max, e);
    const double ratio_max = std::sqrt(eps_max) + 0.25;

    struct RawRatioCandidate
    {
        int eig_index = -1;
        double beta_ratio = 0.0;
        double ez_ratio = 0.0;
    };

    std::vector<RawRatioCandidate> raw;
    const int n = result.eig.n;
    for (int i = 0; i < n; ++i)
    {
        if (!std::isfinite(result.eig.lambda_re[i]))
            continue;
        if (std::abs(result.eig.lambda_im[i]) > 1e-7)
            continue;
        if (result.eig.lambda_re[i] <= 1e-10)
            continue;

        const double beta = std::sqrt(result.eig.lambda_re[i]);
        const double ratio = beta / k0;
        if (ratio <= 1e-6 || ratio > ratio_max)
            continue;

        double et = 0.0;
        double ez = 0.0;
        for (int r = 0; r < result.sys.nt; ++r)
        {
            const double v = result.eig.VRcol[(size_t)i * (size_t)n + (size_t)r];
            et += v * v;
        }
        for (int r = 0; r < result.sys.nz; ++r)
        {
            const double v = result.eig.VRcol[(size_t)i * (size_t)n + (size_t)(result.sys.nt + r)];
            ez += v * v;
        }
        const double den = et + ez;
        const double ratio_ez = (den > 0.0) ? (ez / den) : 0.0;
        raw.push_back({i, ratio, ratio_ez});
    }
    std::sort(raw.begin(), raw.end(), [](const RawRatioCandidate &a, const RawRatioCandidate &b)
              { return a.beta_ratio < b.beta_ratio; });

    for (const RawRatioCandidate &cand : raw)
    {
        if (!result.candidates.empty() &&
            std::abs(cand.beta_ratio - result.candidates.back().beta_ratio) < 1e-8)
        {
            if (cand.ez_ratio > result.candidates.back().ez_ratio)
            {
                result.candidates.back().eig_index = cand.eig_index;
                result.candidates.back().ez_ratio = cand.ez_ratio;
            }
            continue;
        }

        RatioCandidate rep;
        rep.candidate_rank = static_cast<int>(result.candidates.size()) + 1;
        rep.eig_index = cand.eig_index;
        rep.beta_ratio = cand.beta_ratio;
        rep.ez_ratio = cand.ez_ratio;
        result.candidates.push_back(rep);
    }

    return result;
}

/******************************************************************************/
/* FUNCAO: match_ratio_to_reference                                           */
/* DESCRICAO: Associa resultados FEM a referencias analiticas usando criterio */
/* de proximidade em beta/k0 para comparar com Tabela 9 e Tabela 10.          */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; br_over_lambda: const std::vector<double> &; b:     */
/* double; ref_ratio: const std::vector<double> &; debug_candidates: bool;    */
/* backend: ElementAssemblyBackend.                                           */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<SelectedRatioPoint> match_ratio_to_reference(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    const std::vector<double> &br_over_lambda,
    double b,
    const std::vector<double> &ref_ratio,
    bool debug_candidates,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr,
    const std::filesystem::path *linop_dir = nullptr,
    const std::string &prefix_base = "",
    const helmvec3_output::CaseDirs *out_dirs = nullptr)
{
    std::vector<SelectedRatioPoint> out;
    out.reserve(br_over_lambda.size());

    for (int i = 0; i < (int)br_over_lambda.size(); ++i)
    {
        if (!std::isfinite(ref_ratio[i]))
        {
            out.push_back(SelectedRatioPoint{});
            continue;
        }

        const double s = br_over_lambda[i];
        const double k0 = 2.0 * M_PI * s / b;
        const std::string point_prefix =
            prefix_base.empty() ? std::string() : (prefix_base + "_br" + label_number(s));
        auto solve = solve_beta_point(
            mesh,
            eps,
            mu,
            k0,
            backend,
            perf,
            linop_dir,
            prefix_base.empty() ? nullptr : &point_prefix);
        if (solve.candidates.empty())
        {
            out.push_back(SelectedRatioPoint{});
            continue;
        }

        int best = 0;
        double best_err = std::abs(solve.candidates[0].beta_ratio - ref_ratio[i]);
        for (int j = 1; j < (int)solve.candidates.size(); ++j)
        {
            const double e = std::abs(solve.candidates[j].beta_ratio - ref_ratio[i]);
            if (e < best_err)
            {
                best_err = e;
                best = j;
            }
        }

        SelectedRatioPoint selected;
        selected.beta_ratio = solve.candidates[best].beta_ratio;
        selected.selected_candidate_rank = solve.candidates[best].candidate_rank;
        selected.selected_eig_index = solve.candidates[best].eig_index;
        selected.ez_ratio = solve.candidates[best].ez_ratio;
        selected.match_status = "matched";
        if (out_dirs != nullptr && !point_prefix.empty())
        {
            selected.spatial = helmvec3_field::export_mode(
                *out_dirs,
                mesh,
                solve.sys,
                solve.eig.VRcol,
                selected.selected_eig_index,
                point_prefix);
        }
        out.push_back(selected);

        if (debug_candidates)
        {
            std::cout << "  [debug] s=" << s << " cands:";
            for (const RatioCandidate &c : solve.candidates)
                std::cout << " " << c.beta_ratio;
            std::cout << "\n";
        }
    }

    return out;
}

/******************************************************************************/
/* FUNCAO: trace_ratio_branch                                                 */
/* DESCRICAO: Rastreia uma rama de dispersao beta/k0 ao longo da variacao de parametro do problema. */
/* Esta rotina implementa continuidade modal numerica entre pontos de         */
/* amostragem consecutivos em b_r/lambda_0.                                   */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; br_over_lambda: const std::vector<double> &; b:     */
/* double; seed_ratio: double; backend: ElementAssemblyBackend.               */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<SelectedRatioPoint> trace_ratio_branch(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    const std::vector<double> &br_over_lambda,
    double b,
    double seed_ratio,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr,
    const std::filesystem::path *linop_dir = nullptr,
    const std::string &prefix_base = "",
    const helmvec3_output::CaseDirs *out_dirs = nullptr)
{
    std::vector<SelectedRatioPoint> out;
    out.reserve(br_over_lambda.size());

    double prev = seed_ratio;
    for (double s : br_over_lambda)
    {
        const double k0 = 2.0 * M_PI * s / b;
        const std::string point_prefix =
            prefix_base.empty() ? std::string() : (prefix_base + "_br" + label_number(s));
        auto solve = solve_beta_point(
            mesh,
            eps,
            mu,
            k0,
            backend,
            perf,
            linop_dir,
            prefix_base.empty() ? nullptr : &point_prefix);
        if (solve.candidates.empty())
        {
            out.push_back(SelectedRatioPoint{});
            continue;
        }

        int best = 0;
        double best_err = std::abs(solve.candidates[0].beta_ratio - prev);
        for (int i = 1; i < (int)solve.candidates.size(); ++i)
        {
            const double e = std::abs(solve.candidates[i].beta_ratio - prev);
            if (e < best_err)
            {
                best_err = e;
                best = i;
            }
        }

        SelectedRatioPoint selected;
        selected.beta_ratio = solve.candidates[best].beta_ratio;
        selected.selected_candidate_rank = solve.candidates[best].candidate_rank;
        selected.selected_eig_index = solve.candidates[best].eig_index;
        selected.ez_ratio = solve.candidates[best].ez_ratio;
        selected.match_status = "tracked_branch";
        if (out_dirs != nullptr && !point_prefix.empty())
        {
            selected.spatial = helmvec3_field::export_mode(
                *out_dirs,
                mesh,
                solve.sys,
                solve.eig.VRcol,
                selected.selected_eig_index,
                point_prefix);
        }
        out.push_back(selected);
        prev = selected.beta_ratio;
    }
    return out;
}

struct Fig12CliConfig
{
    helmvec2::CoupledCliOptions cli;
    int nx = 10;
    int ny = 5;
};

struct Fig13CliConfig
{
    helmvec2::CoupledCliOptions cli;
    double d13_over_a = 0.20;
    int nx = 10;
    int ny = 5;
};

void print_output_dirs(const helmvec3_output::CaseDirs &out_dirs)
{
    std::cout << "[output] root_dir=\"" << out_dirs.root.string() << "\"\n";
    std::cout << "[output] csv_dir=\"" << out_dirs.csv.string() << "\"\n";
    std::cout << "[output] vtk_dir=\"" << out_dirs.vtk.string() << "\"\n";
    std::cout << "[output] img_dir=\"" << out_dirs.img.string() << "\"\n";
    std::cout << "[output] linop_dir=\"" << out_dirs.linop.string() << "\"\n";
}

Fig12CliConfig parse_fig12_cli(int argc, char **argv)
{
    Fig12CliConfig cfg;
    cfg.cli = helmvec2::parse_coupled_cli_options(argc, argv);

    if (!cfg.cli.positionals.empty())
        cfg.nx = std::atoi(cfg.cli.positionals[0].c_str());
    if (cfg.cli.positionals.size() >= 2)
        cfg.ny = std::atoi(cfg.cli.positionals[1].c_str());
    if (cfg.cli.positionals.size() >= 3)
    {
        const bool legacy_debug = (std::atoi(cfg.cli.positionals[2].c_str()) != 0);
        if (legacy_debug)
        {
            cfg.cli.debug_local_blocks = true;
            cfg.cli.debug_candidates = true;
        }
    }

    return cfg;
}

Fig13CliConfig parse_fig13_cli(int argc, char **argv)
{
    Fig13CliConfig cfg;
    cfg.cli = helmvec2::parse_coupled_cli_options(argc, argv);

    if (!cfg.cli.positionals.empty())
        cfg.d13_over_a = std::atof(cfg.cli.positionals[0].c_str());
    if (cfg.cli.positionals.size() >= 2)
        cfg.nx = std::atoi(cfg.cli.positionals[1].c_str());
    if (cfg.cli.positionals.size() >= 3)
        cfg.ny = std::atoi(cfg.cli.positionals[2].c_str());
    if (cfg.cli.positionals.size() >= 4)
    {
        const bool legacy_debug = (std::atoi(cfg.cli.positionals[3].c_str()) != 0);
        if (legacy_debug)
        {
            cfg.cli.debug_local_blocks = true;
            cfg.cli.debug_candidates = true;
        }
    }

    return cfg;
}

std::vector<Table10Block> make_table10_blocks()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {
        {0.0, {nan, 0.03, 0.52, 0.70, 0.79, 0.83, 0.88}, {nan, 0.04, 0.56, 0.71, 0.78, 0.83, 0.87}},
        {0.167, {nan, 0.21, 0.60, 0.72, 0.82, 0.88, 0.91}, {nan, 0.18, 0.59, 0.74, 0.81, 0.87, 0.90}},
        {0.286, {nan, 0.51, 0.78, 0.90, 0.99, 1.03, 1.10}, {nan, 0.44, 0.74, 0.88, 1.05, 1.03, 1.09}},
        {0.375, {nan, 0.68, 0.91, 1.05, 1.13, 1.20, 1.25}, {nan, 0.66, 0.90, 1.03, 1.11, 1.18, 1.23}},
        {0.5, {0.40, 0.90, 1.10, 1.20, 1.25, 1.30, 1.35}, {0.42, 0.89, 1.09, 1.19, 1.24, 1.31, 1.35}},
        {0.6, {0.70, 1.02, 1.18, 1.23, 1.31, 1.38, 1.41}, {0.67, 1.03, 1.19, 1.27, 1.33, 1.37, 1.40}},
        {0.8, {0.90, 1.18, 1.29, 1.38, 1.41, 1.43, 1.44}, {0.91, 1.18, 1.30, 1.37, 1.42, 1.44, 1.47}},
    };
}

} // namespace

int run_helmvec3_fig12_rect(int argc, char **argv)
{
    timing::Breakdown perf;
    timing::Stopwatch total_watch;
    const double a = 1.0;
    const double b = 0.45;
    const double d12 = 0.5 * b;
    const double eps_fill = 2.45;

    Fig12CliConfig cfg;
    try
    {
        cfg = parse_fig12_cli(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./helmvec3_fig12_rect [nx [ny [debug]]]"
                  << " [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    const auto out_dirs = helmvec3_output::ensure_case_dirs("fig12_rect");
    execution_log::ExecutionLogScope log_scope((out_dirs.root / "run.log").string());
    if (!log_scope.active())
    {
        std::cerr << "Aviso: nao foi possivel abrir run.log em "
                  << (out_dirs.root / "run.log")
                  << " (" << log_scope.error_message() << ")\n";
    }

    Mesh2D mesh = make_rect_mesh(a, b, cfg.nx, cfg.ny);
    std::vector<double> mu(mesh.tris.size(), 1.0);
    print_output_dirs(out_dirs);

    const std::vector<double> br_over_lambda_9 = {0.2, 0.3, 0.4, 0.5, 0.6};
    const std::vector<double> ref_ana_9 = {0.48, 1.00, 1.18, 1.26, 1.30};
    const std::vector<double> ref_hvec3_9 = {0.47, 1.01, 1.17, 1.28, 1.35};

    auto eps12 = helmvec23::eps_step_y(mesh, d12, eps_fill, 1.0);
    if (cfg.cli.debug_local_blocks)
    {
        const double k0_sample = 2.0 * M_PI * br_over_lambda_9.front() / b;
        print_first_triangle_closed_form_debug(mesh, k0_sample, eps12, mu);
    }

    auto ratio9 = match_ratio_to_reference(
        mesh,
        eps12,
        mu,
        br_over_lambda_9,
        b,
        ref_ana_9,
        cfg.cli.debug_candidates,
        cfg.cli.backend,
        &perf,
        &out_dirs.linop,
        "helmvec3_fig12_rect_figure12",
        &out_dirs);

    std::vector<helmvec3_output::Table9CsvRecord> table9_records;
    table9_records.reserve(br_over_lambda_9.size());

    std::cout << "[2.2.4] beta from given k0 (Figure 12)\n";
    std::cout << "a=" << a << " b=" << b << " d=" << d12 << " eps_fill=" << eps_fill
              << " nx=" << cfg.nx << " ny=" << cfg.ny << " tris=" << mesh.tris.size()
              << " backend=" << element_assembly_backend_name(cfg.cli.backend) << "\n";
    std::cout << "br/lambda0  beta/k0(FEM)  Analytic(ref)  HELMVEC3(ref)\n";
    for (int i = 0; i < (int)br_over_lambda_9.size(); ++i)
    {
        std::cout << br_over_lambda_9[i]
                  << "  " << ratio9[i].beta_ratio
                  << "  " << ref_ana_9[i]
                  << "  " << ref_hvec3_9[i] << "\n";
        const double err_ana =
            std::isfinite(ratio9[i].beta_ratio) ? error_metrics::absolute_relative_error_percent(ref_ana_9[i], ratio9[i].beta_ratio)
                                                : std::numeric_limits<double>::quiet_NaN();
        const double err_hvec3 =
            std::isfinite(ratio9[i].beta_ratio) ? error_metrics::absolute_relative_error_percent(ref_hvec3_9[i], ratio9[i].beta_ratio)
                                                : std::numeric_limits<double>::quiet_NaN();
        table9_records.push_back(
            {
                br_over_lambda_9[i],
                ratio9[i].beta_ratio,
                ref_ana_9[i],
                ref_hvec3_9[i],
                ratio9[i].selected_candidate_rank,
                ratio9[i].selected_eig_index,
                ratio9[i].ez_ratio,
                err_ana,
                err_hvec3,
                ratio9[i].match_status,
                ratio9[i].spatial.field_status,
                ratio9[i].spatial.et_fields_csv_file,
                ratio9[i].spatial.ez_fields_csv_file,
                ratio9[i].spatial.et_vtk_file,
                ratio9[i].spatial.ez_vtk_file,
            });
    }

    perf.total_ms = total_watch.elapsed_ms();
    perf.post_ms = std::max(0.0, perf.total_ms - perf.assembly_ms - perf.solve_ms);

    const auto table9_csv_path = out_dirs.csv / "helmvec3_fig12_rect_table9.csv";
    if (!helmvec3_output::write_table9_csv(table9_csv_path, table9_records))
        std::cerr << "Aviso: falha ao escrever " << table9_csv_path << "\n";

    run_timing_dispersion_csv::Record timing_record;
    timing_record.case_label = "helmvec3_fig12_rect";
    timing_record.geometry = "rect";
    timing_record.backend = element_assembly_backend_name(cfg.cli.backend);
    timing_record.debug_local_blocks = cfg.cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cfg.cli.debug_candidates ? 1 : 0;
    timing_record.a = a;
    timing_record.b = b;
    timing_record.d12 = d12;
    timing_record.d13_preview_over_a = 0.0;
    timing_record.eps_fill = eps_fill;
    timing_record.nx = cfg.nx;
    timing_record.ny = cfg.ny;
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.table9_sample_count = static_cast<int>(table9_records.size());
    timing_record.preview_sample_count = 0;
    timing_record.table10_block_count = 0;
    timing_record.table10_row_count = 0;
    timing_record.assembly_ms = perf.assembly_ms;
    timing_record.solve_ms = perf.solve_ms;
    timing_record.post_ms = perf.post_ms;
    timing_record.total_ms = perf.total_ms;
    const auto timing_csv_path = out_dirs.root / "run_timing.csv";
    if (!run_timing_dispersion_csv::write_csv(timing_csv_path.string(), timing_record))
        std::cerr << "Aviso: falha ao escrever " << timing_csv_path << "\n";

    timing::print_breakdown("helmvec3_fig12_rect", perf);
    return 0;
}

int run_helmvec3_fig13_rect(int argc, char **argv)
{
    timing::Breakdown perf;
    timing::Stopwatch total_watch;
    const double a = 1.0;
    const double b = 0.45;
    const double d12 = 0.5 * b;
    const double eps_fill = 2.45;

    Fig13CliConfig cfg;
    try
    {
        cfg = parse_fig13_cli(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./helmvec3_fig13_rect [d_over_a_preview [nx [ny [debug]]]]"
                  << " [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    const auto out_dirs = helmvec3_output::ensure_case_dirs("fig13_rect");
    execution_log::ExecutionLogScope log_scope((out_dirs.root / "run.log").string());
    if (!log_scope.active())
    {
        std::cerr << "Aviso: nao foi possivel abrir run.log em "
                  << (out_dirs.root / "run.log")
                  << " (" << log_scope.error_message() << ")\n";
    }

    Mesh2D mesh = make_rect_mesh(a, b, cfg.nx, cfg.ny);
    std::vector<double> mu(mesh.tris.size(), 1.0);
    print_output_dirs(out_dirs);

    const std::vector<double> br_over_lambda_9 = {0.2, 0.3, 0.4, 0.5, 0.6};
    const std::vector<double> br_over_lambda_10 = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const std::vector<Table10Block> table10 = make_table10_blocks();

    auto eps13_single = helmvec23::eps_step_x(mesh, cfg.d13_over_a * a, eps_fill, 1.0);
    if (cfg.cli.debug_local_blocks)
    {
        const double k0_sample = 2.0 * M_PI * br_over_lambda_9.front() / b;
        print_first_triangle_closed_form_debug(mesh, k0_sample, eps13_single, mu);
    }

    auto ratio13_preview = trace_ratio_branch(
        mesh,
        eps13_single,
        mu,
        br_over_lambda_9,
        b,
        0.5,
        cfg.cli.backend,
        &perf,
        &out_dirs.linop,
        std::string("helmvec3_fig13_rect_preview_da") + label_number(cfg.d13_over_a),
        &out_dirs);
    std::vector<helmvec3_output::PreviewCsvRecord> preview_records;
    preview_records.reserve(br_over_lambda_9.size());

    std::cout << "[2.2.4] beta from given k0 (Figure 13, single d/a preview)\n";
    std::cout << "a=" << a << " b=" << b << " eps_fill=" << eps_fill
              << " d/a=" << cfg.d13_over_a
              << " nx=" << cfg.nx << " ny=" << cfg.ny
              << " tris=" << mesh.tris.size()
              << " backend=" << element_assembly_backend_name(cfg.cli.backend) << "\n";
    std::cout << "br/lambda0  beta/k0(FEM branch)\n";
    for (int i = 0; i < (int)br_over_lambda_9.size(); ++i)
    {
        std::cout << br_over_lambda_9[i] << "  " << ratio13_preview[i].beta_ratio << "\n";
        preview_records.push_back(
            {
                cfg.d13_over_a,
                br_over_lambda_9[i],
                ratio13_preview[i].beta_ratio,
                ratio13_preview[i].selected_candidate_rank,
                ratio13_preview[i].selected_eig_index,
                ratio13_preview[i].ez_ratio,
                ratio13_preview[i].match_status,
                ratio13_preview[i].spatial.field_status,
                ratio13_preview[i].spatial.et_fields_csv_file,
                ratio13_preview[i].spatial.ez_fields_csv_file,
                ratio13_preview[i].spatial.et_vtk_file,
                ratio13_preview[i].spatial.ez_vtk_file,
            });
    }

    std::cout << "\n[2.2.4] Figure 13 / Table 10 validation\n";
    std::cout << "d/a  br/lambda0  beta/k0(FEM matched)  Analytical(ref)  HELMVEC3(ref)\n";
    std::vector<helmvec3_output::Table10CsvRecord> table10_records;
    for (const auto &blk : table10)
    {
        if (cfg.cli.debug_candidates)
            std::cout << "  [debug] Figure13 block d/a=" << blk.d_over_a << "\n";

        auto eps13 = helmvec23::eps_step_x(mesh, blk.d_over_a * a, eps_fill, 1.0);
        auto fem = match_ratio_to_reference(
            mesh,
            eps13,
            mu,
            br_over_lambda_10,
            b,
            blk.analytic_beta_over_k0,
            cfg.cli.debug_candidates,
            cfg.cli.backend,
            &perf,
            &out_dirs.linop,
            std::string("helmvec3_fig13_rect_table10_da") + label_number(blk.d_over_a),
            &out_dirs);

        for (int i = 0; i < (int)br_over_lambda_10.size(); ++i)
        {
            if (!std::isfinite(blk.analytic_beta_over_k0[i]) || !std::isfinite(blk.helmvec3_beta_over_k0[i]))
                continue;

            std::cout << blk.d_over_a
                      << "  " << br_over_lambda_10[i]
                      << "  " << fem[i].beta_ratio
                      << "  " << blk.analytic_beta_over_k0[i]
                      << "  " << blk.helmvec3_beta_over_k0[i] << "\n";

            const double err_ana =
                std::isfinite(fem[i].beta_ratio) ? error_metrics::absolute_relative_error_percent(blk.analytic_beta_over_k0[i], fem[i].beta_ratio)
                                                 : std::numeric_limits<double>::quiet_NaN();
            const double err_hvec3 =
                std::isfinite(fem[i].beta_ratio) ? error_metrics::absolute_relative_error_percent(blk.helmvec3_beta_over_k0[i], fem[i].beta_ratio)
                                                 : std::numeric_limits<double>::quiet_NaN();
            table10_records.push_back(
                {
                    blk.d_over_a,
                    br_over_lambda_10[i],
                    fem[i].beta_ratio,
                    blk.analytic_beta_over_k0[i],
                    blk.helmvec3_beta_over_k0[i],
                    fem[i].selected_candidate_rank,
                    fem[i].selected_eig_index,
                    fem[i].ez_ratio,
                    err_ana,
                    err_hvec3,
                    fem[i].match_status,
                    fem[i].spatial.field_status,
                    fem[i].spatial.et_fields_csv_file,
                    fem[i].spatial.ez_fields_csv_file,
                    fem[i].spatial.et_vtk_file,
                    fem[i].spatial.ez_vtk_file,
                });
        }
    }

    perf.total_ms = total_watch.elapsed_ms();
    perf.post_ms = std::max(0.0, perf.total_ms - perf.assembly_ms - perf.solve_ms);

    const auto preview_csv_path = out_dirs.csv / "helmvec3_fig13_rect_preview.csv";
    const auto table10_csv_path = out_dirs.csv / "helmvec3_fig13_rect_table10.csv";
    if (!helmvec3_output::write_preview_csv(preview_csv_path, preview_records))
        std::cerr << "Aviso: falha ao escrever " << preview_csv_path << "\n";
    if (!helmvec3_output::write_table10_csv(table10_csv_path, table10_records))
        std::cerr << "Aviso: falha ao escrever " << table10_csv_path << "\n";

    run_timing_dispersion_csv::Record timing_record;
    timing_record.case_label = "helmvec3_fig13_rect";
    timing_record.geometry = "rect";
    timing_record.backend = element_assembly_backend_name(cfg.cli.backend);
    timing_record.debug_local_blocks = cfg.cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cfg.cli.debug_candidates ? 1 : 0;
    timing_record.a = a;
    timing_record.b = b;
    timing_record.d12 = d12;
    timing_record.d13_preview_over_a = cfg.d13_over_a;
    timing_record.eps_fill = eps_fill;
    timing_record.nx = cfg.nx;
    timing_record.ny = cfg.ny;
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.table9_sample_count = 0;
    timing_record.preview_sample_count = static_cast<int>(preview_records.size());
    timing_record.table10_block_count = static_cast<int>(table10.size());
    timing_record.table10_row_count = static_cast<int>(table10_records.size());
    timing_record.assembly_ms = perf.assembly_ms;
    timing_record.solve_ms = perf.solve_ms;
    timing_record.post_ms = perf.post_ms;
    timing_record.total_ms = perf.total_ms;
    const auto timing_csv_path = out_dirs.root / "run_timing.csv";
    if (!run_timing_dispersion_csv::write_csv(timing_csv_path.string(), timing_record))
        std::cerr << "Aviso: falha ao escrever " << timing_csv_path << "\n";

    timing::print_breakdown("helmvec3_fig13_rect", perf);
    return 0;
}
