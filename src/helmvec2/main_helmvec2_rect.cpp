/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/main_helmvec2_rect.cpp                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter k0 dado beta.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.3, Eq.    */
/* (108)-(109), Fig. 11, Tabela 8.                                            */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC2 do         */
/* apendice em FORTRAN. A montagem global correspondente a Eq. (119) ocorre em*/
/* src/helmvec2/helmvec2_coupled_system.cpp. Comentarios priorizam didatica,  */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/error_metrics.hpp"
#include "core/execution_log.hpp"
#include "core/mesh2d.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/lapack_eig.hpp"
#include "core/run_timing_wavenumber_csv.hpp"
#include "core/spectral_csv.hpp"
#include "core/timing_utils.hpp"
#include "explicit/tri2d_coupled_explicit.hpp"
#include "helmvec23_shared.hpp"
#include "coupled_cli_options.hpp"
#include "helmvec2_case_output.hpp"
#include "helmvec2_coupled_system.hpp"
#include "helmvec2_field_output.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

struct ModeCand
{
    int eig_index = -1;
    double k0 = 0.0;
    double ez_ratio = 0.0; // ||Ez||^2 / (||Et||^2 + ||Ez||^2)
};

struct PhysicalRepresentative
{
    int candidate_rank = 0;
    int eig_index = -1;
    double k0 = 0.0;
    double ez_ratio = 0.0;
};

/******************************************************************************/
/* FUNCAO: print_block_3x3                                                    */
/* DESCRICAO: Imprime uma matriz 3x3 em formato compacto para depuracao       */
/* didatica dos blocos locais da formulacao acoplada.                         */
/* ENTRADA: name: const char *; M: const double[3][3].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static void print_block_3x3(const char *name, const double M[3][3])
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

/******************************************************************************/
/* FUNCAO: print_first_triangle_closed_form_debug                             */
/* DESCRICAO: Imprime os blocos locais closed-form do primeiro triangulo para */
/* a formulacao de 2.2.3, mostrando tanto as Eq. (120)-(125) quanto a forma  */
/* rearranjada usada no sistema global Eq. (119).                             */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps: const std::vector<double>*/
/* &; mu: const std::vector<double> &.                                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static void print_first_triangle_closed_form_debug(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps,
    const std::vector<double> &mu)
{
    if (mesh.tris.empty())
        return;

    const Tri &t = mesh.tris.front();
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    const auto blk125 = explicit_tri2d::tri2d_wavenumber_closed_form_eq_120_125(
        tg,
        beta,
        eps.front(),
        mu.front());
    const auto blk119 = explicit_tri2d::tri2d_wavenumber_rearranged_closed_form_eq_119(
        tg,
        beta,
        eps.front(),
        mu.front());

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form de 2.2.3\n";
    std::cout << "  beta_amostra=" << beta << "\n";
    std::cout << "  artigo Eq. (120)-(125):\n";
    print_block_3x3("Sel_tt", blk125.Sel_tt);
    print_block_3x3("Sel_tz", blk125.Sel_tz);
    print_block_3x3("Sel_zt", blk125.Sel_zt);
    print_block_3x3("Sel_zz", blk125.Sel_zz);
    print_block_3x3("Tel_tt", blk125.Tel_tt);
    print_block_3x3("Tel_zz", blk125.Tel_zz);

    std::cout << "  rearranjo usado no codigo para Eq. (119):\n";
    print_block_3x3("A_tt", blk119.A_tt);
    print_block_3x3("A_tz", blk119.A_tz);
    print_block_3x3("A_zt", blk119.A_zt);
    print_block_3x3("A_zz", blk119.A_zz);
    print_block_3x3("B_tt", blk119.B_tt);
    print_block_3x3("B_zz", blk119.B_zz);
}

/******************************************************************************/
/* FUNCAO: collect_mode_candidates                                            */
/* DESCRICAO: Extrai candidatos fisicos do espectro de A x = lambda B x,      */
/* usando lambda = k0^2 (Secao 2.2.3), e estima energia relativa no bloco Ez. */
/* ENTRADA: res: const GenEigGeneralResult &; nt: int; nz: int; imag_tol:     */
/* double.                                                                    */
/* SAIDA: std::vector<ModeCand>.                                              */
/******************************************************************************/
static std::vector<ModeCand> collect_mode_candidates(
    const GenEigGeneralResult &res,
    int nt,
    int nz,
    double imag_tol = 1e-7)
{
    std::vector<ModeCand> out;
    const int n = res.n;
    for (int i = 0; i < res.n; ++i)
    {
        if (!std::isfinite(res.lambda_re[i]))
            continue;
        if (std::abs(res.lambda_im[i]) > imag_tol)
            continue;
        if (res.lambda_re[i] <= 1e-10)
            continue;

        const double k0 = std::sqrt(res.lambda_re[i]);
        double et = 0.0, ez = 0.0;
        for (int r = 0; r < nt; ++r)
        {
            const double v = res.VRcol[(size_t)i * n + r];
            et += v * v;
        }
        for (int r = 0; r < nz; ++r)
        {
            const double v = res.VRcol[(size_t)i * n + (nt + r)];
            ez += v * v;
        }
        const double den = et + ez;
        const double ratio = (den > 0.0) ? (ez / den) : 0.0;
        out.push_back({i, k0, ratio});
    }
    std::sort(out.begin(), out.end(), [](const ModeCand &a, const ModeCand &b)
              { return a.k0 < b.k0; });
    return out;
}

/******************************************************************************/
/* FUNCAO: unique_physical_representatives                                    */
/* DESCRICAO: Condensa candidatos fisicos quase degenerados em uma lista      */
/* unica por raiz de k0, preservando um representante real para exportacao    */
/* espacial. Entre candidatos com mesmo k0, escolhe o de maior Ez-ratio.      */
/* ENTRADA: all_modes: const std::vector<ModeCand> &; k0_min_phys: double;    */
/* uniq_tol: double.                                                          */
/* SAIDA: std::vector<PhysicalRepresentative>.                                */
/******************************************************************************/
static std::vector<PhysicalRepresentative> unique_physical_representatives(
    const std::vector<ModeCand> &all_modes,
    double k0_min_phys,
    double uniq_tol = 1e-8)
{
    std::vector<PhysicalRepresentative> reps;
    for (int candidate_rank = 0; candidate_rank < (int)all_modes.size(); ++candidate_rank)
    {
        const ModeCand &mode = all_modes[(size_t)candidate_rank];
        if (mode.k0 <= k0_min_phys)
            continue;

        if (!reps.empty() && std::abs(mode.k0 - reps.back().k0) < uniq_tol)
        {
            if (mode.ez_ratio > reps.back().ez_ratio)
            {
                reps.back().candidate_rank = candidate_rank + 1;
                reps.back().eig_index = mode.eig_index;
                reps.back().k0 = mode.k0;
                reps.back().ez_ratio = mode.ez_ratio;
            }
            continue;
        }

        reps.push_back(
            {
                candidate_rank + 1,
                mode.eig_index,
                mode.k0,
                mode.ez_ratio,
            });
    }
    return reps;
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
    // Secao 2.2.3.5 / Figura 11: quadrado, parte superior eps=1.0 e inferior eps=1.5.
    // Objetivo: reproduzir Tabela 8 com beta fixo e extrair k0L dos modos LSM.
    const double L = 1.0;
    int nx = 6, ny = 6; // 72 triangles
    double beta = 10.0; // beta*L = 10
    helmvec2::CoupledCliOptions cli;
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso: ./helmvec2_rect [beta [nx ny [debug]]]"
           << " [--backend closed-form|gauss]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--beta BETA] [--nx NX] [--ny NY]"
           << " (nao misture com os posicionais principais)\n";
        os << "Compatibilidade: os posicionais principais continuam aceitos, mas estao deprecated;"
           << " prefira os aliases nomeados acima.\n";
        os << "Compatibilidade: o debug posicional legado [debug] continua aceito, mas esta deprecated;"
           << " prefira --debug-local-blocks e/ou --debug-candidates.\n";
    };
    if (helmvec2::coupled_cli_requests_help(argc, argv))
    {
        print_usage(std::cout);
        return 0;
    }
    try
    {
        cli = helmvec2::parse_coupled_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const bool has_named_primary_args =
        cli.beta_was_provided ||
        cli.nx_was_provided ||
        cli.ny_was_provided;
    const bool used_legacy_positionals =
        !has_named_primary_args && !cli.positionals.empty();
    const bool used_legacy_debug =
        !has_named_primary_args && cli.positionals.size() >= 4;

    if (cli.d_over_a_preview_was_provided)
    {
        std::cerr << "Erro: helmvec2_rect nao aceita --d-over-a-preview; use --beta.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_primary_args && !cli.positionals.empty())
    {
        std::cerr << "Erro: nao misture aliases nomeados principais com argumentos posicionais principais em helmvec2_rect; escolha apenas um estilo de chamada.\n";
        print_usage(std::cerr);
        return 2;
    }

    if (has_named_primary_args)
    {
        if (cli.beta_was_provided)
            beta = cli.beta;
        if (cli.nx_was_provided)
            nx = cli.nx;
        if (cli.ny_was_provided)
            ny = cli.ny;
    }
    else
    {
        if (!cli.positionals.empty())
            beta = std::atof(cli.positionals[0].c_str());
        if (cli.positionals.size() >= 3)
        {
            nx = std::atoi(cli.positionals[1].c_str());
            ny = std::atoi(cli.positionals[2].c_str());
        }
        if (cli.positionals.size() >= 4)
        {
            const bool legacy_debug = (std::atoi(cli.positionals[3].c_str()) != 0);
            if (legacy_debug)
            {
                cli.debug_local_blocks = true;
                cli.debug_candidates = true;
            }
        }
    }

    const auto out_dirs = helmvec2_output::ensure_case_dirs("rect");
    execution_log::ExecutionLogScope log_scope((out_dirs.root / "run.log").string());
    if (!log_scope.active())
    {
        std::cerr << "Aviso: nao foi possivel abrir run.log em "
                  << (out_dirs.root / "run.log")
                  << " (" << log_scope.error_message() << ")\n";
    }
    if (used_legacy_positionals)
    {
        std::cerr << "Aviso: os argumentos posicionais principais de helmvec2_rect continuam aceitos por compatibilidade, mas estao deprecated; prefira --beta/--nx/--ny.\n";
    }
    if (used_legacy_debug)
    {
        std::cerr << "Aviso: o debug posicional legado de helmvec2_rect continua aceito por compatibilidade, mas esta deprecated; prefira --debug-local-blocks e/ou --debug-candidates.\n";
    }

    auto mesh = make_rect_mesh(L, L, nx, ny);
    auto eps = helmvec23::eps_step_y(mesh, 0.5 * L, 1.5, 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);

    std::cout << "[output] root_dir=" << out_dirs.root << "\n";
    std::cout << "[output] csv_dir=" << out_dirs.csv << "\n";
    std::cout << "[output] vtk_dir=" << out_dirs.vtk << "\n";
    std::cout << "[output] img_dir=" << out_dirs.img << "\n";
    std::cout << "[output] linop_dir=" << out_dirs.linop << "\n";

    if (cli.debug_local_blocks)
        print_first_triangle_closed_form_debug(mesh, beta, eps, mu);

    timing::Stopwatch stage;
    auto sys = tp3485::build_eq119_helmvec2_system_E(mesh, beta, eps, mu, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    auto eig = generalized_eigs_real_vec(sys.A, sys.B);
    perf.solve_ms += stage.elapsed_ms();
    if (!spectral_csv::write_general_problem_exports(
            out_dirs.linop,
            "helmvec2_rect",
            sys.A,
            sys.B,
            eig) ||
        !spectral_csv::write_dense_crs_csv(out_dirs.linop / "helmvec2_rect_A_tt_crs.csv", sys.A_tt) ||
        !spectral_csv::write_rect_like_crs_csv(
            out_dirs.linop / "helmvec2_rect_A_tz_crs.csv",
            sys.A_tz.nr,
            sys.A_tz.nc,
            [&](int i, int j)
            { return sys.A_tz(i, j); }) ||
        !spectral_csv::write_rect_like_crs_csv(
            out_dirs.linop / "helmvec2_rect_A_zt_crs.csv",
            sys.A_zt.nr,
            sys.A_zt.nc,
            [&](int i, int j)
            { return sys.A_zt(i, j); }) ||
        !spectral_csv::write_dense_crs_csv(out_dirs.linop / "helmvec2_rect_A_zz_crs.csv", sys.A_zz) ||
        !spectral_csv::write_dense_crs_csv(out_dirs.linop / "helmvec2_rect_B_tt_crs.csv", sys.B_tt) ||
        !spectral_csv::write_dense_crs_csv(out_dirs.linop / "helmvec2_rect_B_zz_crs.csv", sys.B_zz))
    {
        std::cerr << "Aviso: falha ao escrever artefatos espectrais em " << out_dirs.linop << "\n";
        return 4;
    }
    auto all_modes = collect_mode_candidates(eig, sys.nt, sys.nz);

    double eps_max = 1.0;
    for (double e : eps)
        eps_max = std::max(eps_max, e);
    const double k0_min_phys = beta / std::sqrt(eps_max);
    // Regiao propagante: beta < k0*sqrt(eps_r_max)  => k0 > beta/sqrt(eps_r_max).
    // Esse filtro remove raizes nao propagantes para o meio estratificado da Fig. 11.

    // Mantem o espectro bruto para auditoria e, em seguida, cria uma lista
    // unica de representantes fisicos para o matching da Tabela 8.
    std::vector<helmvec2_output::CandidateCsvRecord> candidate_records;
    candidate_records.reserve(all_modes.size());
    for (int candidate_rank = 0; candidate_rank < (int)all_modes.size(); ++candidate_rank)
    {
        const auto &m = all_modes[(size_t)candidate_rank];
        const bool passes = (m.k0 > k0_min_phys);
        candidate_records.push_back(
            {
                candidate_rank + 1,
                m.eig_index,
                m.k0,
                m.k0 * L,
                m.ez_ratio,
                passes ? 1 : 0,
            });
    }
    const std::vector<PhysicalRepresentative> physical_modes =
        unique_physical_representatives(all_modes, k0_min_phys);
    std::vector<double> k0_phys;
    k0_phys.reserve(physical_modes.size());
    for (const auto &mode : physical_modes)
        k0_phys.push_back(mode.k0);

    const std::vector<double> ref_helmvec2 = {8.8150, 9.4430, 10.3500, 11.1410, 11.2890, 11.4246, 12.1460, 12.5894, 12.8237, 12.9987};
    const std::vector<double> ref_hayata = {8.8093, 9.3896, 10.2752, 11.1030, 11.2677, 11.4501, 11.9882, 12.6686, 12.8092, 12.9575};

    std::cout << "[2.2.3] wave-number from given beta (Figure 11)\n";
    std::cout << "L=" << L << " beta=" << beta << " beta*L=" << beta * L
              << " nx=" << nx << " ny=" << ny
              << " tris=" << mesh.tris.size()
              << " k0_min_phys~" << k0_min_phys
              << " backend=" << element_assembly_backend_name(cli.backend) << "\n\n";

    std::cout << "first raw roots (k0L), after physical filter:\n";
    for (int i = 0; i < (int)k0_phys.size() && i < 14; ++i)
        std::cout << " " << (i + 1) << "  " << k0_phys[i] << "\n";
    std::cout << "\n";

    std::cout << "mode  k0L(FEM matched)  HELMVEC2(ref)  Hayata(ref)\n";
    std::vector<char> used(physical_modes.size(), 0);
    std::vector<helmvec2_output::ModeCsvRecord> mode_records;
    const int N = std::min<int>(10, ref_helmvec2.size());
    for (int i = 0; i < N; ++i)
    {
        // Modo de validacao: escolhe a raiz calculada nao usada mais proxima da
        // familia modal publicada do HELMVEC2 (Tabela 8), preservando ordenacao.
        const int matched_index =
            helmvec23::pick_closest_unused_index(ref_helmvec2[i], k0_phys, used);
        const double k0L =
            (matched_index >= 0) ? k0_phys[(size_t)matched_index] : std::numeric_limits<double>::quiet_NaN();
        std::cout << (i + 1) << "  " << k0L
                  << "  " << ref_helmvec2[i]
                  << "  " << ref_hayata[i] << "\n";
        const double err_helmvec2 =
            std::isfinite(k0L) ? error_metrics::absolute_relative_error_percent(ref_helmvec2[i], k0L)
                               : std::numeric_limits<double>::quiet_NaN();
        const double err_hayata =
            std::isfinite(k0L) ? error_metrics::absolute_relative_error_percent(ref_hayata[i], k0L)
                               : std::numeric_limits<double>::quiet_NaN();
        helmvec2_output::ModeCsvRecord rec;
        rec.mode = i + 1;
        rec.matched_candidate_rank =
            (matched_index >= 0) ? physical_modes[(size_t)matched_index].candidate_rank : 0;
        rec.matched_eig_index =
            (matched_index >= 0) ? physical_modes[(size_t)matched_index].eig_index : 0;
        rec.k0_fem_matched =
            std::isfinite(k0L) ? k0L / L : std::numeric_limits<double>::quiet_NaN();
        rec.k0l_fem_matched = k0L;
        rec.ez_ratio =
            (matched_index >= 0) ? physical_modes[(size_t)matched_index].ez_ratio
                                 : std::numeric_limits<double>::quiet_NaN();
        rec.ref_helmvec2 = ref_helmvec2[i];
        rec.ref_hayata = ref_hayata[i];
        rec.error_percent_helmvec2 = err_helmvec2;
        rec.error_percent_hayata = err_hayata;
        rec.match_status = (matched_index >= 0) ? "matched" : "missing_candidate";

        if (matched_index >= 0)
        {
            const auto artifacts = helmvec2_field::export_mode(
                out_dirs,
                mesh,
                sys,
                eig.VRcol,
                physical_modes[(size_t)matched_index].eig_index,
                i + 1,
                physical_modes[(size_t)matched_index].candidate_rank);
            rec.field_status = artifacts.field_status;
            rec.et_fields_csv_file = artifacts.et_fields_csv_file;
            rec.ez_fields_csv_file = artifacts.ez_fields_csv_file;
            rec.et_vtk_file = artifacts.et_vtk_file;
            rec.ez_vtk_file = artifacts.ez_vtk_file;
        }

        mode_records.push_back(
            rec);
    }

    if (cli.debug_candidates)
    {
        std::cout << "\n[debug] candidate modes after k0_min filter:\n";
        std::cout << "k0L  ez_ratio\n";
        for (const auto &m : all_modes)
        {
            if (m.k0 <= k0_min_phys)
                continue;
            std::cout << m.k0 * L << "  " << m.ez_ratio << "\n";
        }
    }

    stage.reset();
    const auto modes_csv_path = out_dirs.csv / "helmvec2_rect_modes.csv";
    const auto candidates_csv_path = out_dirs.csv / "helmvec2_rect_candidates.csv";
    if (!helmvec2_output::write_modes_csv(modes_csv_path, mode_records))
    {
        std::cerr << "Aviso: falha ao escrever " << modes_csv_path << "\n";
    }
    if (!helmvec2_output::write_candidates_csv(candidates_csv_path, candidate_records))
    {
        std::cerr << "Aviso: falha ao escrever " << candidates_csv_path << "\n";
    }
    perf.post_ms += stage.elapsed_ms();

    run_timing_wavenumber_csv::Record timing_record;
    timing_record.case_label = "helmvec2_rect";
    timing_record.geometry = "rect";
    timing_record.backend = element_assembly_backend_name(cli.backend);
    timing_record.debug_local_blocks = cli.debug_local_blocks ? 1 : 0;
    timing_record.debug_candidates = cli.debug_candidates ? 1 : 0;
    timing_record.beta = beta;
    timing_record.L = L;
    timing_record.betaL = beta * L;
    timing_record.nx = nx;
    timing_record.ny = ny;
    timing_record.mesh_nodes = static_cast<int>(mesh.nodes.size());
    timing_record.mesh_tris = static_cast<int>(mesh.tris.size());
    timing_record.eps_top = 1.0;
    timing_record.eps_bottom = 1.5;
    timing_record.mu_r = 1.0;
    timing_record.k0_min_phys = k0_min_phys;
    timing_record.candidate_count_raw = static_cast<int>(all_modes.size());
    timing_record.candidate_count_phys = static_cast<int>(k0_phys.size());
    timing_record.matched_mode_count = static_cast<int>(mode_records.size());
    timing_record.assembly_ms = perf.assembly_ms;
    timing_record.solve_ms = perf.solve_ms;
    timing_record.post_ms = perf.post_ms;
    timing_record.total_ms = total_watch.elapsed_ms();
    const auto timing_csv_path = out_dirs.root / "run_timing.csv";
    if (!run_timing_wavenumber_csv::write_csv(timing_csv_path.string(), timing_record))
    {
        std::cerr << "Aviso: falha ao escrever " << timing_csv_path << "\n";
    }

    perf.total_ms = timing_record.total_ms;
    timing::print_breakdown("helmvec2_rect", perf);
    return 0;
}
