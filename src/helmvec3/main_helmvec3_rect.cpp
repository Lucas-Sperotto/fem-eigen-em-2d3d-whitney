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
#include <cctype>
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
    int raw_row_index = -1;
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
    std::vector<helmvec3_output::RawSpectrumCsvRecord> raw_spectrum;
    std::vector<helmvec3_output::MatrixAuditCsvRecord> matrix_audit;
};

struct MatrixBlockAuditMetrics
{
    double P_fro = 0.0;
    double Q_fro = 0.0;
    double P_tt_fro = 0.0;
    double P_zz_fro = 0.0;
    double Q_tt_fro = 0.0;
    double Q_tz_fro = 0.0;
    double Q_zt_fro = 0.0;
    double Q_zz_fro = 0.0;
    double P_tt_asym_rel = 0.0;
    double P_zz_asym_rel = 0.0;
    double Q_tt_asym_rel = 0.0;
    double Q_zz_asym_rel = 0.0;
    double Q_tz_Q_zt_transpose_mismatch_rel = 0.0;
    double block_norm_max = 0.0;
    double P_tt_scale_rel = 0.0;
    double P_zz_scale_rel = 0.0;
    double Q_tt_scale_rel = 0.0;
    double Q_tz_scale_rel = 0.0;
    double Q_zt_scale_rel = 0.0;
    double Q_zz_scale_rel = 0.0;
};

bool env_flag_enabled(const char *name)
{
    const char *raw = std::getenv(name);
    if (raw == nullptr)
        return false;

    std::string value(raw);
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char ch)
        { return static_cast<char>(std::tolower(ch)); });
    return !(value.empty() || value == "0" || value == "false" || value == "off" || value == "no");
}

std::string normalize_token(std::string token)
{
    std::transform(
        token.begin(),
        token.end(),
        token.begin(),
        [](unsigned char ch)
        { return static_cast<char>(std::tolower(ch)); });
    return token;
}

CoupledBetaDiagVariant env_beta_diag_variant()
{
    const char *raw = std::getenv("TP3485_HELMVEC3_DIAG_BETA_VARIANT");
    if (raw == nullptr || std::string(raw).empty())
        return CoupledBetaDiagVariant::Baseline;

    const std::string token = normalize_token(std::string(raw));
    if (token == "baseline")
        return CoupledBetaDiagVariant::Baseline;
    if (token == "diag_eq141_blend_half_qtt" || token == "eq141_blend_half" || token == "eq141_half")
        return CoupledBetaDiagVariant::DiagEq141BlendHalfQtt;
    if (token == "diag_eq141_eps_mass_qtt" || token == "eq141" || token == "eq141_qtt")
        return CoupledBetaDiagVariant::DiagEq141EpsMassQtt;
    if (token == "diag_eq142_doc_qzz" || token == "eq142" || token == "eq142_qzz")
        return CoupledBetaDiagVariant::DiagEq142DocQzz;
    if (token == "diag_scale_pzz_double" || token == "pzz_double")
        return CoupledBetaDiagVariant::DiagScalePzzDouble;
    if (token == "diag_scale_qzz_half" || token == "qzz_half")
        return CoupledBetaDiagVariant::DiagScaleQzzHalf;
    if (token == "diag_scale_coupling_double" || token == "coupling_double")
        return CoupledBetaDiagVariant::DiagScaleCouplingDouble;
    if (token == "diag_eq141_eps_mass_qtt_plus_pzz_double" || token == "eq141_qtt_plus_pzz_double")
        return CoupledBetaDiagVariant::DiagEq141EpsMassQttPlusPzzDouble;
    if (token == "diag_eq141_eps_mass_qtt_plus_qzz_half" || token == "eq141_qtt_plus_qzz_half")
        return CoupledBetaDiagVariant::DiagEq141EpsMassQttPlusQzzHalf;

    throw std::runtime_error(
        "Valor invalido em TP3485_HELMVEC3_DIAG_BETA_VARIANT: " + std::string(raw));
}

bool nearly_equal(double lhs, double rhs, double tol = 1.0e-12)
{
    return std::abs(lhs - rhs) <= tol;
}

std::string figure13_priority_point_key(double d_over_a, double br_over_lambda0)
{
    if (nearly_equal(d_over_a, 0.167) && nearly_equal(br_over_lambda0, 0.5))
        return "P1";
    if (nearly_equal(d_over_a, 0.286) && nearly_equal(br_over_lambda0, 0.5))
        return "P2";
    if (nearly_equal(d_over_a, 0.5) && nearly_equal(br_over_lambda0, 0.4))
        return "P3";
    return "";
}

bool should_export_raw_spectrum_diag(double d_over_a, double br_over_lambda0)
{
    return !figure13_priority_point_key(d_over_a, br_over_lambda0).empty();
}

double real_positive_beta_ratio(double lambda_real, double k0)
{
    if (!std::isfinite(lambda_real) || lambda_real <= 0.0 || !std::isfinite(k0) || k0 == 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(lambda_real) / k0;
}

void mark_selected_raw_spectrum(
    std::vector<helmvec3_output::RawSpectrumCsvRecord> &rows,
    int selected_eig_index)
{
    for (auto &row : rows)
    {
        if (row.solver_index == selected_eig_index)
            row.selected_after_matching = 1;
    }
}

void mark_selected_matrix_audit(
    std::vector<helmvec3_output::MatrixAuditCsvRecord> &rows,
    int selected_eig_index)
{
    for (auto &row : rows)
    {
        if (row.solver_index == selected_eig_index)
            row.selected_after_matching = 1;
    }
}

double vector_l2_norm(const std::vector<double> &x)
{
    double acc = 0.0;
    for (double value : x)
        acc += value * value;
    return std::sqrt(acc);
}

double dense_frobenius_norm(const DenseMat &mat)
{
    double acc = 0.0;
    for (double value : mat.a)
        acc += value * value;
    return std::sqrt(acc);
}

double rect_frobenius_norm(const DenseRectMat &mat)
{
    double acc = 0.0;
    for (double value : mat.a)
        acc += value * value;
    return std::sqrt(acc);
}

double dense_relative_asymmetry(const DenseMat &mat)
{
    double diff_acc = 0.0;
    for (int i = 0; i < mat.n; ++i)
    {
        for (int j = 0; j < mat.n; ++j)
        {
            const double diff = mat(i, j) - mat(j, i);
            diff_acc += diff * diff;
        }
    }
    const double den = std::max(dense_frobenius_norm(mat), 1.0e-30);
    return std::sqrt(diff_acc) / den;
}

double rect_transpose_mismatch_rel(const DenseRectMat &lhs, const DenseRectMat &rhs)
{
    if (lhs.nr != rhs.nc || lhs.nc != rhs.nr)
        return std::numeric_limits<double>::infinity();

    double diff_acc = 0.0;
    for (int i = 0; i < lhs.nr; ++i)
    {
        for (int j = 0; j < lhs.nc; ++j)
        {
            const double diff = lhs(i, j) - rhs(j, i);
            diff_acc += diff * diff;
        }
    }
    const double den = std::max(std::max(rect_frobenius_norm(lhs), rect_frobenius_norm(rhs)), 1.0e-30);
    return std::sqrt(diff_acc) / den;
}

void dense_matvec(const DenseMat &mat, const std::vector<double> &x, std::vector<double> &y)
{
    y.assign((size_t)mat.n, 0.0);
    for (int i = 0; i < mat.n; ++i)
    {
        double acc = 0.0;
        for (int j = 0; j < mat.n; ++j)
            acc += mat(i, j) * x[(size_t)j];
        y[(size_t)i] = acc;
    }
}

void binding_vector_parts(
    const GenEigGeneralResult &eig,
    const spectral_csv::GeneralEigenColumnBinding &binding,
    std::vector<double> &real_part,
    std::vector<double> &imag_part)
{
    real_part.assign((size_t)eig.n, 0.0);
    imag_part.assign((size_t)eig.n, 0.0);
    for (int dof = 0; dof < eig.n; ++dof)
    {
        real_part[(size_t)dof] = eig.VRcol[(size_t)binding.real_column * (size_t)eig.n + (size_t)dof];
        if (binding.imag_column >= 0 && binding.imag_sign != 0)
            imag_part[(size_t)dof] =
                binding.imag_sign * eig.VRcol[(size_t)binding.imag_column * (size_t)eig.n + (size_t)dof];
    }
}

struct GeneralizedResidualAudit
{
    double abs_norm = 0.0;
    double rel_norm = 0.0;
};

GeneralizedResidualAudit compute_generalized_residual(
    const DenseMat &P,
    const DenseMat &Q,
    const GenEigGeneralResult &eig,
    const spectral_csv::GeneralEigenColumnBinding &binding,
    double P_fro,
    double Q_fro)
{
    std::vector<double> xr;
    std::vector<double> xi;
    binding_vector_parts(eig, binding, xr, xi);

    std::vector<double> Pxr;
    std::vector<double> Pxi;
    std::vector<double> Qxr;
    std::vector<double> Qxi;
    dense_matvec(P, xr, Pxr);
    dense_matvec(P, xi, Pxi);
    dense_matvec(Q, xr, Qxr);
    dense_matvec(Q, xi, Qxi);

    double abs_acc = 0.0;
    for (int i = 0; i < eig.n; ++i)
    {
        const double rr =
            Pxr[(size_t)i] - binding.lambda_real * Qxr[(size_t)i] + binding.lambda_imag * Qxi[(size_t)i];
        const double ri =
            Pxi[(size_t)i] - binding.lambda_real * Qxi[(size_t)i] - binding.lambda_imag * Qxr[(size_t)i];
        abs_acc += rr * rr + ri * ri;
    }

    const double abs_norm = std::sqrt(abs_acc);
    const double x_norm = std::hypot(vector_l2_norm(xr), vector_l2_norm(xi));
    const double den =
        (P_fro + std::hypot(binding.lambda_real, binding.lambda_imag) * Q_fro) * std::max(x_norm, 1.0e-30);

    GeneralizedResidualAudit out;
    out.abs_norm = abs_norm;
    out.rel_norm = (den > 0.0) ? (abs_norm / den) : abs_norm;
    return out;
}

MatrixBlockAuditMetrics compute_matrix_block_audit_metrics(const CoupledBetaSystem &sys)
{
    MatrixBlockAuditMetrics out;
    out.P_fro = dense_frobenius_norm(sys.P);
    out.Q_fro = dense_frobenius_norm(sys.Q);
    out.P_tt_fro = dense_frobenius_norm(sys.P_tt);
    out.P_zz_fro = dense_frobenius_norm(sys.P_zz);
    out.Q_tt_fro = dense_frobenius_norm(sys.Q_tt);
    out.Q_tz_fro = rect_frobenius_norm(sys.Q_tz);
    out.Q_zt_fro = rect_frobenius_norm(sys.Q_zt);
    out.Q_zz_fro = dense_frobenius_norm(sys.Q_zz);
    out.P_tt_asym_rel = dense_relative_asymmetry(sys.P_tt);
    out.P_zz_asym_rel = dense_relative_asymmetry(sys.P_zz);
    out.Q_tt_asym_rel = dense_relative_asymmetry(sys.Q_tt);
    out.Q_zz_asym_rel = dense_relative_asymmetry(sys.Q_zz);
    out.Q_tz_Q_zt_transpose_mismatch_rel = rect_transpose_mismatch_rel(sys.Q_tz, sys.Q_zt);
    out.block_norm_max = std::max(
        std::max(std::max(out.P_tt_fro, out.P_zz_fro), std::max(out.Q_tt_fro, out.Q_tz_fro)),
        std::max(std::max(out.Q_zt_fro, out.Q_zz_fro), 1.0e-30));
    out.P_tt_scale_rel = out.P_tt_fro / out.block_norm_max;
    out.P_zz_scale_rel = out.P_zz_fro / out.block_norm_max;
    out.Q_tt_scale_rel = out.Q_tt_fro / out.block_norm_max;
    out.Q_tz_scale_rel = out.Q_tz_fro / out.block_norm_max;
    out.Q_zt_scale_rel = out.Q_zt_fro / out.block_norm_max;
    out.Q_zz_scale_rel = out.Q_zz_fro / out.block_norm_max;
    return out;
}

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
    CoupledBetaDiagVariant diag_variant,
    timing::Breakdown *perf = nullptr,
    const std::filesystem::path *linop_dir = nullptr,
    const std::string *artifact_prefix = nullptr,
    bool collect_raw_spectrum_diag = false,
    bool collect_matrix_audit_diag = false)
{
    timing::Stopwatch stage;
    BetaPointSolve result;
    result.sys = tp3485::build_eq136_helmvec3_system_E(mesh, k0, eps, mu, backend, diag_variant);
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
        int raw_row_index = -1;
        double beta_ratio = 0.0;
        double ez_ratio = 0.0;
    };

    const int n = result.eig.n;
    std::vector<int> solver_to_raw_row_index((size_t)n, -1);
    std::vector<int> solver_to_matrix_row_index((size_t)n, -1);
    const bool collect_any_diag = collect_raw_spectrum_diag || collect_matrix_audit_diag;
    const MatrixBlockAuditMetrics block_metrics =
        collect_matrix_audit_diag ? compute_matrix_block_audit_metrics(result.sys) : MatrixBlockAuditMetrics{};
    if (collect_any_diag)
    {
        const auto bindings = spectral_csv::general_eigen_bindings(result.eig);
        if (collect_raw_spectrum_diag)
            result.raw_spectrum.reserve(bindings.size());
        if (collect_matrix_audit_diag)
            result.matrix_audit.reserve(bindings.size());
        for (size_t rank = 0; rank < bindings.size(); ++rank)
        {
            const auto &binding = bindings[rank];
            std::vector<double> xr;
            std::vector<double> xi;
            binding_vector_parts(result.eig, binding, xr, xi);

            double et = 0.0;
            double ez = 0.0;
            for (int dof = 0; dof < result.sys.nt; ++dof)
            {
                const double real_part = xr[(size_t)dof];
                const double imag_part = xi[(size_t)dof];
                et += real_part * real_part + imag_part * imag_part;
            }
            for (int dof = 0; dof < result.sys.nz; ++dof)
            {
                const int col_offset = result.sys.nt + dof;
                const double real_part = xr[(size_t)col_offset];
                const double imag_part = xi[(size_t)col_offset];
                ez += real_part * real_part + imag_part * imag_part;
            }
            const double den = et + ez;
            const double ratio_ez = (den > 0.0) ? (ez / den) : 0.0;
            if (collect_raw_spectrum_diag)
            {
                helmvec3_output::RawSpectrumCsvRecord row;
                row.section = "2.2.4";
                row.case_label = "Figure13_Table10";
                row.ordered_rank = static_cast<int>(rank) + 1;
                row.solver_index = binding.original_index;
                row.lambda_real = binding.lambda_real;
                row.lambda_imag = binding.lambda_imag;
                row.beta_ratio_if_real_positive = real_positive_beta_ratio(binding.lambda_real, k0);
                row.filter_reason = "pending_filter";
                row.et_energy = et;
                row.ez_energy = ez;
                row.ez_ratio = ratio_ez;
                result.raw_spectrum.push_back(row);
                if (binding.original_index >= 0 && binding.original_index < n)
                    solver_to_raw_row_index[(size_t)binding.original_index] =
                        static_cast<int>(result.raw_spectrum.size()) - 1;
            }
            if (collect_matrix_audit_diag)
            {
                const GeneralizedResidualAudit residual =
                    compute_generalized_residual(result.sys.P, result.sys.Q, result.eig, binding, block_metrics.P_fro, block_metrics.Q_fro);
                helmvec3_output::MatrixAuditCsvRecord row;
                row.section = "2.2.4";
                row.case_label = "Figure13_Table10";
                row.backend = element_assembly_backend_name(backend);
                row.ordered_rank = static_cast<int>(rank) + 1;
                row.solver_index = binding.original_index;
                row.lambda_real = binding.lambda_real;
                row.lambda_imag = binding.lambda_imag;
                row.beta_ratio_if_real_positive = real_positive_beta_ratio(binding.lambda_real, k0);
                row.filter_reason = "pending_filter";
                row.et_energy = et;
                row.ez_energy = ez;
                row.ez_ratio = ratio_ez;
                row.residual_abs = residual.abs_norm;
                row.residual_rel = residual.rel_norm;
                row.P_fro = block_metrics.P_fro;
                row.Q_fro = block_metrics.Q_fro;
                row.P_tt_fro = block_metrics.P_tt_fro;
                row.P_zz_fro = block_metrics.P_zz_fro;
                row.Q_tt_fro = block_metrics.Q_tt_fro;
                row.Q_tz_fro = block_metrics.Q_tz_fro;
                row.Q_zt_fro = block_metrics.Q_zt_fro;
                row.Q_zz_fro = block_metrics.Q_zz_fro;
                row.P_tt_asym_rel = block_metrics.P_tt_asym_rel;
                row.P_zz_asym_rel = block_metrics.P_zz_asym_rel;
                row.Q_tt_asym_rel = block_metrics.Q_tt_asym_rel;
                row.Q_zz_asym_rel = block_metrics.Q_zz_asym_rel;
                row.Q_tz_Q_zt_transpose_mismatch_rel = block_metrics.Q_tz_Q_zt_transpose_mismatch_rel;
                row.block_norm_max = block_metrics.block_norm_max;
                row.P_tt_scale_rel = block_metrics.P_tt_scale_rel;
                row.P_zz_scale_rel = block_metrics.P_zz_scale_rel;
                row.Q_tt_scale_rel = block_metrics.Q_tt_scale_rel;
                row.Q_tz_scale_rel = block_metrics.Q_tz_scale_rel;
                row.Q_zt_scale_rel = block_metrics.Q_zt_scale_rel;
                row.Q_zz_scale_rel = block_metrics.Q_zz_scale_rel;
                result.matrix_audit.push_back(row);
                if (binding.original_index >= 0 && binding.original_index < n)
                    solver_to_matrix_row_index[(size_t)binding.original_index] =
                        static_cast<int>(result.matrix_audit.size()) - 1;
            }
        }
    }

    std::vector<RawRatioCandidate> raw;
    for (int i = 0; i < n; ++i)
    {
        const int raw_row_index =
            (i >= 0 && i < (int)solver_to_raw_row_index.size()) ? solver_to_raw_row_index[(size_t)i] : -1;
        auto *raw_row = (raw_row_index >= 0) ? &result.raw_spectrum[(size_t)raw_row_index] : nullptr;
        const int matrix_row_index =
            (i >= 0 && i < (int)solver_to_matrix_row_index.size()) ? solver_to_matrix_row_index[(size_t)i] : -1;
        auto *matrix_row =
            (matrix_row_index >= 0) ? &result.matrix_audit[(size_t)matrix_row_index] : nullptr;
        if (!std::isfinite(result.eig.lambda_re[i]))
        {
            if (raw_row != nullptr)
                raw_row->filter_reason = "discard_lambda_real_non_finite";
            if (matrix_row != nullptr)
                matrix_row->filter_reason = "discard_lambda_real_non_finite";
            continue;
        }
        if (std::abs(result.eig.lambda_im[i]) > 1e-7)
        {
            if (raw_row != nullptr)
                raw_row->filter_reason = "discard_lambda_imag_above_tol";
            if (matrix_row != nullptr)
                matrix_row->filter_reason = "discard_lambda_imag_above_tol";
            continue;
        }
        if (result.eig.lambda_re[i] <= 1e-10)
        {
            if (raw_row != nullptr)
                raw_row->filter_reason = "discard_lambda_real_non_positive";
            if (matrix_row != nullptr)
                matrix_row->filter_reason = "discard_lambda_real_non_positive";
            continue;
        }

        const double beta = std::sqrt(result.eig.lambda_re[i]);
        const double ratio = beta / k0;
        if (ratio <= 1e-6)
        {
            if (raw_row != nullptr)
                raw_row->filter_reason = "discard_beta_ratio_too_small";
            if (matrix_row != nullptr)
                matrix_row->filter_reason = "discard_beta_ratio_too_small";
            continue;
        }
        if (ratio > ratio_max)
        {
            if (raw_row != nullptr)
                raw_row->filter_reason = "discard_beta_ratio_above_physical_max";
            if (matrix_row != nullptr)
                matrix_row->filter_reason = "discard_beta_ratio_above_physical_max";
            continue;
        }
        if (raw_row != nullptr)
        {
            raw_row->filter_reason = "kept_physical_pre_dedup";
            raw_row->kept_after_physical_filter = 1;
        }
        if (matrix_row != nullptr)
        {
            matrix_row->filter_reason = "kept_physical_pre_dedup";
            matrix_row->kept_after_physical_filter = 1;
        }
        double ratio_ez = 0.0;
        if (raw_row != nullptr)
        {
            ratio_ez = raw_row->ez_ratio;
        }
        else if (matrix_row != nullptr)
        {
            ratio_ez = matrix_row->ez_ratio;
        }
        else
        {
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
            ratio_ez = (den > 0.0) ? (ez / den) : 0.0;
        }
        raw.push_back({i, raw_row_index >= 0 ? raw_row_index : matrix_row_index, ratio, ratio_ez});
    }
    std::sort(raw.begin(), raw.end(), [](const RawRatioCandidate &a, const RawRatioCandidate &b)
              { return a.beta_ratio < b.beta_ratio; });

    for (const RawRatioCandidate &cand : raw)
    {
        if (!result.candidates.empty() &&
            std::abs(cand.beta_ratio - result.candidates.back().beta_ratio) < 1e-8)
        {
            if (cand.raw_row_index >= 0)
                result.raw_spectrum[(size_t)cand.raw_row_index].filter_reason = "discard_dedup_nonrepresentative";
            if (cand.raw_row_index >= 0 && cand.raw_row_index < (int)result.matrix_audit.size())
                result.matrix_audit[(size_t)cand.raw_row_index].filter_reason = "discard_dedup_nonrepresentative";
            if (cand.ez_ratio > result.candidates.back().ez_ratio)
            {
                const int prev_row_index = result.candidates.back().raw_row_index;
                if (prev_row_index >= 0)
                {
                    if (prev_row_index < (int)result.raw_spectrum.size())
                    {
                        auto &prev = result.raw_spectrum[(size_t)prev_row_index];
                        prev.kept_after_dedup = 0;
                        prev.candidate_rank_after_dedup = 0;
                        prev.filter_reason = "discard_dedup_replaced_by_higher_ez";
                    }
                    if (prev_row_index < (int)result.matrix_audit.size())
                    {
                        auto &prev = result.matrix_audit[(size_t)prev_row_index];
                        prev.kept_after_dedup = 0;
                        prev.candidate_rank_after_dedup = 0;
                        prev.filter_reason = "discard_dedup_replaced_by_higher_ez";
                    }
                }
                result.candidates.back().eig_index = cand.eig_index;
                result.candidates.back().raw_row_index = cand.raw_row_index;
                result.candidates.back().ez_ratio = cand.ez_ratio;
                if (cand.raw_row_index >= 0)
                {
                    if (cand.raw_row_index < (int)result.raw_spectrum.size())
                    {
                        auto &kept = result.raw_spectrum[(size_t)cand.raw_row_index];
                        kept.kept_after_dedup = 1;
                        kept.candidate_rank_after_dedup = result.candidates.back().candidate_rank;
                        kept.filter_reason = "kept_physical_dedup_representative";
                    }
                    if (cand.raw_row_index < (int)result.matrix_audit.size())
                    {
                        auto &kept = result.matrix_audit[(size_t)cand.raw_row_index];
                        kept.kept_after_dedup = 1;
                        kept.candidate_rank_after_dedup = result.candidates.back().candidate_rank;
                        kept.filter_reason = "kept_physical_dedup_representative";
                    }
                }
            }
            continue;
        }

        RatioCandidate rep;
        rep.candidate_rank = static_cast<int>(result.candidates.size()) + 1;
        rep.eig_index = cand.eig_index;
        rep.raw_row_index = cand.raw_row_index;
        rep.beta_ratio = cand.beta_ratio;
        rep.ez_ratio = cand.ez_ratio;
        result.candidates.push_back(rep);
        if (cand.raw_row_index >= 0)
        {
            if (cand.raw_row_index < (int)result.raw_spectrum.size())
            {
                auto &kept = result.raw_spectrum[(size_t)cand.raw_row_index];
                kept.kept_after_dedup = 1;
                kept.candidate_rank_after_dedup = rep.candidate_rank;
                kept.filter_reason = "kept_physical_unique";
            }
            if (cand.raw_row_index < (int)result.matrix_audit.size())
            {
                auto &kept = result.matrix_audit[(size_t)cand.raw_row_index];
                kept.kept_after_dedup = 1;
                kept.candidate_rank_after_dedup = rep.candidate_rank;
                kept.filter_reason = "kept_physical_unique";
            }
        }
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
    CoupledBetaDiagVariant diag_variant,
    timing::Breakdown *perf = nullptr,
    const std::filesystem::path *linop_dir = nullptr,
    const std::string &prefix_base = "",
    const helmvec3_output::CaseDirs *out_dirs = nullptr,
    bool export_raw_spectrum_diag = false,
    bool export_matrix_audit_diag = false,
    double raw_spectrum_d_over_a = std::numeric_limits<double>::quiet_NaN())
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
        const bool should_export_diag =
            (export_raw_spectrum_diag || export_matrix_audit_diag) &&
            out_dirs != nullptr &&
            !point_prefix.empty() &&
            std::isfinite(raw_spectrum_d_over_a) &&
            should_export_raw_spectrum_diag(raw_spectrum_d_over_a, s);
        auto solve = solve_beta_point(
            mesh,
            eps,
            mu,
            k0,
            backend,
            diag_variant,
            perf,
            linop_dir,
            prefix_base.empty() ? nullptr : &point_prefix,
            should_export_diag,
            should_export_diag && export_matrix_audit_diag);
        const std::string priority_point =
            should_export_diag ? figure13_priority_point_key(raw_spectrum_d_over_a, s) : std::string();
        if (solve.candidates.empty())
        {
            if (should_export_diag)
            {
                std::vector<helmvec3_output::RawSpectrumCsvRecord> rows = solve.raw_spectrum;
                for (auto &row : rows)
                {
                    row.priority_point = priority_point;
                    row.d_over_a = raw_spectrum_d_over_a;
                    row.br_over_lambda0 = s;
                    row.match_target_beta_over_k0 = ref_ratio[i];
                }
                const auto raw_csv_path = out_dirs->csv / (point_prefix + "_raw_spectrum.csv");
                if (!helmvec3_output::write_raw_spectrum_csv(raw_csv_path, rows))
                    throw std::runtime_error("Falha ao escrever CSV de espectro bruto em " + raw_csv_path.string());
                if (export_matrix_audit_diag)
                {
                    std::vector<helmvec3_output::MatrixAuditCsvRecord> matrix_rows = solve.matrix_audit;
                    for (auto &row : matrix_rows)
                    {
                        row.priority_point = priority_point;
                        row.d_over_a = raw_spectrum_d_over_a;
                        row.br_over_lambda0 = s;
                        row.match_target_beta_over_k0 = ref_ratio[i];
                    }
                    const auto matrix_csv_path = out_dirs->csv / (point_prefix + "_matrix_audit.csv");
                    if (!helmvec3_output::write_matrix_audit_csv(matrix_csv_path, matrix_rows))
                        throw std::runtime_error("Falha ao escrever CSV de auditoria matricial em " + matrix_csv_path.string());
                }
            }
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
        if (should_export_diag)
        {
            mark_selected_raw_spectrum(solve.raw_spectrum, selected.selected_eig_index);
            std::vector<helmvec3_output::RawSpectrumCsvRecord> rows = solve.raw_spectrum;
            for (auto &row : rows)
            {
                row.priority_point = priority_point;
                row.d_over_a = raw_spectrum_d_over_a;
                row.br_over_lambda0 = s;
                row.match_target_beta_over_k0 = ref_ratio[i];
            }
            const auto raw_csv_path = out_dirs->csv / (point_prefix + "_raw_spectrum.csv");
            if (!helmvec3_output::write_raw_spectrum_csv(raw_csv_path, rows))
                throw std::runtime_error("Falha ao escrever CSV de espectro bruto em " + raw_csv_path.string());
            if (export_matrix_audit_diag)
            {
                mark_selected_matrix_audit(solve.matrix_audit, selected.selected_eig_index);
                std::vector<helmvec3_output::MatrixAuditCsvRecord> matrix_rows = solve.matrix_audit;
                for (auto &row : matrix_rows)
                {
                    row.priority_point = priority_point;
                    row.d_over_a = raw_spectrum_d_over_a;
                    row.br_over_lambda0 = s;
                    row.match_target_beta_over_k0 = ref_ratio[i];
                }
                const auto matrix_csv_path = out_dirs->csv / (point_prefix + "_matrix_audit.csv");
                if (!helmvec3_output::write_matrix_audit_csv(matrix_csv_path, matrix_rows))
                    throw std::runtime_error("Falha ao escrever CSV de auditoria matricial em " + matrix_csv_path.string());
            }
        }
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
    CoupledBetaDiagVariant diag_variant,
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
            diag_variant,
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
    bool used_legacy_positionals = false;
    bool used_legacy_debug = false;
};

struct Fig13CliConfig
{
    helmvec2::CoupledCliOptions cli;
    double d13_over_a = 0.20;
    int nx = 10;
    int ny = 5;
    bool used_legacy_positionals = false;
    bool used_legacy_debug = false;
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

    const bool has_named_primary_args =
        cfg.cli.nx_was_provided ||
        cfg.cli.ny_was_provided;

    if (cfg.cli.beta_was_provided || cfg.cli.d_over_a_preview_was_provided)
    {
        throw std::runtime_error(
            "helmvec3_fig12_rect aceita apenas --nx/--ny como aliases nomeados principais");
    }

    if (has_named_primary_args && !cfg.cli.positionals.empty())
    {
        throw std::runtime_error(
            "nao misture aliases nomeados principais com argumentos posicionais principais em helmvec3_fig12_rect; escolha apenas um estilo de chamada");
    }

    if (has_named_primary_args)
    {
        if (cfg.cli.nx_was_provided)
            cfg.nx = cfg.cli.nx;
        if (cfg.cli.ny_was_provided)
            cfg.ny = cfg.cli.ny;
    }
    else
    {
        if (cfg.cli.positionals.size() > 3)
        {
            throw std::runtime_error(
                "argumentos posicionais em excesso em helmvec3_fig12_rect; use no maximo [nx [ny [debug]]]");
        }
        cfg.used_legacy_positionals = !cfg.cli.positionals.empty();
        try
        {
            if (!cfg.cli.positionals.empty())
                cfg.nx = helmvec2::parse_positive_coupled_cli_int(cfg.cli.positionals[0], "nx");
            if (cfg.cli.positionals.size() >= 2)
                cfg.ny = helmvec2::parse_positive_coupled_cli_int(cfg.cli.positionals[1], "ny");
            if (cfg.cli.positionals.size() >= 3)
            {
                cfg.used_legacy_debug = true;
                const bool legacy_debug =
                    (helmvec2::parse_coupled_cli_int(cfg.cli.positionals[2], "debug") != 0);
                if (legacy_debug)
                {
                    cfg.cli.debug_local_blocks = true;
                    cfg.cli.debug_candidates = true;
                }
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                std::string("erro nos argumentos posicionais de helmvec3_fig12_rect: ") + e.what());
        }
    }

    return cfg;
}

Fig13CliConfig parse_fig13_cli(int argc, char **argv)
{
    Fig13CliConfig cfg;
    cfg.cli = helmvec2::parse_coupled_cli_options(argc, argv);

    const bool has_named_primary_args =
        cfg.cli.d_over_a_preview_was_provided ||
        cfg.cli.nx_was_provided ||
        cfg.cli.ny_was_provided;

    if (cfg.cli.beta_was_provided)
    {
        throw std::runtime_error(
            "helmvec3_fig13_rect nao aceita --beta; use --d-over-a-preview");
    }

    if (has_named_primary_args && !cfg.cli.positionals.empty())
    {
        throw std::runtime_error(
            "nao misture aliases nomeados principais com argumentos posicionais principais em helmvec3_fig13_rect; escolha apenas um estilo de chamada");
    }

    if (has_named_primary_args)
    {
        if (cfg.cli.d_over_a_preview_was_provided)
            cfg.d13_over_a = cfg.cli.d_over_a_preview;
        if (cfg.cli.nx_was_provided)
            cfg.nx = cfg.cli.nx;
        if (cfg.cli.ny_was_provided)
            cfg.ny = cfg.cli.ny;
    }
    else
    {
        if (cfg.cli.positionals.size() > 4)
        {
            throw std::runtime_error(
                "argumentos posicionais em excesso em helmvec3_fig13_rect; use no maximo [d_over_a_preview [nx [ny [debug]]]]");
        }
        cfg.used_legacy_positionals = !cfg.cli.positionals.empty();
        try
        {
            if (!cfg.cli.positionals.empty())
                cfg.d13_over_a = helmvec2::parse_cli_real(cfg.cli.positionals[0], "d_over_a_preview");
            if (cfg.cli.positionals.size() >= 2)
                cfg.nx = helmvec2::parse_positive_coupled_cli_int(cfg.cli.positionals[1], "nx");
            if (cfg.cli.positionals.size() >= 3)
                cfg.ny = helmvec2::parse_positive_coupled_cli_int(cfg.cli.positionals[2], "ny");
            if (cfg.cli.positionals.size() >= 4)
            {
                cfg.used_legacy_debug = true;
                const bool legacy_debug =
                    (helmvec2::parse_coupled_cli_int(cfg.cli.positionals[3], "debug") != 0);
                if (legacy_debug)
                {
                    cfg.cli.debug_local_blocks = true;
                    cfg.cli.debug_candidates = true;
                }
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                std::string("erro nos argumentos posicionais de helmvec3_fig13_rect: ") + e.what());
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
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso: ./helmvec3_fig12_rect [nx [ny [debug]]]"
           << " [--backend closed-form|gauss]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--nx NX] [--ny NY]"
           << " (nao misture com os posicionais principais)\n";
        os << "Compatibilidade: os posicionais principais continuam aceitos, mas estao deprecated;"
           << " prefira os aliases nomeados acima.\n";
        os << "Compatibilidade: o debug posicional legado [debug] continua aceito, mas esta deprecated;"
           << " prefira --debug-local-blocks e/ou --debug-candidates.\n";
    };

    Fig12CliConfig cfg;
    if (helmvec2::coupled_cli_requests_help(argc, argv))
    {
        print_usage(std::cout);
        return 0;
    }
    try
    {
        cfg = parse_fig12_cli(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        print_usage(std::cerr);
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
    if (cfg.used_legacy_positionals)
    {
        std::cerr << "Aviso: os argumentos posicionais principais de helmvec3_fig12_rect continuam aceitos por compatibilidade, mas estao deprecated; prefira --nx/--ny.\n";
    }
    if (cfg.used_legacy_debug)
    {
        std::cerr << "Aviso: o debug posicional legado de helmvec3_fig12_rect continua aceito por compatibilidade, mas esta deprecated; prefira --debug-local-blocks e/ou --debug-candidates.\n";
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
        CoupledBetaDiagVariant::Baseline,
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
    const auto print_usage = [](std::ostream &os)
    {
        os << "Uso: ./helmvec3_fig13_rect [d_over_a_preview [nx [ny [debug]]]]"
           << " [--backend closed-form|gauss]"
           << " [--debug-local-blocks] [--debug-candidates]\n";
        os << "Aliases nomeados: [--d-over-a-preview VAL] [--nx NX] [--ny NY]"
           << " (nao misture com os posicionais principais)\n";
        os << "Compatibilidade: os posicionais principais continuam aceitos, mas estao deprecated;"
           << " prefira os aliases nomeados acima.\n";
        os << "Compatibilidade: o debug posicional legado [debug] continua aceito, mas esta deprecated;"
           << " prefira --debug-local-blocks e/ou --debug-candidates.\n";
    };

    Fig13CliConfig cfg;
    if (helmvec2::coupled_cli_requests_help(argc, argv))
    {
        print_usage(std::cout);
        return 0;
    }
    try
    {
        cfg = parse_fig13_cli(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    }

    const auto out_dirs = helmvec3_output::ensure_case_dirs("fig13_rect");
    const bool export_raw_spectrum_diag = env_flag_enabled("TP3485_HELMVEC3_DIAG_RAW_SPECTRUM");
    const bool export_matrix_audit_diag = env_flag_enabled("TP3485_HELMVEC3_DIAG_MATRIX_AUDIT");
    const CoupledBetaDiagVariant diag_variant = env_beta_diag_variant();
    execution_log::ExecutionLogScope log_scope((out_dirs.root / "run.log").string());
    if (!log_scope.active())
    {
        std::cerr << "Aviso: nao foi possivel abrir run.log em "
                  << (out_dirs.root / "run.log")
                  << " (" << log_scope.error_message() << ")\n";
    }
    if (cfg.used_legacy_positionals)
    {
        std::cerr << "Aviso: os argumentos posicionais principais de helmvec3_fig13_rect continuam aceitos por compatibilidade, mas estao deprecated; prefira --d-over-a-preview/--nx/--ny.\n";
    }
    if (cfg.used_legacy_debug)
    {
        std::cerr << "Aviso: o debug posicional legado de helmvec3_fig13_rect continua aceito por compatibilidade, mas esta deprecated; prefira --debug-local-blocks e/ou --debug-candidates.\n";
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
        diag_variant,
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
            diag_variant,
            &perf,
            &out_dirs.linop,
            std::string("helmvec3_fig13_rect_table10_da") + label_number(blk.d_over_a),
            &out_dirs,
            export_raw_spectrum_diag,
            export_matrix_audit_diag,
            blk.d_over_a);

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
