/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec3/helmvec3_case_output.hpp                             */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de saida para a familia HELMVEC3.         */
/*****************************************************************************/

#pragma once

#include "core/output_paths.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace helmvec3_output
{

struct CaseDirs
{
    std::filesystem::path root;
    std::filesystem::path csv;
    std::filesystem::path vtk;
    std::filesystem::path img;
    std::filesystem::path linop;
};

struct Table9CsvRecord
{
    double br_over_lambda0 = 0.0;
    double beta_over_k0_fem = 0.0;
    double beta_over_k0_analytic = 0.0;
    double beta_over_k0_helmvec3 = 0.0;
    int selected_candidate_rank = 0;
    int selected_eig_index = -1;
    double ez_ratio = 0.0;
    double error_percent_analytic = 0.0;
    double error_percent_helmvec3 = 0.0;
    std::string match_status;
    std::string field_status;
    std::string et_fields_csv_file;
    std::string ez_fields_csv_file;
    std::string et_vtk_file;
    std::string ez_vtk_file;
};

struct PreviewCsvRecord
{
    double d_over_a_preview = 0.0;
    double br_over_lambda0 = 0.0;
    double beta_over_k0_fem_branch = 0.0;
    int selected_candidate_rank = 0;
    int selected_eig_index = -1;
    double ez_ratio = 0.0;
    std::string branch_status;
    std::string field_status;
    std::string et_fields_csv_file;
    std::string ez_fields_csv_file;
    std::string et_vtk_file;
    std::string ez_vtk_file;
};

struct Table10CsvRecord
{
    double d_over_a = 0.0;
    double br_over_lambda0 = 0.0;
    double beta_over_k0_fem_matched = 0.0;
    double beta_over_k0_analytic = 0.0;
    double beta_over_k0_helmvec3 = 0.0;
    int selected_candidate_rank = 0;
    int selected_eig_index = -1;
    double ez_ratio = 0.0;
    double error_percent_analytic = 0.0;
    double error_percent_helmvec3 = 0.0;
    std::string match_status;
    std::string field_status;
    std::string et_fields_csv_file;
    std::string ez_fields_csv_file;
    std::string et_vtk_file;
    std::string ez_vtk_file;
};

struct RawSpectrumCsvRecord
{
    std::string section;
    std::string case_label;
    std::string priority_point;
    double d_over_a = 0.0;
    double br_over_lambda0 = 0.0;
    double match_target_beta_over_k0 = 0.0;
    int ordered_rank = 0;
    int solver_index = -1;
    double lambda_real = 0.0;
    double lambda_imag = 0.0;
    double beta_ratio_if_real_positive = 0.0;
    std::string filter_reason;
    int kept_after_physical_filter = 0;
    int kept_after_dedup = 0;
    int candidate_rank_after_dedup = 0;
    int selected_after_matching = 0;
    double et_energy = 0.0;
    double ez_energy = 0.0;
    double ez_ratio = 0.0;
};

struct MatrixAuditCsvRecord
{
    std::string section;
    std::string case_label;
    std::string priority_point;
    std::string backend;
    double d_over_a = 0.0;
    double br_over_lambda0 = 0.0;
    double match_target_beta_over_k0 = 0.0;
    int ordered_rank = 0;
    int solver_index = -1;
    double lambda_real = 0.0;
    double lambda_imag = 0.0;
    double beta_ratio_if_real_positive = 0.0;
    std::string filter_reason;
    int kept_after_physical_filter = 0;
    int kept_after_dedup = 0;
    int candidate_rank_after_dedup = 0;
    int selected_after_matching = 0;
    double et_energy = 0.0;
    double ez_energy = 0.0;
    double ez_ratio = 0.0;
    double residual_abs = 0.0;
    double residual_rel = 0.0;
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

inline std::string format_float(double value)
{
    std::ostringstream oss;
    oss << std::setprecision(16) << value;
    return oss.str();
}

inline CaseDirs ensure_case_dirs(const std::string &case_name)
{
    CaseDirs dirs;
    dirs.root = output_paths::ensure_case_dir("helmvec3/" + case_name);
    dirs.csv = dirs.root / "csv";
    dirs.vtk = dirs.root / "vtk";
    dirs.img = dirs.root / "img";
    dirs.linop = dirs.root / "linop";
    std::error_code ec;
    std::filesystem::create_directories(dirs.csv, ec);
    std::filesystem::create_directories(dirs.vtk, ec);
    std::filesystem::create_directories(dirs.img, ec);
    std::filesystem::create_directories(dirs.linop, ec);
    return dirs;
}

inline bool write_table9_csv(
    const std::filesystem::path &path,
    const std::vector<Table9CsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "br_over_lambda0,beta_over_k0_fem,beta_over_k0_analytic,"
           "beta_over_k0_helmvec3,selected_candidate_rank,selected_eig_index,"
           "ez_ratio,error_percent_analytic,error_percent_helmvec3,match_status,"
           "field_status,et_fields_csv_file,ez_fields_csv_file,et_vtk_file,"
           "ez_vtk_file\n";

    for (const Table9CsvRecord &rec : records)
    {
        out << rec.br_over_lambda0 << ","
            << rec.beta_over_k0_fem << ","
            << rec.beta_over_k0_analytic << ","
            << rec.beta_over_k0_helmvec3 << ","
            << rec.selected_candidate_rank << ","
            << rec.selected_eig_index << ","
            << rec.ez_ratio << ","
            << rec.error_percent_analytic << ","
            << rec.error_percent_helmvec3 << ","
            << rec.match_status << ","
            << rec.field_status << ","
            << rec.et_fields_csv_file << ","
            << rec.ez_fields_csv_file << ","
            << rec.et_vtk_file << ","
            << rec.ez_vtk_file << "\n";
    }

    return static_cast<bool>(out);
}

inline bool write_preview_csv(
    const std::filesystem::path &path,
    const std::vector<PreviewCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "d_over_a_preview,br_over_lambda0,beta_over_k0_fem_branch,"
           "selected_candidate_rank,selected_eig_index,ez_ratio,branch_status,"
           "field_status,et_fields_csv_file,ez_fields_csv_file,et_vtk_file,"
           "ez_vtk_file\n";

    for (const PreviewCsvRecord &rec : records)
    {
        out << rec.d_over_a_preview << ","
            << rec.br_over_lambda0 << ","
            << rec.beta_over_k0_fem_branch << ","
            << rec.selected_candidate_rank << ","
            << rec.selected_eig_index << ","
            << rec.ez_ratio << ","
            << rec.branch_status << ","
            << rec.field_status << ","
            << rec.et_fields_csv_file << ","
            << rec.ez_fields_csv_file << ","
            << rec.et_vtk_file << ","
            << rec.ez_vtk_file << "\n";
    }

    return static_cast<bool>(out);
}

inline bool write_table10_csv(
    const std::filesystem::path &path,
    const std::vector<Table10CsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "d_over_a,br_over_lambda0,beta_over_k0_fem_matched,"
           "beta_over_k0_analytic,beta_over_k0_helmvec3,selected_candidate_rank,"
           "selected_eig_index,ez_ratio,error_percent_analytic,"
           "error_percent_helmvec3,match_status,field_status,et_fields_csv_file,"
           "ez_fields_csv_file,et_vtk_file,ez_vtk_file\n";

    for (const Table10CsvRecord &rec : records)
    {
        out << rec.d_over_a << ","
            << rec.br_over_lambda0 << ","
            << rec.beta_over_k0_fem_matched << ","
            << rec.beta_over_k0_analytic << ","
            << rec.beta_over_k0_helmvec3 << ","
            << rec.selected_candidate_rank << ","
            << rec.selected_eig_index << ","
            << rec.ez_ratio << ","
            << rec.error_percent_analytic << ","
            << rec.error_percent_helmvec3 << ","
            << rec.match_status << ","
            << rec.field_status << ","
            << rec.et_fields_csv_file << ","
            << rec.ez_fields_csv_file << ","
            << rec.et_vtk_file << ","
            << rec.ez_vtk_file << "\n";
    }

    return static_cast<bool>(out);
}

inline bool write_raw_spectrum_csv(
    const std::filesystem::path &path,
    const std::vector<RawSpectrumCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "section,case,priority_point,d_over_a,br_over_lambda0,"
           "match_target_beta_over_k0,ordered_rank,solver_index,lambda_real,"
           "lambda_imag,beta_ratio_if_real_positive,filter_reason,"
           "kept_after_physical_filter,kept_after_dedup,"
           "candidate_rank_after_dedup,selected_after_matching,et_energy,"
           "ez_energy,ez_ratio\n";

    for (const RawSpectrumCsvRecord &rec : records)
    {
        out << rec.section << ","
            << rec.case_label << ","
            << rec.priority_point << ","
            << rec.d_over_a << ","
            << rec.br_over_lambda0 << ","
            << rec.match_target_beta_over_k0 << ","
            << rec.ordered_rank << ","
            << rec.solver_index << ","
            << rec.lambda_real << ","
            << rec.lambda_imag << ","
            << rec.beta_ratio_if_real_positive << ","
            << rec.filter_reason << ","
            << rec.kept_after_physical_filter << ","
            << rec.kept_after_dedup << ","
            << rec.candidate_rank_after_dedup << ","
            << rec.selected_after_matching << ","
            << rec.et_energy << ","
            << rec.ez_energy << ","
            << rec.ez_ratio << "\n";
    }

    return static_cast<bool>(out);
}

inline bool write_matrix_audit_csv(
    const std::filesystem::path &path,
    const std::vector<MatrixAuditCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "section,case,priority_point,backend,d_over_a,br_over_lambda0,"
           "match_target_beta_over_k0,ordered_rank,solver_index,lambda_real,"
           "lambda_imag,beta_ratio_if_real_positive,filter_reason,"
           "kept_after_physical_filter,kept_after_dedup,"
           "candidate_rank_after_dedup,selected_after_matching,et_energy,"
           "ez_energy,ez_ratio,residual_abs,residual_rel,P_fro,Q_fro,"
           "P_tt_fro,P_zz_fro,Q_tt_fro,Q_tz_fro,Q_zt_fro,Q_zz_fro,"
           "P_tt_asym_rel,P_zz_asym_rel,Q_tt_asym_rel,Q_zz_asym_rel,"
           "Q_tz_Q_zt_transpose_mismatch_rel,block_norm_max,P_tt_scale_rel,"
           "P_zz_scale_rel,Q_tt_scale_rel,Q_tz_scale_rel,Q_zt_scale_rel,"
           "Q_zz_scale_rel\n";

    for (const MatrixAuditCsvRecord &rec : records)
    {
        out << rec.section << ","
            << rec.case_label << ","
            << rec.priority_point << ","
            << rec.backend << ","
            << rec.d_over_a << ","
            << rec.br_over_lambda0 << ","
            << rec.match_target_beta_over_k0 << ","
            << rec.ordered_rank << ","
            << rec.solver_index << ","
            << rec.lambda_real << ","
            << rec.lambda_imag << ","
            << rec.beta_ratio_if_real_positive << ","
            << rec.filter_reason << ","
            << rec.kept_after_physical_filter << ","
            << rec.kept_after_dedup << ","
            << rec.candidate_rank_after_dedup << ","
            << rec.selected_after_matching << ","
            << rec.et_energy << ","
            << rec.ez_energy << ","
            << rec.ez_ratio << ","
            << rec.residual_abs << ","
            << rec.residual_rel << ","
            << rec.P_fro << ","
            << rec.Q_fro << ","
            << rec.P_tt_fro << ","
            << rec.P_zz_fro << ","
            << rec.Q_tt_fro << ","
            << rec.Q_tz_fro << ","
            << rec.Q_zt_fro << ","
            << rec.Q_zz_fro << ","
            << rec.P_tt_asym_rel << ","
            << rec.P_zz_asym_rel << ","
            << rec.Q_tt_asym_rel << ","
            << rec.Q_zz_asym_rel << ","
            << rec.Q_tz_Q_zt_transpose_mismatch_rel << ","
            << rec.block_norm_max << ","
            << rec.P_tt_scale_rel << ","
            << rec.P_zz_scale_rel << ","
            << rec.Q_tt_scale_rel << ","
            << rec.Q_tz_scale_rel << ","
            << rec.Q_zt_scale_rel << ","
            << rec.Q_zz_scale_rel << "\n";
    }

    return static_cast<bool>(out);
}

} // namespace helmvec3_output
