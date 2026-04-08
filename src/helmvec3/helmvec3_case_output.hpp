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

} // namespace helmvec3_output
