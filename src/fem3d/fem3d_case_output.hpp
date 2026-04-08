/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_case_output.hpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de saida para os casos 3D FEM3D.          */
/*****************************************************************************/

#pragma once

#include "core/error_metrics.hpp"
#include "core/output_paths.hpp"
#include "fem3d_compare.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace fem3d_output
{

struct SolverDirs
{
    std::filesystem::path root;
};

struct CaseDirs
{
    std::filesystem::path solver_root;
    std::filesystem::path root;
    std::filesystem::path csv;
    std::filesystem::path linop;
    std::filesystem::path vtk;
};

struct ModeCsvRecord
{
    int reference_index = 0;
    std::string mode_label;
    double k0_analytic = 0.0;
    double k0_fem = std::numeric_limits<double>::quiet_NaN();
    double ref_paper = 0.0;
    double error_percent_analytic = std::numeric_limits<double>::quiet_NaN();
    double error_percent_ref_paper = std::numeric_limits<double>::quiet_NaN();
    int matched_eig_index = -1;
    std::string match_status;
    std::string field_status;
    std::string fields_csv_file;
    std::string vtk_file;
};

inline SolverDirs ensure_solver_dirs(const std::string &solver_name)
{
    SolverDirs dirs;
    dirs.root = output_paths::ensure_case_dir(solver_name);
    return dirs;
}

inline CaseDirs ensure_case_dirs(const std::string &solver_name, const std::string &case_name)
{
    CaseDirs dirs;
    dirs.solver_root = output_paths::ensure_case_dir(solver_name);
    dirs.root = output_paths::ensure_case_dir(solver_name + "/" + case_name);
    dirs.csv = dirs.root / "csv";
    dirs.linop = dirs.root / "linop";
    dirs.vtk = dirs.root / "vtk";
    std::error_code ec;
    std::filesystem::create_directories(dirs.csv, ec);
    std::filesystem::create_directories(dirs.linop, ec);
    std::filesystem::create_directories(dirs.vtk, ec);
    return dirs;
}

inline std::vector<ModeCsvRecord> build_mode_records(
    const std::vector<fem3d::RefRow> &ref_rows,
    const std::vector<double> &k0_matched)
{
    std::vector<ModeCsvRecord> records;
    records.reserve(ref_rows.size());

    for (size_t i = 0; i < ref_rows.size(); ++i)
    {
        const double k0_fem =
            i < k0_matched.size() ? k0_matched[i] : std::numeric_limits<double>::quiet_NaN();

        ModeCsvRecord row;
        row.reference_index = static_cast<int>(i) + 1;
        row.mode_label = ref_rows[i].mode ? ref_rows[i].mode : "";
        row.k0_analytic = ref_rows[i].analytical;
        row.k0_fem = k0_fem;
        row.ref_paper = ref_rows[i].ref_paper;
        row.error_percent_analytic =
            error_metrics::absolute_relative_error_percent(ref_rows[i].analytical, k0_fem);
        row.error_percent_ref_paper =
            error_metrics::absolute_relative_error_percent(ref_rows[i].ref_paper, k0_fem);
        row.match_status = std::isfinite(k0_fem) ? "matched" : "missing";
        records.push_back(row);
    }

    return records;
}

inline bool write_modes_csv(
    const std::filesystem::path &path,
    const std::vector<ModeCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "reference_index,mode_label,k0_analytic,k0_fem,ref_paper,"
           "error_percent_analytic,error_percent_ref_paper,matched_eig_index,"
           "match_status,field_status,fields_csv_file,vtk_file\n";

    for (const ModeCsvRecord &row : records)
    {
        out << row.reference_index << ","
            << row.mode_label << ","
            << row.k0_analytic << ","
            << row.k0_fem << ","
            << row.ref_paper << ","
            << row.error_percent_analytic << ","
            << row.error_percent_ref_paper << ","
            << row.matched_eig_index << ","
            << row.match_status << ","
            << row.field_status << ","
            << row.fields_csv_file << ","
            << row.vtk_file << "\n";
    }

    return static_cast<bool>(out);
}

} // namespace fem3d_output
