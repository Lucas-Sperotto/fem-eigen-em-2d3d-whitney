/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_case_output.hpp                                */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de saida para a familia HELMVEC1.         */
/*****************************************************************************/

#pragma once

#include "core/output_paths.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace helmvec1_output
{

struct CaseDirs
{
    std::filesystem::path root;
    std::filesystem::path csv;
    std::filesystem::path img;
};

struct ModeCsvRecord
{
    std::string formulation;
    std::string dominant_block;
    std::string component_label;
    std::string family;
    std::string mode_label;
    int positive_rank = 0;
    int eig_index = 0;
    std::string m;
    std::string n;
    std::string p;
    std::string ar_m;
    std::string b_m;
    std::string r_m;
    std::string r1_m;
    std::string r2_m;
    std::string kc_fem;
    std::string kc_ana;
    std::string kc_ar_fem;
    std::string kc_ar_ana;
    std::string kc_r_fem;
    std::string kc_r_ana;
    std::string kc_r1_fem;
    std::string kc_r1_ana;
    std::string error_percent;
    double rho_abs = 0.0;
    double edge_energy = 0.0;
    double scalar_energy = 0.0;
    double dominant_energy_ratio = 0.0;
    std::string match_space;
    std::string match_method;
    std::string mode_status;
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
    dirs.root = output_paths::ensure_case_dir("helmvec1/" + case_name);
    dirs.csv = dirs.root / "csv";
    dirs.img = dirs.root / "img";
    std::error_code ec;
    std::filesystem::create_directories(dirs.csv, ec);
    std::filesystem::create_directories(dirs.img, ec);
    return dirs;
}

inline bool write_modes_csv(
    const std::filesystem::path &path,
    const std::vector<ModeCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "formulation,dominant_block,component_label,family,mode_label,"
           "positive_rank,eig_index,m,n,p,ar_m,b_m,r_m,r1_m,r2_m,"
           "kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,kc_r_fem,kc_r_ana,"
           "kc_r1_fem,kc_r1_ana,error_percent,rho_abs,edge_energy,scalar_energy,"
           "dominant_energy_ratio,match_space,match_method,mode_status\n";

    for (const ModeCsvRecord &rec : records)
    {
        out << rec.formulation << ","
            << rec.dominant_block << ","
            << rec.component_label << ","
            << rec.family << ","
            << rec.mode_label << ","
            << rec.positive_rank << ","
            << rec.eig_index << ","
            << rec.m << ","
            << rec.n << ","
            << rec.p << ","
            << rec.ar_m << ","
            << rec.b_m << ","
            << rec.r_m << ","
            << rec.r1_m << ","
            << rec.r2_m << ","
            << rec.kc_fem << ","
            << rec.kc_ana << ","
            << rec.kc_ar_fem << ","
            << rec.kc_ar_ana << ","
            << rec.kc_r_fem << ","
            << rec.kc_r_ana << ","
            << rec.kc_r1_fem << ","
            << rec.kc_r1_ana << ","
            << rec.error_percent << ","
            << rec.rho_abs << ","
            << rec.edge_energy << ","
            << rec.scalar_energy << ","
            << rec.dominant_energy_ratio << ","
            << rec.match_space << ","
            << rec.match_method << ","
            << rec.mode_status << "\n";
    }

    return static_cast<bool>(out);
}

} // namespace helmvec1_output
