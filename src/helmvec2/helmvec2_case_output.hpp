/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/helmvec2_case_output.hpp                             */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de saida para a familia HELMVEC2.         */
/*****************************************************************************/

#pragma once

#include "core/output_paths.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace helmvec2_output
{

struct CaseDirs
{
    std::filesystem::path root;
    std::filesystem::path csv;
    std::filesystem::path img;
};

struct ModeCsvRecord
{
    int mode = 0;
    int matched_candidate_rank = 0;
    double k0_fem_matched = 0.0;
    double k0l_fem_matched = 0.0;
    double ref_helmvec2 = 0.0;
    double ref_hayata = 0.0;
    double error_percent_helmvec2 = 0.0;
    double error_percent_hayata = 0.0;
    std::string match_status;
};

struct CandidateCsvRecord
{
    int candidate_rank = 0;
    int eig_index = 0;
    double k0 = 0.0;
    double k0l = 0.0;
    double ez_ratio = 0.0;
    int passes_physical_filter = 0;
};

inline CaseDirs ensure_case_dirs(const std::string &case_name)
{
    CaseDirs dirs;
    dirs.root = output_paths::ensure_case_dir("helmvec2/" + case_name);
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
    out << "mode,matched_candidate_rank,k0_fem_matched,k0l_fem_matched,"
           "ref_helmvec2,ref_hayata,error_percent_helmvec2,error_percent_hayata,"
           "match_status\n";

    for (const ModeCsvRecord &rec : records)
    {
        out << rec.mode << ","
            << rec.matched_candidate_rank << ","
            << rec.k0_fem_matched << ","
            << rec.k0l_fem_matched << ","
            << rec.ref_helmvec2 << ","
            << rec.ref_hayata << ","
            << rec.error_percent_helmvec2 << ","
            << rec.error_percent_hayata << ","
            << rec.match_status << "\n";
    }

    return static_cast<bool>(out);
}

inline bool write_candidates_csv(
    const std::filesystem::path &path,
    const std::vector<CandidateCsvRecord> &records)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "candidate_rank,eig_index,k0,k0l,ez_ratio,passes_physical_filter\n";

    for (const CandidateCsvRecord &rec : records)
    {
        out << rec.candidate_rank << ","
            << rec.eig_index << ","
            << rec.k0 << ","
            << rec.k0l << ","
            << rec.ez_ratio << ","
            << rec.passes_physical_filter << "\n";
    }

    return static_cast<bool>(out);
}

} // namespace helmvec2_output
