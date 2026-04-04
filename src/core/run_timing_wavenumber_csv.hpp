/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/run_timing_wavenumber_csv.hpp                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Escrita estruturada dos tempos de execucao dos casos HELMVEC2.  */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iomanip>
#include <string>

namespace run_timing_wavenumber_csv
{

struct Record
{
    std::string case_label;
    std::string geometry;
    std::string backend;
    int debug_local_blocks = 0;
    int debug_candidates = 0;
    double beta = 0.0;
    double L = 0.0;
    double betaL = 0.0;
    int nx = 0;
    int ny = 0;
    int mesh_nodes = 0;
    int mesh_tris = 0;
    double eps_top = 0.0;
    double eps_bottom = 0.0;
    double mu_r = 1.0;
    double k0_min_phys = 0.0;
    int candidate_count_raw = 0;
    int candidate_count_phys = 0;
    int matched_mode_count = 0;
    double assembly_ms = 0.0;
    double solve_ms = 0.0;
    double post_ms = 0.0;
    double total_ms = 0.0;
};

inline bool write_csv(const std::string &path, const Record &record)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "case_label,geometry,backend,debug_local_blocks,debug_candidates,"
           "beta,L,betaL,nx,ny,mesh_nodes,mesh_tris,eps_top,eps_bottom,mu_r,"
           "k0_min_phys,candidate_count_raw,candidate_count_phys,matched_mode_count,"
           "assembly_ms,solve_ms,post_ms,total_ms\n";

    out << record.case_label << ","
        << record.geometry << ","
        << record.backend << ","
        << record.debug_local_blocks << ","
        << record.debug_candidates << ","
        << record.beta << ","
        << record.L << ","
        << record.betaL << ","
        << record.nx << ","
        << record.ny << ","
        << record.mesh_nodes << ","
        << record.mesh_tris << ","
        << record.eps_top << ","
        << record.eps_bottom << ","
        << record.mu_r << ","
        << record.k0_min_phys << ","
        << record.candidate_count_raw << ","
        << record.candidate_count_phys << ","
        << record.matched_mode_count << ","
        << record.assembly_ms << ","
        << record.solve_ms << ","
        << record.post_ms << ","
        << record.total_ms << "\n";

    return static_cast<bool>(out);
}

} // namespace run_timing_wavenumber_csv
