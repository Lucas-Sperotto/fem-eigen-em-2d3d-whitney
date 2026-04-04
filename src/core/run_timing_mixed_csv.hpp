/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/run_timing_mixed_csv.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Escrita estruturada dos tempos de execucao dos casos HELMVEC1.  */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iomanip>
#include <string>

namespace run_timing_mixed_csv
{

struct Record
{
    std::string case_label;
    std::string geometry;
    std::string backend;
    int debug_local_blocks = 0;
    int debug_candidates = 0;
    std::string ar_m;
    std::string b_m;
    std::string r_m;
    std::string r1_m;
    std::string r2_m;
    std::string nx;
    std::string ny;
    std::string nr;
    std::string nt;
    int mesh_nodes = 0;
    int mesh_tris = 0;
    double e_assembly_ms = 0.0;
    double e_solve_ms = 0.0;
    double e_post_ms = 0.0;
    double h_assembly_ms = 0.0;
    double h_solve_ms = 0.0;
    double h_post_ms = 0.0;
    double assembly_ms_total = 0.0;
    double solve_ms_total = 0.0;
    double post_ms_total = 0.0;
    double total_ms = 0.0;
};

inline bool write_csv(const std::string &path, const Record &record)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "case_label,geometry,backend,debug_local_blocks,debug_candidates,"
           "ar_m,b_m,r_m,r1_m,r2_m,nx,ny,nr,nt,mesh_nodes,mesh_tris,"
           "e_assembly_ms,e_solve_ms,e_post_ms,"
           "h_assembly_ms,h_solve_ms,h_post_ms,"
           "assembly_ms_total,solve_ms_total,post_ms_total,total_ms\n";

    out << record.case_label << ","
        << record.geometry << ","
        << record.backend << ","
        << record.debug_local_blocks << ","
        << record.debug_candidates << ","
        << record.ar_m << ","
        << record.b_m << ","
        << record.r_m << ","
        << record.r1_m << ","
        << record.r2_m << ","
        << record.nx << ","
        << record.ny << ","
        << record.nr << ","
        << record.nt << ","
        << record.mesh_nodes << ","
        << record.mesh_tris << ","
        << record.e_assembly_ms << ","
        << record.e_solve_ms << ","
        << record.e_post_ms << ","
        << record.h_assembly_ms << ","
        << record.h_solve_ms << ","
        << record.h_post_ms << ","
        << record.assembly_ms_total << ","
        << record.solve_ms_total << ","
        << record.post_ms_total << ","
        << record.total_ms << "\n";

    return static_cast<bool>(out);
}

} // namespace run_timing_mixed_csv
