/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/run_timing_dispersion_csv.hpp                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Escrita estruturada dos tempos de execucao dos casos HELMVEC3.  */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iomanip>
#include <string>

namespace run_timing_dispersion_csv
{

struct Record
{
    std::string case_label;
    std::string geometry;
    std::string backend;
    int debug_local_blocks = 0;
    int debug_candidates = 0;
    double a = 0.0;
    double b = 0.0;
    double d12 = 0.0;
    double d13_preview_over_a = 0.0;
    double eps_fill = 0.0;
    int nx = 0;
    int ny = 0;
    int mesh_nodes = 0;
    int mesh_tris = 0;
    int table9_sample_count = 0;
    int preview_sample_count = 0;
    int table10_block_count = 0;
    int table10_row_count = 0;
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
           "a,b,d12,d13_preview_over_a,eps_fill,nx,ny,mesh_nodes,mesh_tris,"
           "table9_sample_count,preview_sample_count,table10_block_count,"
           "table10_row_count,assembly_ms,solve_ms,post_ms,total_ms\n";

    out << record.case_label << ","
        << record.geometry << ","
        << record.backend << ","
        << record.debug_local_blocks << ","
        << record.debug_candidates << ","
        << record.a << ","
        << record.b << ","
        << record.d12 << ","
        << record.d13_preview_over_a << ","
        << record.eps_fill << ","
        << record.nx << ","
        << record.ny << ","
        << record.mesh_nodes << ","
        << record.mesh_tris << ","
        << record.table9_sample_count << ","
        << record.preview_sample_count << ","
        << record.table10_block_count << ","
        << record.table10_row_count << ","
        << record.assembly_ms << ","
        << record.solve_ms << ","
        << record.post_ms << ","
        << record.total_ms << "\n";

    return static_cast<bool>(out);
}

} // namespace run_timing_dispersion_csv
