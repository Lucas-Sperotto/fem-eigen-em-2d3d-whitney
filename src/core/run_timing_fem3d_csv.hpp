/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/run_timing_fem3d_csv.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Escrita estruturada dos tempos de execucao dos casos FEM3D.     */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iomanip>
#include <string>

namespace run_timing_fem3d_csv
{

struct Record
{
    std::string solver;
    std::string case_name;
    std::string storage;
    std::string backend;
    int debug_local_blocks = 0;
    int debug_candidates = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int mesh_nodes = 0;
    int mesh_tets = 0;
    int mesh_edges = 0;
    int ndof = 0;
    unsigned long long nnz_lower_s = 0;
    unsigned long long nnz_lower_t = 0;
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
    out << "solver,case_name,storage,backend,debug_local_blocks,debug_candidates,"
           "nx,ny,nz,mesh_nodes,mesh_tets,mesh_edges,ndof,nnz_lower_s,nnz_lower_t,"
           "assembly_ms,solve_ms,post_ms,total_ms\n";

    out << record.solver << ","
        << record.case_name << ","
        << record.storage << ","
        << record.backend << ","
        << record.debug_local_blocks << ","
        << record.debug_candidates << ","
        << record.nx << ","
        << record.ny << ","
        << record.nz << ","
        << record.mesh_nodes << ","
        << record.mesh_tets << ","
        << record.mesh_edges << ","
        << record.ndof << ","
        << record.nnz_lower_s << ","
        << record.nnz_lower_t << ","
        << record.assembly_ms << ","
        << record.solve_ms << ","
        << record.post_ms << ","
        << record.total_ms << "\n";

    return static_cast<bool>(out);
}

} // namespace run_timing_fem3d_csv
