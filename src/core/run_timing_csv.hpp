/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/run_timing_csv.hpp                                       */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Escrita estruturada dos tempos de execucao e da configuracao    */
/* do caso em CSV, para auditoria e comparacao entre rodadas.                 */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iomanip>
#include <string>

namespace run_timing_csv
{

struct Record
{
    std::string case_label;
    std::string geometry;
    std::string backend;
    int legacy_cli_used = 0;
    int debug_local_blocks = 0;
    int debug_candidates = 0;
    int frequency_was_provided = 0;
    double frequency_input_hz = 0.0;
    std::string frequency_policy;
    double frequency_used_hz = 0.0;
    double eps_r = 1.0;
    double mu_r = 1.0;
    double k_medium = 0.0;
    double reference_tm_kc = 0.0;
    double reference_tm_ztm_actual = 0.0;
    std::string ar_m;
    std::string b_m;
    std::string r_m;
    std::string r1_m;
    std::string r2_m;
    std::string nx;
    std::string ny;
    std::string nr;
    std::string nt;
    int nmodos = 0;
    int mesh_nodes = 0;
    int mesh_tris = 0;
    double te_assembly_ms = 0.0;
    double te_solve_ms = 0.0;
    double tm_assembly_ms = 0.0;
    double tm_solve_ms = 0.0;
    double assembly_ms_total = 0.0;
    double solve_ms_total = 0.0;
    double post_ms = 0.0;
    double total_ms = 0.0;
};

inline bool write_csv(const std::string &path, const Record &record)
{
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "case_label,geometry,backend,legacy_cli_used,debug_local_blocks,debug_candidates,"
           "frequency_was_provided,frequency_input_hz,frequency_policy,frequency_used_hz,"
           "eps_r,mu_r,k_medium,reference_tm_kc,reference_tm_ztm_actual,"
           "ar_m,b_m,r_m,r1_m,r2_m,nx,ny,nr,nt,nmodos,mesh_nodes,mesh_tris,"
           "te_assembly_ms,te_solve_ms,tm_assembly_ms,tm_solve_ms,"
           "assembly_ms_total,solve_ms_total,post_ms,total_ms\n";

    out << record.case_label << ","
        << record.geometry << ","
        << record.backend << ","
        << record.legacy_cli_used << ","
        << record.debug_local_blocks << ","
        << record.debug_candidates << ","
        << record.frequency_was_provided << ","
        << record.frequency_input_hz << ","
        << record.frequency_policy << ","
        << record.frequency_used_hz << ","
        << record.eps_r << ","
        << record.mu_r << ","
        << record.k_medium << ","
        << record.reference_tm_kc << ","
        << record.reference_tm_ztm_actual << ","
        << record.ar_m << ","
        << record.b_m << ","
        << record.r_m << ","
        << record.r1_m << ","
        << record.r2_m << ","
        << record.nx << ","
        << record.ny << ","
        << record.nr << ","
        << record.nt << ","
        << record.nmodos << ","
        << record.mesh_nodes << ","
        << record.mesh_tris << ","
        << record.te_assembly_ms << ","
        << record.te_solve_ms << ","
        << record.tm_assembly_ms << ","
        << record.tm_solve_ms << ","
        << record.assembly_ms_total << ","
        << record.solve_ms_total << ","
        << record.post_ms << ","
        << record.total_ms << "\n";

    return static_cast<bool>(out);
}

} // namespace run_timing_csv
