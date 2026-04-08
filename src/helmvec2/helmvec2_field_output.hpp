/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/helmvec2_field_output.hpp                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exportacao espacial didatica para modos casados do HELMVEC2.    */
/*****************************************************************************/
/* Observacao: A Eq. (119) usa x = [Et ; Ez]. Por isso, cada modo exportado   */
/* gera dois artefatos espaciais complementares:                              */
/* - o campo vetorial transversal Et por celula;                              */
/* - a componente longitudinal escalar Ez por no.                             */
/*****************************************************************************/

#pragma once

#include "core/io_vtk.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/mesh2d.hpp"
#include "helm10/scalar_mode_post.hpp"
#include "helmvec/edge_case_output.hpp"
#include "helmvec/edge_mode_post.hpp"
#include "helmvec2/helmvec2_case_output.hpp"
#include "helmvec2/helmvec2_coupled_system.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace helmvec2_field
{

struct SpatialArtifacts
{
    std::string field_status;
    std::string et_fields_csv_file;
    std::string ez_fields_csv_file;
    std::string et_vtk_file;
    std::string ez_vtk_file;
};

inline std::vector<double> extract_general_block_mode_as_column_major(
    const std::vector<double> &vr_col,
    int total_ndof,
    int eig_index,
    int block_offset,
    int block_size)
{
    std::vector<double> out((size_t)block_size, 0.0);
    const size_t col_offset = (size_t)eig_index * (size_t)total_ndof;
    for (int i = 0; i < block_size; ++i)
        out[(size_t)i] = vr_col[col_offset + (size_t)(block_offset + i)];
    return out;
}

inline std::vector<double> extract_edge_block_mode_as_column_major(
    const std::vector<double> &vr_col,
    int nt,
    int nz,
    int eig_index)
{
    return extract_general_block_mode_as_column_major(vr_col, nt + nz, eig_index, 0, nt);
}

inline std::vector<double> extract_scalar_block_mode_as_column_major(
    const std::vector<double> &vr_col,
    int nt,
    int nz,
    int eig_index)
{
    return extract_general_block_mode_as_column_major(vr_col, nt + nz, eig_index, nt, nz);
}

inline bool write_scalar_mode_fields_csv(
    const std::filesystem::path &path,
    const Mesh2D &mesh,
    const std::vector<double> &nodal_scalar,
    const std::string &scalar_label)
{
    if (mesh.nodes.size() != nodal_scalar.size())
        return false;

    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "node_id,x_m,y_m," << scalar_label << "\n";
    for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
    {
        const Node2D &node = mesh.nodes[node_id];
        out << node_id << ","
            << node.x << ","
            << node.y << ","
            << nodal_scalar[node_id] << "\n";
    }

    return static_cast<bool>(out);
}

inline std::string make_mode_stem(int mode_index, int candidate_rank)
{
    std::ostringstream oss;
    oss << "helmvec2_rect_mode" << std::setw(2) << std::setfill('0') << mode_index
        << "_cand" << std::setw(2) << std::setfill('0') << candidate_rank;
    return oss.str();
}

inline SpatialArtifacts export_mode(
    const helmvec2_output::CaseDirs &dirs,
    const Mesh2D &mesh,
    const CoupledWaveNumberSystem &sys,
    const std::vector<double> &vr_col,
    int eig_index,
    int mode_index,
    int candidate_rank)
{
    const std::vector<double> edge_block =
        extract_edge_block_mode_as_column_major(vr_col, sys.nt, sys.nz, eig_index);
    const std::vector<double> scalar_block =
        extract_scalar_block_mode_as_column_major(vr_col, sys.nt, sys.nz, eig_index);

    std::vector<double> cell_vx;
    std::vector<double> cell_vy;
    helmvec_post::reconstruct_cell_field_from_edge_mode(
        mesh, sys.edge, edge_block, 0, cell_vx, cell_vy);
    const std::vector<double> nodal_scalar =
        helm10_post::extract_mode_nodal_from_Z(mesh, sys.scal, scalar_block, 0, true);

    const std::string stem = make_mode_stem(mode_index, candidate_rank);
    const std::string et_fields_csv_file = stem + "_Et_fields.csv";
    const std::string ez_fields_csv_file = stem + "_Ez_fields.csv";
    const std::string et_vtk_file = stem + "_Et.vtk";
    const std::string ez_vtk_file = stem + "_Ez.vtk";

    if (!helmvec_output::write_mode_fields_csv(
            dirs.csv / et_fields_csv_file,
            mesh,
            cell_vx,
            cell_vy,
            "E"))
    {
        throw std::runtime_error("Erro ao escrever CSV de Et: " + et_fields_csv_file);
    }
    if (!write_scalar_mode_fields_csv(
            dirs.csv / ez_fields_csv_file,
            mesh,
            nodal_scalar,
            "Ez"))
    {
        throw std::runtime_error("Erro ao escrever CSV de Ez: " + ez_fields_csv_file);
    }

    write_vtk_unstructured_tri_cell_vector(
        (dirs.vtk / et_vtk_file).string(),
        mesh,
        cell_vx,
        cell_vy,
        "Et");
    write_vtk_unstructured_tri_scalar(
        (dirs.vtk / ez_vtk_file).string(),
        mesh,
        nodal_scalar,
        "Ez");

    SpatialArtifacts artifacts;
    artifacts.field_status = "Et_cell_unit_peak_normalized__Ez_nodal_unit_peak_normalized";
    artifacts.et_fields_csv_file = et_fields_csv_file;
    artifacts.ez_fields_csv_file = ez_fields_csv_file;
    artifacts.et_vtk_file = et_vtk_file;
    artifacts.ez_vtk_file = ez_vtk_file;
    return artifacts;
}

} // namespace helmvec2_field
