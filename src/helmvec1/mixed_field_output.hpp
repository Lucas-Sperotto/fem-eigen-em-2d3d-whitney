/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_field_output.hpp                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exportacao espacial didatica para modos da familia HELMVEC1.    */
/*****************************************************************************/
/* Observacao: No sistema misto da Eq. (92), modos dominados pelo bloco edge  */
/* sao visualizados como campos vetoriais transversais (`Et` ou `Ht`),        */
/* enquanto modos dominados pelo bloco scalar sao visualizados como campos     */
/* longitudinais escalares (`Ez` ou `Hz`).                                    */
/*****************************************************************************/

#pragma once

#include "core/io_vtk.hpp"
#include "core/io_vtk_sv.hpp"
#include "core/mesh2d.hpp"
#include "helm10/scalar_mode_post.hpp"
#include "helmvec/edge_mode_post.hpp"
#include "helmvec1/helmvec1_mixed_system.hpp"
#include "helmvec1/mixed_case_output.hpp"
#include "helmvec1/mixed_mode_match.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace helmvec1_field
{

struct SpatialArtifacts
{
    std::string field_data_kind;
    std::string field_status;
    std::string fields_csv_file;
    std::string vtk_file;
};

inline std::string make_mode_stem(
    const std::string &case_name,
    const helmvec1_output::ModeCsvRecord &rec)
{
    std::ostringstream oss;
    oss << "mixed_" << case_name
        << "_" << rec.formulation
        << "_" << rec.component_label
        << "_" << rec.mode_label
        << "_rank" << std::setw(2) << std::setfill('0') << rec.positive_rank;
    return oss.str();
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

inline bool write_edge_mode_fields_csv(
    const std::filesystem::path &path,
    const Mesh2D &mesh,
    const std::vector<double> &cell_vx,
    const std::vector<double> &cell_vy,
    const std::string &component_label)
{
    if (mesh.tris.size() != cell_vx.size() || mesh.tris.size() != cell_vy.size())
        return false;

    const std::string prefix =
        (component_label == "Ht") ? "H" : "E";
    const std::string x_name = prefix + "x";
    const std::string y_name = prefix + "y";
    const std::string mag_name = prefix + "mag";

    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "cell_id,xc_m,yc_m,"
        << x_name << ","
        << y_name << ","
        << mag_name << "\n";
    for (size_t cell_id = 0; cell_id < mesh.tris.size(); ++cell_id)
    {
        const Tri &tri = mesh.tris[cell_id];
        const Node2D &n0 = mesh.nodes[(size_t)tri.v[0]];
        const Node2D &n1 = mesh.nodes[(size_t)tri.v[1]];
        const Node2D &n2 = mesh.nodes[(size_t)tri.v[2]];
        const double xc = (n0.x + n1.x + n2.x) / 3.0;
        const double yc = (n0.y + n1.y + n2.y) / 3.0;
        const double fx = cell_vx[cell_id];
        const double fy = cell_vy[cell_id];
        const double fmag = std::sqrt(fx * fx + fy * fy);
        out << cell_id << ","
            << xc << ","
            << yc << ","
            << fx << ","
            << fy << ","
            << fmag << "\n";
    }

    return static_cast<bool>(out);
}

inline SpatialArtifacts export_edge_mode(
    const std::string &case_name,
    const helmvec1_output::CaseDirs &dirs,
    const Mesh2D &mesh,
    const MixedSystem92 &sys,
    const std::vector<double> &mixed_zcol,
    int eig_index,
    const helmvec1_output::ModeCsvRecord &rec)
{
    const std::vector<double> block_zcol =
        helmvec1_match::extract_edge_block_mode_as_column_major(
            mixed_zcol, sys.nt, sys.nz, eig_index);

    std::vector<double> cell_vx;
    std::vector<double> cell_vy;
    helmvec_post::reconstruct_cell_field_from_edge_mode(
        mesh, sys.edge, block_zcol, 0, cell_vx, cell_vy);

    const std::string stem = make_mode_stem(case_name, rec);
    const std::string fields_csv_file = stem + "_fields.csv";
    const std::string vtk_file = stem + ".vtk";

    if (!write_edge_mode_fields_csv(
            dirs.csv / fields_csv_file,
            mesh,
            cell_vx,
            cell_vy,
            rec.component_label))
    {
        throw std::runtime_error("Erro ao escrever CSV espacial de modo edge: " +
                                 fields_csv_file);
    }

    write_vtk_unstructured_tri_cell_vector(
        (dirs.vtk / vtk_file).string(),
        mesh,
        cell_vx,
        cell_vy,
        rec.component_label);

    SpatialArtifacts artifacts;
    artifacts.field_data_kind = "edge_vector_cell";
    artifacts.field_status = "cell_centroid_unit_peak_normalized";
    artifacts.fields_csv_file = fields_csv_file;
    artifacts.vtk_file = vtk_file;
    return artifacts;
}

inline SpatialArtifacts export_scalar_mode(
    const std::string &case_name,
    const helmvec1_output::CaseDirs &dirs,
    const Mesh2D &mesh,
    const MixedSystem92 &sys,
    const std::vector<double> &mixed_zcol,
    int eig_index,
    const helmvec1_output::ModeCsvRecord &rec)
{
    const std::vector<double> block_zcol =
        helmvec1_match::extract_scalar_block_mode_as_column_major(
            mixed_zcol, sys.nt, sys.nz, eig_index);
    const std::vector<double> nodal_scalar =
        helm10_post::extract_mode_nodal_from_Z(
            mesh, sys.scal, block_zcol, 0, true);

    const std::string stem = make_mode_stem(case_name, rec);
    const std::string fields_csv_file = stem + "_fields.csv";
    const std::string vtk_file = stem + ".vtk";

    if (!write_scalar_mode_fields_csv(
            dirs.csv / fields_csv_file,
            mesh,
            nodal_scalar,
            rec.component_label))
    {
        throw std::runtime_error("Erro ao escrever CSV espacial de modo scalar: " +
                                 fields_csv_file);
    }

    write_vtk_unstructured_tri_scalar(
        (dirs.vtk / vtk_file).string(),
        mesh,
        nodal_scalar,
        rec.component_label);

    SpatialArtifacts artifacts;
    artifacts.field_data_kind = "scalar_nodal";
    artifacts.field_status = "nodal_unit_peak_normalized";
    artifacts.fields_csv_file = fields_csv_file;
    artifacts.vtk_file = vtk_file;
    return artifacts;
}

} // namespace helmvec1_field
