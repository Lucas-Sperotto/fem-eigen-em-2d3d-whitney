/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_field_output.hpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exportacao espacial didatica para modos 3D casados.             */
/*****************************************************************************/

#pragma once

#include "core/io_vtk_3d.hpp"
#include "core/mesh3d.hpp"
#include "edge3d/edge3d_basis.hpp"
#include "edge3d/edge3d_dofs.hpp"
#include "fem3d/fem3d_case_output.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem3d_field
{

struct SpatialArtifacts
{
    std::string field_status;
    std::string fields_csv_file;
    std::string vtk_file;
};

inline std::vector<double> extract_mode_as_column_major(
    const std::vector<double> &zcol,
    int ndof,
    int eig_index)
{
    std::vector<double> out((size_t)ndof, 0.0);
    const size_t col_offset = (size_t)eig_index * (size_t)ndof;
    for (int i = 0; i < ndof; ++i)
        out[(size_t)i] = zcol[col_offset + (size_t)i];
    return out;
}

inline std::string sanitize_mode_label(const std::string &label)
{
    std::string out;
    out.reserve(label.size());
    for (char ch : label)
    {
        const bool is_alnum =
            (ch >= '0' && ch <= '9') ||
            (ch >= 'A' && ch <= 'Z') ||
            (ch >= 'a' && ch <= 'z');
        out.push_back(is_alnum ? ch : '_');
    }
    return out;
}

inline bool write_mode_fields_csv(
    const std::filesystem::path &path,
    const Mesh3D &mesh,
    const std::vector<double> &xc,
    const std::vector<double> &yc,
    const std::vector<double> &zc,
    const std::vector<double> &ex,
    const std::vector<double> &ey,
    const std::vector<double> &ez,
    const std::vector<double> &emag)
{
    if (mesh.tets.size() != ex.size() ||
        mesh.tets.size() != ey.size() ||
        mesh.tets.size() != ez.size() ||
        mesh.tets.size() != emag.size())
        return false;

    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "tet_id,xc_m,yc_m,zc_m,Ex,Ey,Ez,Emag\n";
    for (size_t tet_id = 0; tet_id < mesh.tets.size(); ++tet_id)
    {
        out << tet_id << ","
            << xc[tet_id] << ","
            << yc[tet_id] << ","
            << zc[tet_id] << ","
            << ex[tet_id] << ","
            << ey[tet_id] << ","
            << ez[tet_id] << ","
            << emag[tet_id] << "\n";
    }
    return static_cast<bool>(out);
}

inline SpatialArtifacts export_mode(
    const fem3d_output::CaseDirs &dirs,
    const Mesh3D &mesh,
    const EdgeDofs3D &ed,
    const std::vector<double> &zcol,
    const std::string &artifact_prefix,
    int eig_index,
    int mode_rank,
    const std::string &mode_label)
{
    const std::vector<double> mode = extract_mode_as_column_major(zcol, ed.ndof, eig_index);
    const std::array<double, 4> centroid_lambda{{0.25, 0.25, 0.25, 0.25}};

    std::vector<double> xc(mesh.tets.size(), 0.0);
    std::vector<double> yc(mesh.tets.size(), 0.0);
    std::vector<double> zc(mesh.tets.size(), 0.0);
    std::vector<double> ex(mesh.tets.size(), 0.0);
    std::vector<double> ey(mesh.tets.size(), 0.0);
    std::vector<double> ez(mesh.tets.size(), 0.0);
    std::vector<double> emag(mesh.tets.size(), 0.0);

    double max_mag = 0.0;
    for (size_t tet_id = 0; tet_id < mesh.tets.size(); ++tet_id)
    {
        const Tet &t = mesh.tets[tet_id];
        const TetGeomEdge tg = tet_geom_edge(mesh, t);
        const TetEdges &te = ed.tet_edges[tet_id];

        Vec3d field{};
        for (int m = 0; m < 6; ++m)
        {
            const int global_edge = te.e[(size_t)m];
            const int dof = ed.edge_to_dof[(size_t)global_edge];
            if (dof < 0)
                continue;
            const double coeff = (double)te.sgn[(size_t)m] * mode[(size_t)dof];
            const Vec3d Wm = whitney_W_local_3d(m, tg, centroid_lambda);
            field.x += coeff * Wm.x;
            field.y += coeff * Wm.y;
            field.z += coeff * Wm.z;
        }

        const Node3D &n0 = mesh.nodes[(size_t)t.v[0]];
        const Node3D &n1 = mesh.nodes[(size_t)t.v[1]];
        const Node3D &n2 = mesh.nodes[(size_t)t.v[2]];
        const Node3D &n3 = mesh.nodes[(size_t)t.v[3]];
        xc[tet_id] = 0.25 * (n0.x + n1.x + n2.x + n3.x);
        yc[tet_id] = 0.25 * (n0.y + n1.y + n2.y + n3.y);
        zc[tet_id] = 0.25 * (n0.z + n1.z + n2.z + n3.z);
        ex[tet_id] = field.x;
        ey[tet_id] = field.y;
        ez[tet_id] = field.z;
        emag[tet_id] = std::sqrt(field.x * field.x + field.y * field.y + field.z * field.z);
        max_mag = std::max(max_mag, emag[tet_id]);
    }

    if (max_mag > 0.0)
    {
        for (size_t tet_id = 0; tet_id < mesh.tets.size(); ++tet_id)
        {
            ex[tet_id] /= max_mag;
            ey[tet_id] /= max_mag;
            ez[tet_id] /= max_mag;
            emag[tet_id] /= max_mag;
        }
    }

    std::ostringstream stem;
    stem << artifact_prefix
         << "_mode" << std::setw(2) << std::setfill('0') << mode_rank
         << "_" << sanitize_mode_label(mode_label);
    const std::string fields_csv_file = stem.str() + "_E_fields.csv";
    const std::string vtk_file = stem.str() + "_E.vtk";

    if (!write_mode_fields_csv(
            dirs.csv / fields_csv_file,
            mesh,
            xc,
            yc,
            zc,
            ex,
            ey,
            ez,
            emag))
    {
        throw std::runtime_error("Erro ao escrever CSV de campo 3D: " + fields_csv_file);
    }

    write_vtk_unstructured_tet_cell_vector(
        (dirs.vtk / vtk_file).string(),
        mesh,
        ex,
        ey,
        ez,
        emag,
        "E",
        "Emag");

    SpatialArtifacts artifacts;
    artifacts.field_status = "E_cell_centroid_unit_peak_normalized";
    artifacts.fields_csv_file = fields_csv_file;
    artifacts.vtk_file = vtk_file;
    return artifacts;
}

} // namespace fem3d_field
