/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec3/helmvec3_field_output.hpp                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exportacao espacial didatica para pontos casados do HELMVEC3.   */
/*****************************************************************************/
/* Observacao: A Eq. (136) usa x = [Et ; Ez]. Por isso, cada ponto exportado  */
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
#include "helmvec2/helmvec2_field_output.hpp"
#include "helmvec2/helmvec2_coupled_system.hpp"
#include "helmvec3/helmvec3_case_output.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace helmvec3_field
{

struct SpatialArtifacts
{
    std::string field_status;
    std::string et_fields_csv_file;
    std::string ez_fields_csv_file;
    std::string et_vtk_file;
    std::string ez_vtk_file;
};

inline SpatialArtifacts export_mode(
    const helmvec3_output::CaseDirs &dirs,
    const Mesh2D &mesh,
    const CoupledBetaSystem &sys,
    const std::vector<double> &vr_col,
    int eig_index,
    const std::string &stem)
{
    const std::vector<double> edge_block =
        helmvec2_field::extract_edge_block_mode_as_column_major(vr_col, sys.nt, sys.nz, eig_index);
    const std::vector<double> scalar_block =
        helmvec2_field::extract_scalar_block_mode_as_column_major(vr_col, sys.nt, sys.nz, eig_index);

    std::vector<double> cell_vx;
    std::vector<double> cell_vy;
    helmvec_post::reconstruct_cell_field_from_edge_mode(
        mesh, sys.edge, edge_block, 0, cell_vx, cell_vy);
    const std::vector<double> nodal_scalar =
        helm10_post::extract_mode_nodal_from_Z(mesh, sys.scal, scalar_block, 0, true);

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
    if (!helmvec2_field::write_scalar_mode_fields_csv(
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

} // namespace helmvec3_field
