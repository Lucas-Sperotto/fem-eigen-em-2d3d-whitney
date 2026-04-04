/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/edge_case_output.hpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de saida para a familia HELMVEC.          */
/*****************************************************************************/
/* Observacao: Centraliza a organizacao de diretorios e a escrita do CSV de   */
/* campo por modo, mantendo a mesma logica para rect, circle e coax.          */
/*****************************************************************************/

#pragma once

#include "core/mesh2d.hpp"
#include "core/output_paths.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace helmvec_output
{

struct CaseDirs
{
    std::filesystem::path root;
    std::filesystem::path csv;
    std::filesystem::path vtk;
    std::filesystem::path img;
};

/******************************************************************************/
/* FUNCAO: ensure_case_dirs                                                   */
/* DESCRICAO: Cria a arvore padrao de saida do HELMVEC em out/helmvec/<caso>. */
/* ENTRADA: case_name: const std::string &.                                   */
/* SAIDA: CaseDirs.                                                           */
/******************************************************************************/
inline CaseDirs ensure_case_dirs(const std::string &case_name)
{
    CaseDirs dirs;
    dirs.root = output_paths::ensure_case_dir("helmvec/" + case_name);
    dirs.csv = dirs.root / "csv";
    dirs.vtk = dirs.root / "vtk";
    dirs.img = dirs.root / "img";
    std::error_code ec;
    std::filesystem::create_directories(dirs.csv, ec);
    std::filesystem::create_directories(dirs.vtk, ec);
    std::filesystem::create_directories(dirs.img, ec);
    return dirs;
}

/******************************************************************************/
/* FUNCAO: write_mode_fields_csv                                              */
/* DESCRICAO: Escreve um CSV por modo com o campo reconstruido nos centroides */
/* dos triangulos.                                                            */
/* ENTRADA: path: const std::filesystem::path &; mesh: const Mesh2D &;        */
/* cell_vx: const std::vector<double> &; cell_vy: const std::vector<double> &;*/
/* component_symbol: const std::string &.                                     */
/* SAIDA: bool.                                                               */
/******************************************************************************/
inline bool write_mode_fields_csv(
    const std::filesystem::path &path,
    const Mesh2D &mesh,
    const std::vector<double> &cell_vx,
    const std::vector<double> &cell_vy,
    const std::string &component_symbol)
{
    if (mesh.tris.size() != cell_vx.size() || mesh.tris.size() != cell_vy.size())
        return false;

    std::ofstream out(path);
    if (!out)
        return false;

    const std::string prefix = component_symbol.empty() ? "F" : component_symbol;
    const std::string x_name = prefix + "x";
    const std::string y_name = prefix + "y";
    const std::string mag_name = prefix + "mag";

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

} // namespace helmvec_output
