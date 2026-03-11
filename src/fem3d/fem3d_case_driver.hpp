/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_case_driver.hpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Infraestrutura de comparacao dos casos 3D com tabelas do artigo.*/
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh3d_cylinder.hpp"
#include "core/mesh3d_rect.hpp"
#include "core/mesh3d_sphere.hpp"
#include "fem3d_reference_tables.hpp"
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace fem3d
{
enum class CaseId
{
  air,
  half,
  cyl,
  sphere
};

struct CliDefaults
{
  bool run_air = true;
  bool run_half = false;
  bool run_cyl = false;
  bool run_sphere = false;
};

struct CliOptions
{
  bool run_air = true;
  bool run_half = false;
  bool run_cyl = false;
  bool run_sphere = false;

  bool custom_mesh = false;
  int nx = 0;
  int ny = 0;
  int nz = 0;
};

struct PreparedCase
{
  CaseId id = CaseId::air;
  std::string header;
  Mesh3D mesh;
  std::vector<double> eps_r_tet;
  std::vector<double> mu_r_tet;
  std::vector<RefRow> rows;
};

/******************************************************************************/
/* FUNCAO: print_usage                                                        */
/* DESCRICAO: Imprime no terminal as opcoes de linha de comando e a forma     */
/* correta de executar o programa.                                            */
/* ENTRADA: bin_name: const char *.                                           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_usage(const char *bin_name)
{
  std::cout << "Usage: " << bin_name
            << " [--air|--half|--cyl|--sphere|--all] [--nx N --ny N --nz N]\n";
}

/******************************************************************************/
/* FUNCAO: parse_cli                                                          */
/* DESCRICAO: Analisa os argumentos de linha de comando e monta a configuracao*/
/* de execucao de forma consistente.                                          */
/* ENTRADA: argc: int; argv: char **; defaults: const CliDefaults &;          */
/* bin_name: const char *.                                                    */
/* SAIDA: std::optional<CliOptions>.                                          */
/******************************************************************************/
inline std::optional<CliOptions> parse_cli(
    int argc,
    char **argv,
    const CliDefaults &defaults,
    const char *bin_name)
{
  CliOptions opt;
  opt.run_air = defaults.run_air;
  opt.run_half = defaults.run_half;
  opt.run_cyl = defaults.run_cyl;
  opt.run_sphere = defaults.run_sphere;

  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "--air")
    {
      opt.run_air = true;
      opt.run_half = false;
      opt.run_cyl = false;
      opt.run_sphere = false;
    }
    else if (a == "--half")
    {
      opt.run_air = false;
      opt.run_half = true;
      opt.run_cyl = false;
      opt.run_sphere = false;
    }
    else if (a == "--cyl")
    {
      opt.run_air = false;
      opt.run_half = false;
      opt.run_cyl = true;
      opt.run_sphere = false;
    }
    else if (a == "--sphere")
    {
      opt.run_air = false;
      opt.run_half = false;
      opt.run_cyl = false;
      opt.run_sphere = true;
    }
    else if (a == "--all")
    {
      opt.run_air = true;
      opt.run_half = true;
      opt.run_cyl = true;
      opt.run_sphere = true;
    }
    else if (a == "--nx" && i + 1 < argc)
    {
      opt.nx = std::atoi(argv[++i]);
      opt.custom_mesh = true;
    }
    else if (a == "--ny" && i + 1 < argc)
    {
      opt.ny = std::atoi(argv[++i]);
      opt.custom_mesh = true;
    }
    else if (a == "--nz" && i + 1 < argc)
    {
      opt.nz = std::atoi(argv[++i]);
      opt.custom_mesh = true;
    }
    else if (a == "--help")
    {
      print_usage(bin_name);
      return std::nullopt;
    }
  }
  return opt;
}

/******************************************************************************/
/* FUNCAO: selected_cases                                                     */
/* DESCRICAO: Converte flags da CLI em lista ordenada de casos efetivamente   */
/* selecionados.                                                              */
/* ENTRADA: opt: const CliOptions &.                                          */
/* SAIDA: std::vector<CaseId>.                                                */
/******************************************************************************/
inline std::vector<CaseId> selected_cases(const CliOptions &opt)
{
  std::vector<CaseId> out;
  if (opt.run_air)
    out.push_back(CaseId::air);
  if (opt.run_half)
    out.push_back(CaseId::half);
  if (opt.run_cyl)
    out.push_back(CaseId::cyl);
  if (opt.run_sphere)
    out.push_back(CaseId::sphere);
  return out;
}

/******************************************************************************/
/* FUNCAO: default_grid_for_case                                              */
/* DESCRICAO: Retorna a discretizacao padrao recomendada para cada caso 3D de */
/* validacao (Figura 15, Figura 16, Figura 17 e Tabela 15).                   */
/* ENTRADA: id: CaseId.                                                       */
/* SAIDA: Grid3D.                                                             */
/******************************************************************************/
inline Grid3D default_grid_for_case(CaseId id)
{
  switch (id)
  {
  case CaseId::air:
    return default_grid_fig15();
  case CaseId::half:
    return default_grid_fig16();
  case CaseId::cyl:
    return default_grid_fig17();
  case CaseId::sphere:
    return default_grid_table15();
  }
  return default_grid_fig15();
}

/******************************************************************************/
/* FUNCAO: build_case                                                         */
/* DESCRICAO: Prepara geometria, materiais e tabela de referencia para um caso*/
/* 3D especifico, seguindo exatamente os cenarios de validacao da Secao 3.1.5.*/
/* ENTRADA: id: CaseId; g: const Grid3D &; solver_tag: const char *.          */
/* SAIDA: PreparedCase.                                                       */
/******************************************************************************/
inline PreparedCase build_case(CaseId id, const Grid3D &g, const char *solver_tag)
{
  PreparedCase out;
  out.id = id;

  switch (id)
  {
  case CaseId::air:
  {
    const auto geom = fig15_geom();
    out.mesh = make_rect_tet_mesh(geom.lx, geom.ly, geom.lz, g.nx, g.ny, g.nz);
    out.eps_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.mu_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.rows = table12_rows();
    out.header = std::string("[") + solver_tag + "] Figure 15 / Table 12 - Air-filled rectangular cavity";
    break;
  }
  case CaseId::half:
  {
    const auto geom = fig16_geom();
    out.mesh = make_rect_tet_mesh(geom.lx, geom.ly, geom.lz, g.nx, g.ny, g.nz);
    out.eps_r_tet = make_eps_r_tets_by_z(out.mesh, 0.5 * geom.lz, 1.0, 2.0);
    out.mu_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.rows = table13_rows();
    out.header = std::string("[") + solver_tag + "] Figure 16 / Table 13 - Half-filled rectangular cavity";
    break;
  }
  case CaseId::cyl:
  {
    const auto geom = fig17_geom();
    out.mesh = make_cylinder_tet_mesh_cartesian(geom.radius, geom.height, g.nx, g.ny, g.nz);
    out.eps_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.mu_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.rows = table14_rows();
    out.header = std::string("[") + solver_tag + "] Figure 17 / Table 14 - Air-filled circular cylindrical cavity";
    break;
  }
  case CaseId::sphere:
  {
    const auto geom = table15_geom();
    out.mesh = make_sphere_tet_mesh_cartesian(geom.radius, g.nx, g.ny, g.nz);
    out.eps_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.mu_r_tet.assign(out.mesh.tets.size(), 1.0);
    out.rows = table15_rows();
    out.header = std::string("[") + solver_tag + "] Table 15 - Air-filled spherical cavity";
    break;
  }
  }
  return out;
}

template <typename CaseFn>
/******************************************************************************/
/* FUNCAO: for_each_selected_case                                             */
/* DESCRICAO: Itera sobre os casos ativos e executa o callback de             */
/* processamento para cada um deles.                                          */
/* ENTRADA: opt: const CliOptions &; solver_tag: const char *; run_case:      */
/* CaseFn.                                                                    */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void for_each_selected_case(
    const CliOptions &opt,
    const char *solver_tag,
    CaseFn run_case)
{
  for (CaseId id : selected_cases(opt))
  {
    const Grid3D g = opt.custom_mesh ? Grid3D{opt.nx, opt.ny, opt.nz} : default_grid_for_case(id);
    run_case(build_case(id, g, solver_tag));
  }
}

} // namespace fem3d
