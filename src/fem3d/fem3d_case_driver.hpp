/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_case_driver.hpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Infraestrutura de comparacao dos casos 3D com tabelas do artigo.*/
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/assembly_backend.hpp"
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
  ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm;
  bool debug_local_blocks = false;
  bool debug_candidates = false;

  bool custom_mesh = false;
  int nx = 0;
  int ny = 0;
  int nz = 0;
};

struct PreparedCase
{
  CaseId id = CaseId::air;
  std::string case_name;
  std::string header;
  Grid3D grid;
  Mesh3D mesh;
  std::vector<double> eps_r_tet;
  std::vector<double> mu_r_tet;
  std::vector<RefRow> rows;
};

inline void print_usage(const char *bin_name);
inline void print_single_case_usage(const char *bin_name);

inline int parse_positive_cli_int_strict(const char *text, const char *name)
{
  try
  {
    std::string s(text);
    size_t idx = 0;
    const int value = std::stoi(s, &idx);
    if (idx != s.size())
      throw std::runtime_error(std::string(name) + " invalido: " + text);
    if (value <= 0)
      throw std::runtime_error(std::string(name) + " deve ser > 0");
    return value;
  }
  catch (const std::runtime_error &)
  {
    throw;
  }
  catch (const std::exception &)
  {
    throw std::runtime_error(std::string(name) + " invalido: " + text);
  }
}

[[noreturn]] inline void print_cli_error_and_exit(
    const std::string &message,
    const char *bin_name,
    bool single_case)
{
  std::cerr << "Erro: " << message << "\n";
  if (single_case)
    print_single_case_usage(bin_name);
  else
    print_usage(bin_name);
  std::exit(2);
}

/******************************************************************************/
/* FUNCAO: case_defaults                                                     */
/* DESCRICAO: Retorna a configuracao padrao de selecao para um unico caso 3D.*/
/* ENTRADA: id: CaseId.                                                      */
/* SAIDA: CliDefaults.                                                       */
/******************************************************************************/
inline CliDefaults case_defaults(CaseId id)
{
  CliDefaults out{/*run_air=*/false, /*run_half=*/false, /*run_cyl=*/false, /*run_sphere=*/false};
  switch (id)
  {
  case CaseId::air:
    out.run_air = true;
    break;
  case CaseId::half:
    out.run_half = true;
    break;
  case CaseId::cyl:
    out.run_cyl = true;
    break;
  case CaseId::sphere:
    out.run_sphere = true;
    break;
  }
  return out;
}

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
            << " [--air|--half|--cyl|--sphere|--all] [--nx N --ny N --nz N]"
            << " [--backend closed-form|gauss]"
            << " [--debug-local-blocks] [--debug-candidates]\n";
}

/******************************************************************************/
/* FUNCAO: print_single_case_usage                                           */
/* DESCRICAO: Imprime o uso de um executavel 3D fixo em um unico caso.       */
/* ENTRADA: bin_name: const char *.                                          */
/* SAIDA: sem retorno explicito (void).                                      */
/******************************************************************************/
inline void print_single_case_usage(const char *bin_name)
{
  std::cout << "Usage: " << bin_name
            << " [--nx N --ny N --nz N]"
            << " [--backend closed-form|gauss]"
            << " [--debug-local-blocks] [--debug-candidates]\n";
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
    else if (a == "--nx")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --nx", bin_name, false);
      try
      {
        opt.nx = parse_positive_cli_int_strict(argv[++i], "nx");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, false);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--ny")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --ny", bin_name, false);
      try
      {
        opt.ny = parse_positive_cli_int_strict(argv[++i], "ny");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, false);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--nz")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --nz", bin_name, false);
      try
      {
        opt.nz = parse_positive_cli_int_strict(argv[++i], "nz");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, false);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--backend")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --backend", bin_name, false);
      try
      {
        opt.backend = parse_element_assembly_backend(argv[++i]);
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, false);
      }
    }
    else if (a == "--debug" || a == "--debug-all")
    {
      opt.debug_local_blocks = true;
      opt.debug_candidates = true;
    }
    else if (a == "--debug-local-blocks")
    {
      opt.debug_local_blocks = true;
    }
    else if (a == "--debug-candidates")
    {
      opt.debug_candidates = true;
    }
    else if (a.rfind("--backend=", 0) == 0)
    {
      try
      {
        opt.backend = parse_element_assembly_backend(a.substr(std::string("--backend=").size()));
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, false);
      }
    }
    else if (a == "--help")
    {
      print_usage(bin_name);
      return std::nullopt;
    }
    else if (!a.empty() && a[0] == '-')
    {
      print_cli_error_and_exit("opcao desconhecida: " + a, bin_name, false);
    }
    else
    {
      print_cli_error_and_exit("argumento posicional inesperado: " + a, bin_name, false);
    }
  }
  return opt;
}

/******************************************************************************/
/* FUNCAO: parse_single_case_cli                                             */
/* DESCRICAO: Analisa a CLI de um executavel 3D dedicado a um unico caso.    */
/* ENTRADA: argc: int; argv: char **; fixed_case: CaseId; bin_name: const char*. */
/* SAIDA: std::optional<CliOptions>.                                          */
/******************************************************************************/
inline std::optional<CliOptions> parse_single_case_cli(
    int argc,
    char **argv,
    CaseId fixed_case,
    const char *bin_name)
{
  CliOptions opt;
  const auto defaults = case_defaults(fixed_case);
  opt.run_air = defaults.run_air;
  opt.run_half = defaults.run_half;
  opt.run_cyl = defaults.run_cyl;
  opt.run_sphere = defaults.run_sphere;

  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "--air" || a == "--half" || a == "--cyl" || a == "--sphere" || a == "--all")
    {
      print_cli_error_and_exit(
          std::string(bin_name) +
              " ja representa um caso 3D especifico; nao use flags de selecao de caso.",
          bin_name,
          true);
    }
    else if (a == "--nx")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --nx", bin_name, true);
      try
      {
        opt.nx = parse_positive_cli_int_strict(argv[++i], "nx");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, true);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--ny")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --ny", bin_name, true);
      try
      {
        opt.ny = parse_positive_cli_int_strict(argv[++i], "ny");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, true);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--nz")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --nz", bin_name, true);
      try
      {
        opt.nz = parse_positive_cli_int_strict(argv[++i], "nz");
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, true);
      }
      opt.custom_mesh = true;
    }
    else if (a == "--backend")
    {
      if (i + 1 >= argc)
        print_cli_error_and_exit("faltou valor apos --backend", bin_name, true);
      try
      {
        opt.backend = parse_element_assembly_backend(argv[++i]);
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, true);
      }
    }
    else if (a == "--debug" || a == "--debug-all")
    {
      opt.debug_local_blocks = true;
      opt.debug_candidates = true;
    }
    else if (a == "--debug-local-blocks")
    {
      opt.debug_local_blocks = true;
    }
    else if (a == "--debug-candidates")
    {
      opt.debug_candidates = true;
    }
    else if (a.rfind("--backend=", 0) == 0)
    {
      try
      {
        opt.backend = parse_element_assembly_backend(a.substr(std::string("--backend=").size()));
      }
      catch (const std::exception &e)
      {
        print_cli_error_and_exit(e.what(), bin_name, true);
      }
    }
    else if (a == "--help")
    {
      print_single_case_usage(bin_name);
      return std::nullopt;
    }
    else if (!a.empty() && a[0] == '-')
    {
      print_cli_error_and_exit("opcao desconhecida: " + a, bin_name, true);
    }
    else
    {
      print_cli_error_and_exit("argumento posicional inesperado: " + a, bin_name, true);
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
  out.grid = g;

  switch (id)
  {
  case CaseId::air:
  {
    out.case_name = "air";
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
    out.case_name = "half";
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
    out.case_name = "cyl";
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
    out.case_name = "sphere";
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
