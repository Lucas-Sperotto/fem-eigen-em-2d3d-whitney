/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/output_paths.hpp                                         */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Padronizacao de caminhos de saida para resultados e validacoes. */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Infraestrutura de   */
/* execucao.                                                                  */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once

#include <cstdlib>
#include <filesystem>
#include <string>

namespace output_paths
{
/******************************************************************************/
/* FUNCAO: repo_root                                                          */
/* DESCRICAO: Resolve o diretorio raiz do projeto para construir caminhos de saida independentes do cwd. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::filesystem::path.                                              */
/******************************************************************************/
inline std::filesystem::path repo_root()
{
#ifdef TP3485_SOURCE_DIR
  return std::filesystem::path(TP3485_SOURCE_DIR);
#else
  return std::filesystem::current_path();
#endif
}

/******************************************************************************/
/* FUNCAO: out_root                                                           */
/* DESCRICAO: Retorna e cria (se necessario) a pasta principal out/ para centralizar resultados da execucao. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::filesystem::path.                                              */
/******************************************************************************/
inline std::filesystem::path out_root()
{
  if (const char *env = std::getenv("TP3485_OUT_DIR"))
  {
    if (env[0] != '\0')
      return std::filesystem::path(env);
  }
  return repo_root() / "out";
}

/******************************************************************************/
/* FUNCAO: ensure_case_dir                                                    */
/* DESCRICAO: Garante a existencia da pasta de um caso numerico e devolve seu caminho absoluto. */
/* ENTRADA: relative_case_dir: const std::string &.                           */
/* SAIDA: std::filesystem::path.                                              */
/******************************************************************************/
inline std::filesystem::path ensure_case_dir(const std::string &relative_case_dir)
{
  const auto p = out_root() / relative_case_dir;
  std::error_code ec;
  std::filesystem::create_directories(p, ec);
  return p;
}

/******************************************************************************/
/* FUNCAO: file_in                                                            */
/* DESCRICAO: Monta o caminho completo de um arquivo de saida dentro da pasta de um caso. */
/* ENTRADA: dir: const std::filesystem::path &; name: const std::string &.    */
/* SAIDA: std::string.                                                        */
/******************************************************************************/
inline std::string file_in(const std::filesystem::path &dir, const std::string &name)
{
  return (dir / name).string();
}
} // namespace output_paths

