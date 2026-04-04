/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_mode_match.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios para extrair blocos modais do sistema misto e      */
/* reutilizar os matchers por correlacao de massa das formulacoes escalar e   */
/* vetorial.                                                                  */
/*****************************************************************************/

#pragma once

#include <cstddef>
#include <vector>

namespace helmvec1_match
{

/******************************************************************************/
/* FUNCAO: extract_block_mode_as_column_major                                 */
/* DESCRICAO: Extrai um bloco de um autovetor do sistema misto e o devolve    */
/* como uma matriz coluna-major de uma unica coluna, compativel com os        */
/* matchers ja existentes do HELM10 e do HELMVEC.                             */
/* ENTRADA: mixed_zcol: const std::vector<double> &; total_ndof: int;         */
/* mode_idx: int; block_offset: int; block_size: int.                         */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> extract_block_mode_as_column_major(
    const std::vector<double> &mixed_zcol,
    int total_ndof,
    int mode_idx,
    int block_offset,
    int block_size)
{
    std::vector<double> out((size_t)block_size, 0.0);
    const size_t col_offset = (size_t)mode_idx * (size_t)total_ndof;
    for (int i = 0; i < block_size; ++i)
        out[(size_t)i] = mixed_zcol[col_offset + (size_t)(block_offset + i)];
    return out;
}

/******************************************************************************/
/* FUNCAO: extract_edge_block_mode_as_column_major                            */
/* DESCRICAO: Extrai a parte de aresta do autovetor misto.                    */
/* ENTRADA: mixed_zcol: const std::vector<double> &; nt: int; nz: int;        */
/* mode_idx: int.                                                             */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> extract_edge_block_mode_as_column_major(
    const std::vector<double> &mixed_zcol,
    int nt,
    int nz,
    int mode_idx)
{
    return extract_block_mode_as_column_major(mixed_zcol, nt + nz, mode_idx, 0, nt);
}

/******************************************************************************/
/* FUNCAO: extract_scalar_block_mode_as_column_major                          */
/* DESCRICAO: Extrai a parte escalar do autovetor misto.                      */
/* ENTRADA: mixed_zcol: const std::vector<double> &; nt: int; nz: int;        */
/* mode_idx: int.                                                             */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> extract_scalar_block_mode_as_column_major(
    const std::vector<double> &mixed_zcol,
    int nt,
    int nz,
    int mode_idx)
{
    return extract_block_mode_as_column_major(mixed_zcol, nt + nz, mode_idx, nt, nz);
}

} // namespace helmvec1_match
