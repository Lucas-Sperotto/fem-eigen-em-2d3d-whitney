/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_reference_tables.hpp                              */
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
#include "fem3d_compare.hpp"
#include <algorithm>
#include <vector>

namespace fem3d
{
// Geometrias e malhas padrao dos casos 3D da Secao 3.1.
struct Grid3D
{
  int nx = 0;
  int ny = 0;
  int nz = 0;
};

struct RectGeom
{
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
};

struct CylGeom
{
  double radius = 0.0;
  double height = 0.0;
};

struct SphereGeom
{
  double radius = 0.0;
};

/******************************************************************************/
/* FUNCAO: fig15_geom                                                         */
/* DESCRICAO: Retorna geometria do caso da Figura 15 (cavidade retangular sem dielestrico). */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: RectGeom.                                                           */
/******************************************************************************/
inline RectGeom fig15_geom() { return {1.0, 0.5, 0.75}; }
/******************************************************************************/
/* FUNCAO: fig16_geom                                                         */
/* DESCRICAO: Retorna geometria do caso da Figura 16 (cavidade retangular     */
/* meio preenchida por dielestrico).                                          */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: RectGeom.                                                           */
/******************************************************************************/
inline RectGeom fig16_geom() { return {1.0, 1.0, 1.0}; }
/******************************************************************************/
/* FUNCAO: fig17_geom                                                         */
/* DESCRICAO: Retorna geometria do caso da Figura 17 (cavidade cilindrica     */
/* circular sem preenchimento).                                               */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: CylGeom.                                                            */
/******************************************************************************/
inline CylGeom fig17_geom() { return {0.5, 0.5}; } // diameter 1 cm, height 0.5 cm
/******************************************************************************/
/* FUNCAO: table15_geom                                                       */
/* DESCRICAO: Retorna geometria usada para reproduzir a Tabela 15 (cavidade   */
/* esferica de raio unitario em cm).                                          */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: SphereGeom.                                                         */
/******************************************************************************/
inline SphereGeom table15_geom() { return {1.0}; } // radius 1 cm

/******************************************************************************/
/* FUNCAO: default_grid_fig15                                                 */
/* DESCRICAO: Define discretizacao padrao para a Figura 15, equilibrando      */
/* custo computacional e erro nas primeiras frequencias.                      */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: Grid3D.                                                             */
/******************************************************************************/
inline Grid3D default_grid_fig15() { return {6, 3, 3}; }
/******************************************************************************/
/* FUNCAO: default_grid_fig16                                                 */
/* DESCRICAO: Retorna malha padrao para o caso da Figura 16 (cavidade meio preenchida). */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: Grid3D.                                                             */
/******************************************************************************/
inline Grid3D default_grid_fig16() { return {5, 5, 4}; }
/******************************************************************************/
/* FUNCAO: default_grid_fig17                                                 */
/* DESCRICAO: Retorna malha padrao para o caso da Figura 17 (cilindro circular 3D).   */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: Grid3D.                                                             */
/******************************************************************************/
inline Grid3D default_grid_fig17() { return {7, 7, 4}; }
/******************************************************************************/
/* FUNCAO: default_grid_table15                                               */
/* DESCRICAO: Retorna discretizacao recomendada para reproduzir a Tabela 15 (cavidade esferica). */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: Grid3D.                                                             */
/******************************************************************************/
inline Grid3D default_grid_table15() { return {6, 5, 5}; }

/******************************************************************************/
/* FUNCAO: table12_rows                                                       */
/* DESCRICAO: Fornece linhas de referencia da Tabela 12 (cavidade retangular  */
/* sem dielestrico) para comparacao FEM3D0/FEM3D1.                            */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::vector<RefRow>.                                                */
/******************************************************************************/
inline std::vector<RefRow> table12_rows()
{
  return {
      {"TE101", 5.236, 5.213},
      {"TM110", 7.025, 6.977},
      {"TE011", 7.531, 7.474},
      {"TE201", 7.531, 7.573},
      {"TE111", 8.179, 7.991},
      {"TM111", 8.179, 8.122},
      {"TM210", 8.886, 8.572},
      {"TE102", 8.947, 8.795},
  };
}

/******************************************************************************/
/* FUNCAO: table13_rows                                                       */
/* DESCRICAO: Fornece linhas de referencia da Tabela 13 para validacao automatica do caso 3D meio preenchido. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::vector<RefRow>.                                                */
/******************************************************************************/
inline std::vector<RefRow> table13_rows()
{
  return {
      {"TE2101", 3.538, 3.534},
      {"TE2201", 5.445, 5.440},
      {"TE3102", 5.935, 5.916},
      {"TE2301", 7.503, 7.501},
      {"TE2202", 7.633, 7.560},
      {"TE2103", 8.096, 8.056},
  };
}

/******************************************************************************/
/* FUNCAO: table14_rows                                                       */
/* DESCRICAO: Fornece linhas de referencia da Tabela 14 para validacao automatica do cilindro circular. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::vector<RefRow>.                                                */
/******************************************************************************/
inline std::vector<RefRow> table14_rows()
{
  return {
      {"TM010", 4.810, 4.809},
      {"TE111a", 7.283, 7.202},
      {"TE111b", 7.283, 7.288},
      {"TM110a", 7.650, 7.633},
      {"TM110b", 7.650, 7.724},
      {"TM011", 7.840, 7.940},
      {"TE211a", 8.658, 8.697},
  };
}

/******************************************************************************/
/* FUNCAO: table15_rows                                                       */
/* DESCRICAO: Fornece linhas de referencia da Tabela 15 para validacao automatica da cavidade esferica. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::vector<RefRow>.                                                */
/******************************************************************************/
inline std::vector<RefRow> table15_rows()
{
  return {
      {"TM010", 2.744, 2.799},
      {"TM001", 2.744, 2.802},
      {"TM100", 2.744, 2.811},
      {"TM021", 3.870, 3.948},
      {"TM121e", 3.870, 3.986},
      {"TM121o", 3.870, 3.994},
      {"TM221e", 3.870, 4.038},
      {"TM221o", 3.870, 4.048},
      {"TE001", 4.493, 4.433},
      {"TE111e", 4.493, 4.472},
      {"TE111o", 4.493, 4.549},
  };
}

/******************************************************************************/
/* FUNCAO: scan_limit_for_table                                               */
/* DESCRICAO: Define limite superior de busca espectral adequado para cada tabela de referencia. */
/* ENTRADA: nrows: int.                                                       */
/* SAIDA: int.                                                                */
/******************************************************************************/
inline int scan_limit_for_table(int nrows)
{
  return std::max(100, nrows * 12);
}

} // namespace fem3d
