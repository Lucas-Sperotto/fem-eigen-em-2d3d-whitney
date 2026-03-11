/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_rect.hpp                                          */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica retangular para os casos 3D do     */
/* artigo (Figura 15 e Figura 16) e campos de material por tetraedro.         */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh3d.hpp"
#include <vector>

// Malha retangular estruturada tetraedrizada com particionamento globalmente conforme
// de Freudenthal em 6 tetraedros por celula hexaedrica.
/******************************************************************************/
/* FUNCAO: make_rect_tet_mesh                                                 */
/* DESCRICAO: Gera malha tetraedrica estruturada para paralelepipedo          */
/* retangular com particionamento de Freudenthal (6 tets/celula).             */
/* ENTRADA: lx: double; ly: double; lz: double; nx: int; ny: int; nz: int.    */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_rect_tet_mesh(double lx, double ly, double lz, int nx, int ny, int nz);

// Permissividade relativa por tetraedro (constante por partes) via centroide em z.
/******************************************************************************/
/* FUNCAO: make_eps_r_tets_by_z                                               */
/* DESCRICAO: Atribui eps_r por tetraedro usando corte horizontal em z,       */
/* reproduzindo o caso meio preenchido da Figura 16 (Secao 3.1.5).            */
/* ENTRADA: mesh: const Mesh3D &; z_cut: double; eps_r_below: double;         */
/* eps_r_above: double.                                                       */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> make_eps_r_tets_by_z(
    const Mesh3D &mesh,
    double z_cut,
    double eps_r_below,
    double eps_r_above);
