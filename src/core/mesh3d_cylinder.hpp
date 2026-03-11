/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_cylinder.hpp                                      */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica para cavidade cilindrica circular  */
/* (Figura 17 / Tabela 14), via recorte de malha cartesiana.                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh3d.hpp"

// Malha tetraedrica cartesiana recortada para cavidade cilindrica circular centrada em
// (0,0) no plano x-y e com extensao em z no intervalo [0, height].
//
// A caixa de fundo e [-radius, radius] x [-radius, radius] x [0, height].
// Tetraedros sao mantidos quando o centroide fica dentro do cilindro.
/******************************************************************************/
/* FUNCAO: make_cylinder_tet_mesh_cartesian                                   */
/* DESCRICAO: Gera malha tetraedrica recortada para cavidade cilindrica.      */
/* Usada no caso 3D de validacao da Figura 17 (Secao 3.1.5).                  */
/* ENTRADA: radius: double; height: double; nx: int; ny: int; nz: int.        */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_cylinder_tet_mesh_cartesian(
    double radius,
    double height,
    int nx,
    int ny,
    int nz);
