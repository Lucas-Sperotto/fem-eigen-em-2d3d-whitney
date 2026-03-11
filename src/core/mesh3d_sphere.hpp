/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_sphere.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica para cavidade esferica (Tabela 15),*/
/* via recorte de malha cartesiana.                                           */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh3d.hpp"

// Malha tetraedrica cartesiana recortada para cavidade esferica centrada na origem.
//
// A caixa de fundo e [-radius, radius]^3 e tetraedros sao mantidos quando
// seu centroide pertence ao interior da esfera.
/******************************************************************************/
/* FUNCAO: make_sphere_tet_mesh_cartesian                                     */
/* DESCRICAO: Gera malha tetraedrica recortada para cavidade esferica.        */
/* Usada no caso 3D de validacao da Tabela 15 (Secao 3.1.5).                  */
/* ENTRADA: radius: double; nx: int; ny: int; nz: int.                        */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_sphere_tet_mesh_cartesian(
    double radius,
    int nx,
    int ny,
    int nz);
