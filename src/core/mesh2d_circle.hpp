/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_circle.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo escalar 2D (malha, montagem e identificacao modal).      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh2d.hpp"

// Gera malha circular de raio R
/******************************************************************************/
/* FUNCAO: make_circle_mesh                                                   */
/* DESCRICAO: Gera malha triangular para dominio circular por aneis           */
/* concentricos.                                                              */
/* ENTRADA: R: double; nr: int; nt: int.                                      */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_circle_mesh(double R, int nr, int nt);
