/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_coax.hpp                                          */
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

// Gera malha coaxial entre r1 e r2
/******************************************************************************/
/* FUNCAO: make_coax_mesh                                                     */
/* DESCRICAO: Gera malha triangular para dominio anular coaxial.              */
/* ENTRADA: r1: double; r2: double; nr: int; nt: int.                         */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_coax_mesh(double r1, double r2, int nr, int nt);
