/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_rect.hpp                                          */
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

// Gera malha retangular W x H com nx x ny divisões
/******************************************************************************/
/* FUNCAO: make_rect_mesh                                                     */
/* DESCRICAO: Gera malha triangular estruturada para dominio retangular.      */
/* ENTRADA: a: double; b: double; nx: int; ny: int.                           */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_rect_mesh(double a, double b, int nx, int ny);
