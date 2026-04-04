/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_rect_edge_match_wrapper.hpp                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Wrapper local para evitar colisao de simbolos entre os matchers */
/* escalares e de aresta do caso retangular.                                  */
/*****************************************************************************/

#pragma once

#define kc_rect_analytic helmvec1_edge_rect_kc_rect_analytic
#define mass_inner helmvec1_edge_rect_mass_inner
#define mass_norm helmvec1_edge_rect_mass_norm
#define mass_correlation_abs helmvec1_edge_rect_mass_correlation_abs
#include "edge/mode_match_rect_edge.hpp"
#undef mass_correlation_abs
#undef mass_norm
#undef mass_inner
#undef kc_rect_analytic
