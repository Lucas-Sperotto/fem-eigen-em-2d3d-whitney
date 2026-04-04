/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_circle_edge_match_wrapper.hpp                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Wrapper local para evitar colisao de simbolos entre os matchers */
/* escalares e de aresta do caso circular.                                    */
/*****************************************************************************/

#pragma once

#define Jm helmvec1_edge_circle_Jm
#define Jm_prime helmvec1_edge_circle_Jm_prime
#define find_root_bisection helmvec1_edge_circle_find_root_bisection
#define bessel_roots helmvec1_edge_circle_bessel_roots
#define mass_inner helmvec1_edge_circle_mass_inner
#define mass_norm helmvec1_edge_circle_mass_norm
#define mass_correlation_abs helmvec1_edge_circle_mass_correlation_abs
#include "edge/mode_match_circle_edge.hpp"
#undef mass_correlation_abs
#undef mass_norm
#undef mass_inner
#undef bessel_roots
#undef find_root_bisection
#undef Jm_prime
#undef Jm
