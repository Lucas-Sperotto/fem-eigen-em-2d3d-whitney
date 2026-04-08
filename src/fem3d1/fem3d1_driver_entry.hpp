/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d1/fem3d1_driver_entry.hpp                                */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Entradas publicas compartilhadas do solver FEM3D1.              */
/*****************************************************************************/

#pragma once

#include "fem3d/fem3d_case_driver.hpp"

int run_fem3d1_case_driver(int argc, char **argv, fem3d::CaseId fixed_case, const char *bin_name);
