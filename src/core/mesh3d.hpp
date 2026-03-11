/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d.hpp                                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Estruturas base de malha 3D (nos e tetraedros) usadas pelos     */
/* casos de cavidade da Secao 3.1 (FEM3D0/FEM3D1).                            */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include <vector>

struct Node3D
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  // Marcacao de contorno geometrico, usada para impor BC PEC nos DOFs de aresta.
  bool is_boundary = false;
};

struct Tet
{
  int v[4];
};

struct Mesh3D
{
  // Lista global de nos e conectividade tetraedrica.
  std::vector<Node3D> nodes;
  std::vector<Tet> tets;
};
