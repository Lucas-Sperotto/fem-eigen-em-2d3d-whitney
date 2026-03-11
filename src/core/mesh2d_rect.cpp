/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_rect.cpp                                          */
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

#include "mesh2d.hpp"
#include <cmath>

/******************************************************************************/
/* FUNCAO: make_rect_mesh                                                     */
/* DESCRICAO: Gera malha triangular estruturada para dominio retangular.      */
/* ENTRADA: a: double; b: double; nx: int; ny: int.                           */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_rect_mesh(double a, double b, int nx, int ny)
{
  Mesh2D m;
  m.nodes.reserve((size_t)(nx + 1) * (ny + 1));
  m.tris.reserve((size_t)2 * nx * ny);

  auto node_id = [&](int i, int j)
  { return j * (nx + 1) + i; };

  for (int j = 0; j <= ny; j++)
  {
    double y = b * (double)j / (double)ny;
    for (int i = 0; i <= nx; i++)
    {
      double x = a * (double)i / (double)nx;
      bool bd = (i == 0 || i == nx || j == 0 || j == ny);
      m.nodes.push_back({x, y, bd});
    }
  }

  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      int n00 = node_id(i, j);
      int n10 = node_id(i + 1, j);
      int n01 = node_id(i, j + 1);
      int n11 = node_id(i + 1, j + 1);

      // diagonal n00 -> n11
      m.tris.push_back({{n00, n10, n11}});
      m.tris.push_back({{n00, n11, n01}});
    }
  }
  return m;
}
