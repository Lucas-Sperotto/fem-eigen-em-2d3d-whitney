/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_coax.cpp                                          */
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
/* FUNCAO: make_coax_mesh                                                     */
/* DESCRICAO: Gera malha triangular para dominio anular coaxial.              */
/* ENTRADA: r1: double; r2: double; nr: int; nt: int.                         */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_coax_mesh(double r1, double r2, int nr, int nt)
{
    Mesh2D m;

    // nr camadas radiais entre r1 e r2
    // nt divisões angulares

    auto ring_id = [&](int k, int t)
    {
        return k * nt + t;
    };

    // nós
    for (int k = 0; k <= nr; k++)
    {
        double r = r1 + (r2 - r1) * (double)k / (double)nr;
        for (int t = 0; t < nt; t++)
        {
            double th = 2.0 * M_PI * (double)t / (double)nt;
            double x = r * std::cos(th);
            double y = r * std::sin(th);

            bool boundary = (k == 0 || k == nr);
            m.nodes.push_back({x, y, boundary});
        }
    }

    // triângulos
    for (int k = 0; k < nr; k++)
    {
        for (int t = 0; t < nt; t++)
        {
            int t2 = (t + 1) % nt;

            int a = ring_id(k, t);
            int b = ring_id(k, t2);
            int c = ring_id(k + 1, t);
            int d = ring_id(k + 1, t2);

            m.tris.push_back({{a, c, d}});
            m.tris.push_back({{a, d, b}});
        }
    }

    return m;
}
