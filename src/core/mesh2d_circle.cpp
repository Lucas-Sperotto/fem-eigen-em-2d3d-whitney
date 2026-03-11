/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh2d_circle.cpp                                        */
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
/* FUNCAO: make_circle_mesh                                                   */
/* DESCRICAO: Gera malha triangular para dominio circular por aneis           */
/* concentricos.                                                              */
/* ENTRADA: R: double; nr: int; nt: int.                                      */
/* SAIDA: Mesh2D.                                                             */
/******************************************************************************/
Mesh2D make_circle_mesh(double R, int nr, int nt)
{
    Mesh2D m;

    // nós: centro + (nr anéis) * (nt pontos por anel)
    // anel k=1..nr, raio = k*R/nr
    // centro = id 0
    m.nodes.reserve(1 + (size_t)nr * nt);

    // centro
    m.nodes.push_back({0.0, 0.0, false});

    auto ring_id = [&](int k, int t) -> int
    {
        // k in [1..nr], t in [0..nt-1]
        return 1 + (k - 1) * nt + t;
    };

    for (int k = 1; k <= nr; k++)
    {
        double r = R * (double)k / (double)nr;
        for (int t = 0; t < nt; t++)
        {
            double th = 2.0 * M_PI * (double)t / (double)nt;
            double x = r * std::cos(th);
            double y = r * std::sin(th);
            bool bd = (k == nr); // borda externa
            m.nodes.push_back({x, y, bd});
        }
    }

    // Triangulação:
    // 1) setores do centro para o primeiro anel
    for (int t = 0; t < nt; t++)
    {
        int t2 = (t + 1) % nt;
        int v0 = 0;
        int v1 = ring_id(1, t);
        int v2 = ring_id(1, t2);
        m.tris.push_back({{v0, v1, v2}});
    }

    // 2) entre anéis: cada quadrilátero vira 2 triângulos
    for (int k = 1; k < nr; k++)
    {
        for (int t = 0; t < nt; t++)
        {
            int t2 = (t + 1) % nt;
            int a = ring_id(k, t);
            int b = ring_id(k, t2);
            int c = ring_id(k + 1, t);
            int d = ring_id(k + 1, t2);

            // diagonal a->d
            m.tris.push_back({{a, c, d}});
            m.tris.push_back({{a, d, b}});
        }
    }

    return m;
}
