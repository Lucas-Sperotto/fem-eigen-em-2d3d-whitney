/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_dofs.cpp                                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 2D de elementos de aresta (DOFs, base de Whitney e       */
/* montagem).                                                                 */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "edge_dofs.hpp"

struct PairHash
{
/******************************************************************************/
/* FUNCAO: operator                                                           */
/* DESCRICAO: Implementa operador de hash/comparacao auxiliar para indexacao eficiente de entidades topologicas. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::size_t.                                                        */
/******************************************************************************/
    std::size_t operator()(const std::uint64_t &k) const noexcept { return (std::size_t)k; }
};

/******************************************************************************/
/* FUNCAO: key_undirected                                                     */
/* DESCRICAO: Gera chave canonica (min,max) para identificar aresta sem       */
/* orientacao e evitar duplicacao na montagem de DOFs globais.                */
/* ENTRADA: a: int; b: int.                                                   */
/* SAIDA: std::uint64_t.                                                      */
/******************************************************************************/
static std::uint64_t key_undirected(int a, int b)
{
    int i = (a < b) ? a : b;
    int j = (a < b) ? b : a;
    return ((std::uint64_t)(uint32_t)i << 32) | (uint32_t)j;
}

/******************************************************************************/
/* FUNCAO: build_edge_dofs                                                    */
/* DESCRICAO: Monta e retorna a estrutura principal do modulo, combinando     */
/* malha, materiais e condicoes de contorno.                                  */
/* ENTRADA: m: const Mesh2D &.                                                */
/* SAIDA: EdgeDofs.                                                           */
/******************************************************************************/
EdgeDofs build_edge_dofs(const Mesh2D &m)
{
    EdgeDofs out;
    out.tri_edges.resize(m.tris.size());

    std::unordered_map<std::uint64_t, int, PairHash> map;
    map.reserve(m.tris.size() * 3);

    auto add_edge = [&](int tri_id, int a, int b) -> std::pair<int, int>
    {
        int g0 = (a < b) ? a : b;
        int g1 = (a < b) ? b : a;
        int sgn = (a < b) ? +1 : -1; // orientacao local (a->b) vs global (min->max)

        std::uint64_t k = key_undirected(a, b);
        auto it = map.find(k);
        if (it == map.end())
        {
            int id = (int)out.edges.size();
            map[k] = id;
            Edge e;
            e.n0 = g0;
            e.n1 = g1;
            e.triL = tri_id;
            out.edges.push_back(e);
            return {id, sgn};
        }
        else
        {
            int id = it->second;
            if (out.edges[id].triR < 0)
                out.edges[id].triR = tri_id;
            return {id, sgn};
        }
    };

    for (int tid = 0; tid < (int)m.tris.size(); tid++)
    {
        const Tri &t = m.tris[tid];
        int a = t.v[0], b = t.v[1], c = t.v[2];

        auto e0 = add_edge(tid, a, b);
        auto e1 = add_edge(tid, b, c);
        auto e2 = add_edge(tid, c, a);

        out.tri_edges[tid].e = {e0.first, e1.first, e2.first};
        out.tri_edges[tid].sgn = {e0.second, e1.second, e2.second};
    }

    for (auto &e : out.edges)
    {
        e.is_boundary = (e.triR < 0); // 1 tri => contorno
    }

    out.edge_to_dof.assign(out.edges.size(), -1);
    int cnt = 0;
    for (int i = 0; i < (int)out.edges.size(); i++)
    {
        if (!out.edges[i].is_boundary)
        {
            out.edge_to_dof[i] = cnt++;
        }
    }
    out.ndof = cnt;
    return out;
}
/******************************************************************************/
/* FUNCAO: build_edge_dofs                                                    */
/* DESCRICAO: Monta e retorna a estrutura principal do modulo, combinando     */
/* malha, materiais e condicoes de contorno.                                  */
/* ENTRADA: m: const Mesh2D &; bc: EdgeBC.                                    */
/* SAIDA: EdgeDofs.                                                           */
/******************************************************************************/
EdgeDofs build_edge_dofs(const Mesh2D &m, EdgeBC bc)
{
    EdgeDofs out;
    out.tri_edges.resize(m.tris.size());

    std::unordered_map<std::uint64_t, int, PairHash> map;
    map.reserve(m.tris.size() * 3);

    auto add_edge = [&](int tri_id, int a, int b) -> std::pair<int, int>
    {
        int g0 = (a < b) ? a : b;
        int g1 = (a < b) ? b : a;
        int sgn = (a < b) ? +1 : -1; // orientacao local (a->b) vs global (min->max)

        std::uint64_t k = key_undirected(a, b);
        auto it = map.find(k);
        if (it == map.end())
        {
            int id = (int)out.edges.size();
            map[k] = id;
            Edge e;
            e.n0 = g0;
            e.n1 = g1;
            e.triL = tri_id;
            out.edges.push_back(e);
            return {id, sgn};
        }
        else
        {
            int id = it->second;
            if (out.edges[id].triR < 0)
                out.edges[id].triR = tri_id;
            return {id, sgn};
        }
    };

    for (int tid = 0; tid < (int)m.tris.size(); tid++)
    {
        const Tri &t = m.tris[tid];
        int a = t.v[0], b = t.v[1], c = t.v[2];

        auto e0 = add_edge(tid, a, b);
        auto e1 = add_edge(tid, b, c);
        auto e2 = add_edge(tid, c, a);

        out.tri_edges[tid].e = {e0.first, e1.first, e2.first};
        out.tri_edges[tid].sgn = {e0.second, e1.second, e2.second};
    }

    for (auto &e : out.edges)
    {
        e.is_boundary = (e.triR < 0); // 1 tri => contorno
    }

    out.edge_to_dof.assign(out.edges.size(), -1);
    int cnt = 0;
    for (int i = 0; i < (int)out.edges.size(); i++)
    {
        bool eliminate = (bc == EdgeBC::TE_PEC_TangentialZero) && out.edges[i].is_boundary;
        if (!eliminate)
        {
            out.edge_to_dof[i] = cnt++;
        }
    }
    out.ndof = cnt;

    return out;
}
