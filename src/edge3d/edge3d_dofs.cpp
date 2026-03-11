/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_dofs.cpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 3D de elementos de aresta tetraedricos (Whitney 1-form) e*/
/* montagem.                                                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, integrais*/
/* I1..I10.                                                                   */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "edge3d_dofs.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>

namespace
{
constexpr std::array<std::array<int, 2>, 6> kTetLocalEdges = {{
    {{0, 1}},
    {{0, 2}},
    {{0, 3}},
    {{1, 2}},
    {{1, 3}},
    {{2, 3}},
}};

constexpr std::array<std::array<int, 3>, 4> kTetLocalFaces = {{
    {{0, 1, 2}},
    {{0, 1, 3}},
    {{0, 2, 3}},
    {{1, 2, 3}},
}};

/******************************************************************************/
/* FUNCAO: key_edge_undirected                                                */
/* DESCRICAO: Cria chave canonica de aresta nao orientada para mapear graus de liberdade globais. */
/* ENTRADA: a: int; b: int.                                                   */
/* SAIDA: std::uint64_t.                                                      */
/******************************************************************************/
std::uint64_t key_edge_undirected(int a, int b)
{
  const int i = (a < b) ? a : b;
  const int j = (a < b) ? b : a;
  return ((std::uint64_t)(uint32_t)i << 32) | (uint32_t)j;
}

struct FaceKey
{
  int a = 0;
  int b = 0;
  int c = 0;
  bool operator==(const FaceKey &o) const { return a == o.a && b == o.b && c == o.c; }
};

struct FaceKeyHash
{
/******************************************************************************/
/* FUNCAO: operator                                                           */
/* DESCRICAO: Calcula hash de FaceKey para acelerar mapa de contagem de faces  */
/* (identificacao de fronteira).                                              */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::size_t.                                                        */
/******************************************************************************/
  std::size_t operator()(const FaceKey &k) const noexcept
  {
    std::size_t h = (std::size_t)k.a;
    h = h * 1315423911u + (std::size_t)k.b;
    h = h * 1315423911u + (std::size_t)k.c;
    return h;
  }
};

/******************************************************************************/
/* FUNCAO: make_face_key                                                      */
/* DESCRICAO: Canoniza nos de uma face triangular (ordem crescente) para      */
/* detectar fronteira por contagem de incidencias de faces em tetraedros.     */
/* ENTRADA: i: int; j: int; k: int.                                           */
/* SAIDA: FaceKey.                                                            */
/******************************************************************************/
FaceKey make_face_key(int i, int j, int k)
{
  std::array<int, 3> v = {i, j, k};
  std::sort(v.begin(), v.end());
  return {v[0], v[1], v[2]};
}
} // namespace

/******************************************************************************/
/* FUNCAO: build_edge_dofs_3d                                                 */
/* DESCRICAO: Atalho para montagem de DOFs 3D com condicao PEC padrao         */
/* (remove arestas de fronteira).                                             */
/* ENTRADA: mesh: const Mesh3D &.                                             */
/* SAIDA: EdgeDofs3D.                                                         */
/******************************************************************************/
EdgeDofs3D build_edge_dofs_3d(const Mesh3D &mesh)
{
  return build_edge_dofs_3d(mesh, Edge3DBC::PEC_TangentialZero);
}

/******************************************************************************/
/* FUNCAO: build_edge_dofs_3d                                                 */
/* DESCRICAO: Monta o grafo de arestas globais do tetraedro, registra sinais  */
/* locais/global por elemento e aplica eliminacao de contorno conforme BC.    */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC.                               */
/* SAIDA: EdgeDofs3D.                                                         */
/******************************************************************************/
EdgeDofs3D build_edge_dofs_3d(const Mesh3D &mesh, Edge3DBC bc)
{
  EdgeDofs3D out;
  out.tet_edges.resize(mesh.tets.size());

  std::unordered_map<std::uint64_t, int> edge_map;
  edge_map.reserve(mesh.tets.size() * 6);

  std::unordered_map<FaceKey, int, FaceKeyHash> face_count;
  face_count.reserve(mesh.tets.size() * 4);

  auto add_edge = [&](int a, int b) -> std::pair<int, int>
  {
    const int g0 = (a < b) ? a : b;
    const int g1 = (a < b) ? b : a;
    const int sgn = (a < b) ? +1 : -1;
    const std::uint64_t key = key_edge_undirected(a, b);
    auto it = edge_map.find(key);
    if (it == edge_map.end())
    {
      const int id = (int)out.edges.size();
      edge_map[key] = id;
      out.edges.push_back({g0, g1, false});
      return {id, sgn};
    }
    return {it->second, sgn};
  };

  for (int tid = 0; tid < (int)mesh.tets.size(); ++tid)
  {
    const Tet &t = mesh.tets[tid];

    // Cada face aparece 1 vez na fronteira e 2 vezes no interior.
    for (int lf = 0; lf < 4; ++lf)
    {
      const int i = t.v[kTetLocalFaces[lf][0]];
      const int j = t.v[kTetLocalFaces[lf][1]];
      const int k = t.v[kTetLocalFaces[lf][2]];
      face_count[make_face_key(i, j, k)] += 1;
    }

    for (int le = 0; le < 6; ++le)
    {
      const int i = t.v[kTetLocalEdges[le][0]];
      const int j = t.v[kTetLocalEdges[le][1]];
      const auto eg = add_edge(i, j);
      out.tet_edges[(size_t)tid].e[(size_t)le] = eg.first;
      out.tet_edges[(size_t)tid].sgn[(size_t)le] = eg.second;
    }
  }

  for (const auto &kv : face_count)
  {
    if (kv.second != 1)
      continue;
    const FaceKey &f = kv.first;
    const std::uint64_t e0 = key_edge_undirected(f.a, f.b);
    const std::uint64_t e1 = key_edge_undirected(f.a, f.c);
    const std::uint64_t e2 = key_edge_undirected(f.b, f.c);
    auto it0 = edge_map.find(e0);
    auto it1 = edge_map.find(e1);
    auto it2 = edge_map.find(e2);
    if (it0 != edge_map.end())
      out.edges[(size_t)it0->second].is_boundary = true;
    if (it1 != edge_map.end())
      out.edges[(size_t)it1->second].is_boundary = true;
    if (it2 != edge_map.end())
      out.edges[(size_t)it2->second].is_boundary = true;
  }

  // Mapeia arestas ativas para DOFs; no caso PEC elimina arestas de fronteira
  // para impor Et tangencial nulo no contorno condutor.
  out.edge_to_dof.assign(out.edges.size(), -1);
  int nd = 0;
  for (int eid = 0; eid < (int)out.edges.size(); ++eid)
  {
    const bool eliminate =
        (bc == Edge3DBC::PEC_TangentialZero) &&
        out.edges[(size_t)eid].is_boundary;
    if (!eliminate)
      out.edge_to_dof[(size_t)eid] = nd++;
  }
  out.ndof = nd;
  return out;
}
