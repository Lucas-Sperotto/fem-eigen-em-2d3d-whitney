/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_sphere.cpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica para cavidade esferica por recorte */
/* cartesiano (Tabela 15).                                                    */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "mesh3d_sphere.hpp"
#include "mesh3d_rect.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

/******************************************************************************/
/* FUNCAO: make_sphere_tet_mesh_cartesian                                     */
/* DESCRICAO: Gera malha tetraedrica recortada para cavidade esferica.        */
/* Estrategia: caixa envolvente + recorte por centroide + marcacao de         */
/* contorno aproximada pela distancia ao raio da esfera.                      */
/* ENTRADA: radius: double; nx: int; ny: int; nz: int.                        */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_sphere_tet_mesh_cartesian(
    double radius,
    int nx,
    int ny,
    int nz)
{
  if (!(radius > 0.0))
    throw std::runtime_error("radius deve ser > 0.");
  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::runtime_error("nx, ny, nz devem ser > 0.");

  // Passo 1: monta uma malha de caixa estruturada envolvendo a esfera.
  Mesh3D box = make_rect_tet_mesh(2.0 * radius, 2.0 * radius, 2.0 * radius, nx, ny, nz);
  for (auto &n : box.nodes)
  {
    n.x -= radius;
    n.y -= radius;
    n.z -= radius;
  }

  // Passo 2: recorte por centroide para manter tetraedros dentro da esfera.
  // Essa aproximacao geometrica e suficiente para os estudos comparativos da
  // Tabela 15, mantendo simplicidade de implementacao.
  std::vector<int> keep_tet(box.tets.size(), 0);
  std::vector<int> node_map(box.nodes.size(), -1);
  const double r2 = radius * radius;
  const double tol = 1e-12;

  for (int tid = 0; tid < (int)box.tets.size(); ++tid)
  {
    const Tet &t = box.tets[(size_t)tid];
    const Node3D &n0 = box.nodes[(size_t)t.v[0]];
    const Node3D &n1 = box.nodes[(size_t)t.v[1]];
    const Node3D &n2 = box.nodes[(size_t)t.v[2]];
    const Node3D &n3 = box.nodes[(size_t)t.v[3]];
    const double xc = 0.25 * (n0.x + n1.x + n2.x + n3.x);
    const double yc = 0.25 * (n0.y + n1.y + n2.y + n3.y);
    const double zc = 0.25 * (n0.z + n1.z + n2.z + n3.z);
    if (xc * xc + yc * yc + zc * zc <= r2 + tol)
      keep_tet[(size_t)tid] = 1;
  }

  Mesh3D out;
  out.nodes.reserve(box.nodes.size());
  out.tets.reserve(box.tets.size());

  // Passo 3: renumeracao compacta + marcacao de contorno.
  // O contorno curvo e aproximado selecionando nos proximos ao raio da esfera.
  const double hx = 2.0 * radius / (double)nx;
  const double hy = 2.0 * radius / (double)ny;
  const double hz = 2.0 * radius / (double)nz;
  const double h = std::max({hx, hy, hz});
  const double bd_tol = 0.90 * h;

  for (int tid = 0; tid < (int)box.tets.size(); ++tid)
  {
    if (!keep_tet[(size_t)tid])
      continue;

    Tet tnew = {{-1, -1, -1, -1}};
    const Tet &t = box.tets[(size_t)tid];
    for (int lv = 0; lv < 4; ++lv)
    {
      const int vold = t.v[lv];
      int vnew = node_map[(size_t)vold];
      if (vnew < 0)
      {
        const Node3D &n = box.nodes[(size_t)vold];
        const double rn = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
        Node3D nn = n;
        nn.is_boundary = (std::abs(rn - radius) <= bd_tol) || (rn > radius - bd_tol);
        vnew = (int)out.nodes.size();
        out.nodes.push_back(nn);
        node_map[(size_t)vold] = vnew;
      }
      tnew.v[lv] = vnew;
    }
    out.tets.push_back(tnew);
  }

  return out;
}
