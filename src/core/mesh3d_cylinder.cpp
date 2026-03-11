/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_cylinder.cpp                                      */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica para cavidade cilindrica circular  */
/* por recorte cartesiano (Figura 17 / Tabela 14).                            */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "mesh3d_cylinder.hpp"
#include "mesh3d_rect.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

/******************************************************************************/
/* FUNCAO: make_cylinder_tet_mesh_cartesian                                   */
/* DESCRICAO: Gera malha tetraedrica recortada para cavidade cilindrica.      */
/* Estrategia: gera caixa envolvente e mantem tets cujo centroide pertence    */
/* ao cilindro.                                                               */
/* ENTRADA: radius: double; height: double; nx: int; ny: int; nz: int.        */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_cylinder_tet_mesh_cartesian(
    double radius,
    double height,
    int nx,
    int ny,
    int nz)
{
  if (!(radius > 0.0 && height > 0.0))
    throw std::runtime_error("radius e height devem ser > 0.");
  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::runtime_error("nx, ny, nz devem ser > 0.");

  // Passo 1: inicia com uma malha de caixa estruturada envolvendo o cilindro.
  Mesh3D box = make_rect_tet_mesh(2.0 * radius, 2.0 * radius, height, nx, ny, nz);

  for (auto &n : box.nodes)
  {
    n.x -= radius;
    n.y -= radius;
  }

  // Passo 2: mantem tetraedros cujo centroide fica na secao circular do cilindro.
  // Esse recorte por centroide e uma aproximacao simples/robusta usada apenas
  // na geracao geometrica; a formulacao de FEM permanece inalterada.
  const double tol = 1e-12;
  std::vector<int> keep_tet(box.tets.size(), 0);
  std::vector<int> node_map(box.nodes.size(), -1);

  for (int tid = 0; tid < (int)box.tets.size(); ++tid)
  {
    const Tet &t = box.tets[(size_t)tid];
    const Node3D &n0 = box.nodes[(size_t)t.v[0]];
    const Node3D &n1 = box.nodes[(size_t)t.v[1]];
    const Node3D &n2 = box.nodes[(size_t)t.v[2]];
    const Node3D &n3 = box.nodes[(size_t)t.v[3]];
    const double xc = 0.25 * (n0.x + n1.x + n2.x + n3.x);
    const double yc = 0.25 * (n0.y + n1.y + n2.y + n3.y);
    const double rc2 = xc * xc + yc * yc;
    if (rc2 <= radius * radius + tol)
      keep_tet[(size_t)tid] = 1;
  }

  Mesh3D out;
  out.tets.reserve(box.tets.size());
  out.nodes.reserve(box.nodes.size());

  // Passo 3: reconstrucao de numeracao compacta e marcacao de nos de contorno.
  // Para a parede curva, usamos tolerancia dependente do tamanho da malha.
  const double dx = 2.0 * radius / (double)nx;
  const double dy = 2.0 * radius / (double)ny;
  const double dr = std::sqrt(dx * dx + dy * dy);
  const double bd_tol = 0.75 * dr;

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
        const double r = std::sqrt(n.x * n.x + n.y * n.y);
        const bool on_cap = (std::abs(n.z) < 1e-12) || (std::abs(n.z - height) < 1e-12);
        const bool on_wall = (std::abs(r - radius) <= bd_tol) || (r > radius - bd_tol);
        Node3D nn = n;
        nn.is_boundary = on_cap || on_wall;
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
