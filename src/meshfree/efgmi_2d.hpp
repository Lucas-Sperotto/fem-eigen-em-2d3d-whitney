/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/meshfree/efgmi_2d.hpp                                         */
/* Autor: Codex / OpenAI                                                      */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Interpolantes 2D EFGMI por consistencia para montagem escalar,  */
/* vetorial e acoplada em malhas triangulares usadas apenas como fundo de     */
/* integracao.                                                                */
/*****************************************************************************/
/* Observacao: O modulo usa reproducao linear via matriz de momentos (ordem   */
/* 1) e uma base vetorial do tipo Whitney-like formada a partir dos           */
/* interpolantes nodais EFGMI.                                                */
/*****************************************************************************/

#pragma once

#include "core/mesh2d.hpp"
#include "edge/edge_basis.hpp"
#include "edge/edge_dofs.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace efgmi2d
{

struct Mat3
{
    double a[3][3] = {{0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0}};
};

struct ShapeSample
{
    std::vector<int> node_ids;
    std::vector<double> phi;
    std::vector<double> dphidx;
    std::vector<double> dphidy;
};

struct EdgeBasisSample
{
    std::vector<int> dof_ids;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> curl;
};

struct ScalarFieldSample
{
    double value = 0.0;
    double dx = 0.0;
    double dy = 0.0;
};

struct GlobalEdgeBasisValue
{
    double vx = 0.0;
    double vy = 0.0;
    double curl = 0.0;
};

struct LocalTriangleShapeSample
{
    std::array<double, 3> phi = {0.0, 0.0, 0.0};
    std::array<double, 3> dphidx = {0.0, 0.0, 0.0};
    std::array<double, 3> dphidy = {0.0, 0.0, 0.0};
};

struct LocalTriangleEdgeBasisSample
{
    std::array<double, 3> vx = {0.0, 0.0, 0.0};
    std::array<double, 3> vy = {0.0, 0.0, 0.0};
    std::array<double, 3> curl = {0.0, 0.0, 0.0};
};

struct LocalTriangleConsistency
{
    Mat3 nodal_transform{};
    Mat3 edge_transform{};
    double nodal_cond_est = std::numeric_limits<double>::infinity();
    double edge_cond_est = std::numeric_limits<double>::infinity();
    bool nodal_ok = false;
    bool edge_ok = false;
};

struct EFGMIContext2D
{
    const Mesh2D *mesh = nullptr;
    std::vector<double> nodal_length_scale;
    double support_scale = 2.5;
    double interface_support_scale = 4.0;
    double weight_epsilon = 1.0e-3;
    double snap_tolerance = 1.0e-12;
    int min_support = 6;
    int interface_min_support = 9;
    int max_expand_steps = 8;
    bool interface_truncation = false;
    std::vector<int> tri_material_component;
    std::vector<std::vector<int>> node_material_components;
    std::vector<unsigned char> tri_is_interface;
};

struct EFGMIEdgeContext2D
{
    EFGMIContext2D nodal;
    const EdgeDofs *edge_dofs = nullptr;
    std::vector<std::vector<int>> node_to_edges;
    std::vector<std::vector<int>> edge_support_nodes;
    std::vector<int> edge_min_support;
    std::vector<std::vector<int>> edge_patch_nodes;
    std::vector<std::vector<int>> edge_patch_edges;
    std::vector<std::vector<double>> edge_patch_coeffs;
    std::vector<double> edge_patch_cond_est;
    std::vector<unsigned char> edge_patch_ok;
    std::vector<Mat3> tri_edge_transform;
    std::vector<double> tri_edge_cond_est;
    std::vector<unsigned char> tri_edge_ok;
};

struct TriQuadPoint
{
    std::array<double, 3> lambda{};
    double weight = 0.0;
};

constexpr std::array<TriQuadPoint, 7> kTriQuadP5 = {{
    {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}}, 0.2250000000000000},
    {{{0.4701420641051151, 0.4701420641051151, 0.0597158717897698}}, 0.1323941527885062},
    {{{0.4701420641051151, 0.0597158717897698, 0.4701420641051151}}, 0.1323941527885062},
    {{{0.0597158717897698, 0.4701420641051151, 0.4701420641051151}}, 0.1323941527885062},
    {{{0.1012865073234563, 0.1012865073234563, 0.7974269853530872}}, 0.1259391805448272},
    {{{0.1012865073234563, 0.7974269853530872, 0.1012865073234563}}, 0.1259391805448272},
    {{{0.7974269853530872, 0.1012865073234563, 0.1012865073234563}}, 0.1259391805448272},
}};

struct LineQuadPoint
{
    double s = 0.0;
    double weight = 0.0;
};

// Regra de Gauss-Legendre com 4 pontos mapeada para [0,1] (grau 7).
constexpr std::array<LineQuadPoint, 4> kLineQuadP4 = {{
    {0.0694318442029737, 0.1739274225687269},
    {0.3300094782075719, 0.3260725774312731},
    {0.6699905217924281, 0.3260725774312731},
    {0.9305681557970262, 0.1739274225687269},
}};

inline ShapeSample evaluate_shape(
    const EFGMIContext2D &ctx,
    double x,
    double y,
    bool snap_to_node = false,
    int required_min_support = -1,
    int reference_tri_id = -1);

inline ShapeSample evaluate_shape_on_candidates(
    const EFGMIContext2D &ctx,
    const std::vector<int> &candidate_node_ids,
    double x,
    double y,
    bool snap_to_node = false,
    int required_min_support = -1,
    int reference_tri_id = -1);

inline GlobalEdgeBasisValue evaluate_raw_global_edge_basis_from_sample(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    const ShapeSample &sample);

inline bool build_edge_patch_correction(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    std::vector<double> &coeffs,
    double &cond_est);

inline GlobalEdgeBasisValue evaluate_global_edge_basis(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    double x,
    double y);

inline bool build_triangle_edge_transform(
    const EFGMIEdgeContext2D &ctx,
    int tri_id,
    Mat3 &transform,
    double &cond_est);

inline LocalTriangleEdgeBasisSample apply_local_triangle_edge_transform(
    const LocalTriangleEdgeBasisSample &raw,
    const Mat3 &transform);

inline LocalTriangleEdgeBasisSample evaluate_triangle_edge_basis(
    const EFGMIEdgeContext2D &ctx,
    int tri_id,
    double x,
    double y);

namespace detail
{

constexpr double kConsistencyCondLimit = 1.0e6;

inline Mat3 identity3()
{
    Mat3 out;
    out.a[0][0] = 1.0;
    out.a[1][1] = 1.0;
    out.a[2][2] = 1.0;
    return out;
}

inline std::array<double, 3> mat_vec(const Mat3 &m, const std::array<double, 3> &v)
{
    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out[i] += m.a[i][j] * v[j];
    return out;
}

inline Mat3 mat_mul(const Mat3 &lhs, const Mat3 &rhs)
{
    Mat3 out;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            out.a[i][j] = 0.0;
            for (int k = 0; k < 3; ++k)
                out.a[i][j] += lhs.a[i][k] * rhs.a[k][j];
        }
    }
    return out;
}

inline double dot3(const std::array<double, 3> &lhs, const std::array<double, 3> &rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

inline double det3(const Mat3 &m)
{
    return
        m.a[0][0] * (m.a[1][1] * m.a[2][2] - m.a[1][2] * m.a[2][1]) -
        m.a[0][1] * (m.a[1][0] * m.a[2][2] - m.a[1][2] * m.a[2][0]) +
        m.a[0][2] * (m.a[1][0] * m.a[2][1] - m.a[1][1] * m.a[2][0]);
}

inline double inf_norm3(const Mat3 &m)
{
    double out = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        double row_sum = 0.0;
        for (int j = 0; j < 3; ++j)
            row_sum += std::abs(m.a[i][j]);
        out = std::max(out, row_sum);
    }
    return out;
}

inline bool inverse3(const Mat3 &m, Mat3 &inv_out, double &cond_est)
{
    const double det = det3(m);
    if (!std::isfinite(det) || std::abs(det) <= 1.0e-18)
        return false;

    inv_out.a[0][0] = +(m.a[1][1] * m.a[2][2] - m.a[1][2] * m.a[2][1]) / det;
    inv_out.a[0][1] = -(m.a[0][1] * m.a[2][2] - m.a[0][2] * m.a[2][1]) / det;
    inv_out.a[0][2] = +(m.a[0][1] * m.a[1][2] - m.a[0][2] * m.a[1][1]) / det;
    inv_out.a[1][0] = -(m.a[1][0] * m.a[2][2] - m.a[1][2] * m.a[2][0]) / det;
    inv_out.a[1][1] = +(m.a[0][0] * m.a[2][2] - m.a[0][2] * m.a[2][0]) / det;
    inv_out.a[1][2] = -(m.a[0][0] * m.a[1][2] - m.a[0][2] * m.a[1][0]) / det;
    inv_out.a[2][0] = +(m.a[1][0] * m.a[2][1] - m.a[1][1] * m.a[2][0]) / det;
    inv_out.a[2][1] = -(m.a[0][0] * m.a[2][1] - m.a[0][1] * m.a[2][0]) / det;
    inv_out.a[2][2] = +(m.a[0][0] * m.a[1][1] - m.a[0][1] * m.a[1][0]) / det;

    cond_est = inf_norm3(m) * inf_norm3(inv_out);
    return std::isfinite(cond_est);
}

inline void add_outer_product_scaled(
    Mat3 &dst,
    const std::array<double, 3> &p,
    double alpha)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            dst.a[i][j] += alpha * p[i] * p[j];
}

inline void compact_singular_weight(
    double dx,
    double dy,
    double h,
    double eps,
    double &w,
    double &dwdx,
    double &dwdy)
{
    w = 0.0;
    dwdx = 0.0;
    dwdy = 0.0;
    if (!(h > 0.0))
        return;

    const double r2 = dx * dx + dy * dy;
    const double r = std::sqrt(r2);
    const double q = r / h;
    if (q >= 1.0)
        return;

    // Peso singularizado do tipo Schwartz truncado:
    //   w(q) = exp(-1/(1-q^2)) / (q^2 + eps^2),   0 <= q < 1
    //          0,                                  q >= 1
    // Trata-se de uma bump function compacta C-infinita; a singularizacao
    // regularizada preserva a propriedade interpolante sem explodir em q=0.
    const double one_minus_q2 = 1.0 - q * q;
    if (one_minus_q2 <= 0.0)
        return;

    const double eps2 = eps * eps;
    const double denom = q * q + eps2;
    const double bump = std::exp(-1.0 / one_minus_q2);

    w = bump / denom;
    if (r <= 1.0e-18)
        return;

    const double dbump_dq = bump * ((-2.0 * q) / (one_minus_q2 * one_minus_q2));
    const double dw_dq =
        dbump_dq / denom -
        bump * (2.0 * q) / (denom * denom);
    const double dw_dr = dw_dq / h;
    dwdx = dw_dr * dx / r;
    dwdy = dw_dr * dy / r;
}

inline std::vector<double> compute_nodal_length_scale(const Mesh2D &mesh)
{
    std::vector<double> acc(mesh.nodes.size(), 0.0);
    std::vector<int> cnt(mesh.nodes.size(), 0);
    double global_mean = 0.0;
    int global_cnt = 0;

    const auto add_edge = [&](int a, int b)
    {
        const double dx = mesh.nodes[(size_t)a].x - mesh.nodes[(size_t)b].x;
        const double dy = mesh.nodes[(size_t)a].y - mesh.nodes[(size_t)b].y;
        const double len = std::sqrt(dx * dx + dy * dy);
        if (!(len > 0.0))
            return;
        acc[(size_t)a] += len;
        acc[(size_t)b] += len;
        cnt[(size_t)a] += 1;
        cnt[(size_t)b] += 1;
        global_mean += len;
        global_cnt += 1;
    };

    for (const Tri &tri : mesh.tris)
    {
        add_edge(tri.v[0], tri.v[1]);
        add_edge(tri.v[1], tri.v[2]);
        add_edge(tri.v[2], tri.v[0]);
    }

    const double fallback = (global_cnt > 0) ? (global_mean / global_cnt) : 1.0;
    std::vector<double> out(mesh.nodes.size(), fallback);
    for (size_t i = 0; i < mesh.nodes.size(); ++i)
    {
        if (cnt[i] > 0)
            out[i] = acc[i] / static_cast<double>(cnt[i]);
    }
    return out;
}

struct WeightedSupport
{
    std::vector<int> node_ids;
    std::vector<double> w;
    std::vector<double> dwdx;
    std::vector<double> dwdy;
    int compact_count = 0;
    bool used_residual_fallback = false;
};

struct EdgeSupportCoverage
{
    int tangent_neg = 0;
    int tangent_pos = 0;
    int normal_neg = 0;
    int normal_pos = 0;
    int boundary_interior = 0;
};

inline int parse_env_positive_int(const char *name, int fallback)
{
    const char *raw = std::getenv(name);
    if (raw == nullptr || std::string(raw).empty())
        return fallback;
    try
    {
        const int value = std::stoi(std::string(raw));
        if (value <= 0)
            throw std::runtime_error("nao-positivo");
        return value;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(std::string("Valor invalido em ") + name + ": " + raw);
    }
}

inline double parse_env_positive_double(const char *name, double fallback)
{
    const char *raw = std::getenv(name);
    if (raw == nullptr || std::string(raw).empty())
        return fallback;
    try
    {
        const double value = std::stod(std::string(raw));
        if (!(value > 0.0))
            throw std::runtime_error("nao-positivo");
        return value;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(std::string("Valor invalido em ") + name + ": " + raw);
    }
}

struct MaterialKeyBits
{
    std::uint64_t eps_bits = 0;
    std::uint64_t mu_bits = 0;

    bool operator==(const MaterialKeyBits &other) const
    {
        return eps_bits == other.eps_bits && mu_bits == other.mu_bits;
    }
};

struct MaterialKeyBitsHash
{
    std::size_t operator()(const MaterialKeyBits &key) const noexcept
    {
        const std::size_t h1 = std::hash<std::uint64_t>{}(key.eps_bits);
        const std::size_t h2 = std::hash<std::uint64_t>{}(key.mu_bits);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

inline std::uint64_t double_to_bits(double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline MaterialKeyBits make_material_key(double eps_r, double mu_r)
{
    return {double_to_bits(eps_r), double_to_bits(mu_r)};
}

inline std::uint64_t undirected_edge_key(int a, int b)
{
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
}

inline bool tri_is_interface_region(const EFGMIContext2D &ctx, int tri_id)
{
    return tri_id >= 0 &&
           (size_t)tri_id < ctx.tri_is_interface.size() &&
           ctx.tri_is_interface[(size_t)tri_id] != 0;
}

inline int required_min_support_for_triangle(
    const EFGMIContext2D &ctx,
    int tri_id,
    int requested_min_support)
{
    if (requested_min_support > 0)
        return requested_min_support;
    if (tri_is_interface_region(ctx, tri_id))
        return std::max(ctx.min_support, ctx.interface_min_support);
    return ctx.min_support;
}

inline double support_scale_for_triangle(
    const EFGMIContext2D &ctx,
    int tri_id)
{
    if (tri_is_interface_region(ctx, tri_id))
        return std::max(ctx.support_scale, ctx.interface_support_scale);
    return ctx.support_scale;
}

inline bool node_belongs_to_material_component(
    const EFGMIContext2D &ctx,
    int node_id,
    int tri_id)
{
    if (!ctx.interface_truncation)
        return true;
    if (tri_id < 0 || (size_t)tri_id >= ctx.tri_material_component.size())
        return true;
    if (node_id < 0 || (size_t)node_id >= ctx.node_material_components.size())
        return true;

    const int component = ctx.tri_material_component[(size_t)tri_id];
    if (component < 0)
        return true;

    const auto &components = ctx.node_material_components[(size_t)node_id];
    return std::find(components.begin(), components.end(), component) != components.end();
}

inline void build_material_component_metadata(
    EFGMIContext2D &ctx,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    if (ctx.mesh == nullptr)
        throw std::runtime_error("EFGMIContext2D sem mesh associado.");

    const Mesh2D &mesh = *ctx.mesh;
    if (eps_r_tri.size() != mesh.tris.size() || mu_r_tri.size() != mesh.tris.size())
        throw std::runtime_error("EFGMI: perfil material incompatível com mesh.tris.");
    if (mesh.tris.empty())
        return;

    std::vector<MaterialKeyBits> tri_keys(mesh.tris.size());
    std::unordered_map<MaterialKeyBits, int, MaterialKeyBitsHash> unique_keys;
    unique_keys.reserve(mesh.tris.size());
    for (size_t tid = 0; tid < mesh.tris.size(); ++tid)
    {
        tri_keys[tid] = make_material_key(eps_r_tri[tid], mu_r_tri[tid]);
        unique_keys.emplace(tri_keys[tid], 1);
    }
    if (unique_keys.size() <= 1)
        return;

    ctx.interface_truncation = true;
    ctx.tri_material_component.assign(mesh.tris.size(), -1);
    ctx.tri_is_interface.assign(mesh.tris.size(), 0u);
    ctx.node_material_components.assign(mesh.nodes.size(), {});

    std::vector<std::vector<int>> tri_neighbors(mesh.tris.size());
    std::unordered_map<std::uint64_t, int> edge_owner;
    edge_owner.reserve(mesh.tris.size() * 2);

    for (size_t tid = 0; tid < mesh.tris.size(); ++tid)
    {
        const Tri &tri = mesh.tris[tid];
        for (int e = 0; e < 3; ++e)
        {
            const int a = tri.v[e];
            const int b = tri.v[(e + 1) % 3];
            const std::uint64_t key = undirected_edge_key(a, b);
            const auto [it, inserted] = edge_owner.emplace(key, static_cast<int>(tid));
            if (!inserted)
            {
                const int other_tid = it->second;
                tri_neighbors[tid].push_back(other_tid);
                tri_neighbors[(size_t)other_tid].push_back(static_cast<int>(tid));
                if (!(tri_keys[tid] == tri_keys[(size_t)other_tid]))
                {
                    ctx.tri_is_interface[tid] = 1u;
                    ctx.tri_is_interface[(size_t)other_tid] = 1u;
                }
            }
        }
    }

    int component_count = 0;
    std::vector<int> stack;
    stack.reserve(mesh.tris.size());
    for (size_t seed_tid = 0; seed_tid < mesh.tris.size(); ++seed_tid)
    {
        if (ctx.tri_material_component[seed_tid] >= 0)
            continue;
        ctx.tri_material_component[seed_tid] = component_count;
        stack.push_back(static_cast<int>(seed_tid));
        while (!stack.empty())
        {
            const int tid = stack.back();
            stack.pop_back();
            for (const int nbr_tid : tri_neighbors[(size_t)tid])
            {
                if (ctx.tri_material_component[(size_t)nbr_tid] >= 0)
                    continue;
                if (!(tri_keys[(size_t)tid] == tri_keys[(size_t)nbr_tid]))
                    continue;
                ctx.tri_material_component[(size_t)nbr_tid] = component_count;
                stack.push_back(nbr_tid);
            }
        }
        ++component_count;
    }

    for (size_t tid = 0; tid < mesh.tris.size(); ++tid)
    {
        const int component = ctx.tri_material_component[tid];
        const Tri &tri = mesh.tris[tid];
        for (int local = 0; local < 3; ++local)
        {
            auto &node_components = ctx.node_material_components[(size_t)tri.v[local]];
            if (std::find(node_components.begin(), node_components.end(), component) == node_components.end())
                node_components.push_back(component);
        }
    }
}

inline int boundary_interior_sign(
    const Mesh2D &mesh,
    const Edge &edge,
    double mx,
    double my,
    double nx,
    double ny)
{
    const int tri_id = (edge.triL >= 0) ? edge.triL : edge.triR;
    if (tri_id < 0 || tri_id >= (int)mesh.tris.size())
        return 0;

    const Tri &tri = mesh.tris[(size_t)tri_id];
    const Node2D &a = mesh.nodes[(size_t)tri.v[0]];
    const Node2D &b = mesh.nodes[(size_t)tri.v[1]];
    const Node2D &c = mesh.nodes[(size_t)tri.v[2]];
    const double cx = (a.x + b.x + c.x) / 3.0;
    const double cy = (a.y + b.y + c.y) / 3.0;
    const double proj = (cx - mx) * nx + (cy - my) * ny;
    if (proj > 0.0)
        return +1;
    if (proj < 0.0)
        return -1;
    return 0;
}

inline EdgeSupportCoverage evaluate_edge_support_coverage(
    const Mesh2D &mesh,
    const Edge &edge,
    const std::vector<int> &node_ids,
    const std::vector<double> &nodal_length_scale)
{
    EdgeSupportCoverage out;
    const Node2D &n0 = mesh.nodes[(size_t)edge.n0];
    const Node2D &n1 = mesh.nodes[(size_t)edge.n1];
    const double ex = n1.x - n0.x;
    const double ey = n1.y - n0.y;
    const double L = std::sqrt(ex * ex + ey * ey);
    if (!(L > 0.0))
        return out;

    const double tx = ex / L;
    const double ty = ey / L;
    const double nx = -ty;
    const double ny = tx;
    const double mx = 0.5 * (n0.x + n1.x);
    const double my = 0.5 * (n0.y + n1.y);
    const double h_edge =
        0.5 * (nodal_length_scale[(size_t)edge.n0] + nodal_length_scale[(size_t)edge.n1]);
    const double tangent_tol = 0.20 * L;
    const double normal_tol = std::max(1.0e-12, 0.10 * std::max(L, h_edge));
    const int interior_sign = boundary_interior_sign(mesh, edge, mx, my, nx, ny);

    for (const int node_id : node_ids)
    {
        const Node2D &node = mesh.nodes[(size_t)node_id];
        const double rx = node.x - mx;
        const double ry = node.y - my;
        const double tangent_proj = rx * tx + ry * ty;
        const double normal_proj = rx * nx + ry * ny;

        if (tangent_proj <= -tangent_tol)
            out.tangent_neg += 1;
        if (tangent_proj >= tangent_tol)
            out.tangent_pos += 1;
        if (normal_proj <= -normal_tol)
            out.normal_neg += 1;
        if (normal_proj >= normal_tol)
            out.normal_pos += 1;
        if (interior_sign != 0 && interior_sign * normal_proj >= normal_tol)
            out.boundary_interior += 1;
    }

    return out;
}

inline bool edge_support_has_coverage(
    const Mesh2D &mesh,
    const Edge &edge,
    const std::vector<int> &node_ids,
    const std::vector<double> &nodal_length_scale,
    int required_support)
{
    if ((int)node_ids.size() < required_support)
        return false;

    const EdgeSupportCoverage coverage =
        evaluate_edge_support_coverage(mesh, edge, node_ids, nodal_length_scale);
    if (coverage.tangent_neg <= 0 || coverage.tangent_pos <= 0)
        return false;

    if (edge.is_boundary)
        return coverage.boundary_interior >= 2;

    return coverage.normal_neg > 0 && coverage.normal_pos > 0;
}

inline bool solve_dense_linear_system(
    const std::vector<double> &A_in,
    const std::vector<double> &b_in,
    std::vector<double> &x_out,
    double &cond_est)
{
    const int n = static_cast<int>(b_in.size());
    cond_est = std::numeric_limits<double>::infinity();
    if (n <= 0 || (int)A_in.size() != n * n)
        return false;

    std::vector<double> A = A_in;
    std::vector<double> x = b_in;
    double max_pivot = 0.0;
    double min_pivot = std::numeric_limits<double>::infinity();

    for (int col = 0; col < n; ++col)
    {
        int piv = col;
        double best = std::abs(A[(size_t)col * n + col]);
        for (int row = col + 1; row < n; ++row)
        {
            const double cand = std::abs(A[(size_t)row * n + col]);
            if (cand > best)
            {
                best = cand;
                piv = row;
            }
        }
        if (!(best > 1.0e-14))
            return false;
        if (piv != col)
        {
            for (int c = col; c < n; ++c)
                std::swap(A[(size_t)col * n + c], A[(size_t)piv * n + c]);
            std::swap(x[(size_t)col], x[(size_t)piv]);
        }

        const double pivot = A[(size_t)col * n + col];
        max_pivot = std::max(max_pivot, std::abs(pivot));
        min_pivot = std::min(min_pivot, std::abs(pivot));

        for (int row = col + 1; row < n; ++row)
        {
            const double factor = A[(size_t)row * n + col] / pivot;
            A[(size_t)row * n + col] = 0.0;
            if (std::abs(factor) <= 1.0e-18)
                continue;
            for (int c = col + 1; c < n; ++c)
                A[(size_t)row * n + c] -= factor * A[(size_t)col * n + c];
            x[(size_t)row] -= factor * x[(size_t)col];
        }
    }

    x_out.assign((size_t)n, 0.0);
    for (int row = n - 1; row >= 0; --row)
    {
        double rhs = x[(size_t)row];
        for (int c = row + 1; c < n; ++c)
            rhs -= A[(size_t)row * n + c] * x_out[(size_t)c];
        const double pivot = A[(size_t)row * n + row];
        if (!(std::abs(pivot) > 1.0e-14))
            return false;
        x_out[(size_t)row] = rhs / pivot;
    }

    cond_est = max_pivot / std::max(1.0e-30, min_pivot);
    return std::isfinite(cond_est);
}

inline WeightedSupport collect_support(
    const EFGMIContext2D &ctx,
    double x,
    double y,
    double growth,
    bool allow_residual_fallback,
    const std::vector<int> *candidate_node_ids = nullptr,
    int required_min_support = -1,
    int reference_tri_id = -1)
{
    WeightedSupport out;
    const Mesh2D &mesh = *ctx.mesh;
    const int min_support =
        (required_min_support > 0) ? required_min_support : ctx.min_support;
    const double support_scale = support_scale_for_triangle(ctx, reference_tri_id);
    out.node_ids.reserve(mesh.nodes.size());
    out.w.reserve(mesh.nodes.size());
    out.dwdx.reserve(mesh.nodes.size());
    out.dwdy.reserve(mesh.nodes.size());

    const auto consider_node = [&](int node_id)
    {
        if (!node_belongs_to_material_component(ctx, node_id, reference_tri_id))
            return;
        const Node2D &node = mesh.nodes[(size_t)node_id];
        const double dx = x - node.x;
        const double dy = y - node.y;
        const double h = std::max(1.0e-12, growth * support_scale * ctx.nodal_length_scale[(size_t)node_id]);

        double w = 0.0;
        double dwdx = 0.0;
        double dwdy = 0.0;
        compact_singular_weight(dx, dy, h, ctx.weight_epsilon, w, dwdx, dwdy);
        if (w <= 0.0)
            return;

        out.node_ids.push_back(node_id);
        out.w.push_back(w);
        out.dwdx.push_back(dwdx);
        out.dwdy.push_back(dwdy);
    };

    if (candidate_node_ids != nullptr)
    {
        for (const int node_id : *candidate_node_ids)
            consider_node(node_id);
    }
    else
    {
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
            consider_node(static_cast<int>(node_id));
    }
    out.compact_count = static_cast<int>(out.node_ids.size());

    if ((int)out.node_ids.size() >= min_support || !allow_residual_fallback)
        return out;

    std::vector<std::pair<double, int>> nearest;
    nearest.reserve((candidate_node_ids != nullptr) ? candidate_node_ids->size() : mesh.nodes.size());
    const auto append_nearest = [&](int node_id)
    {
        if (!node_belongs_to_material_component(ctx, node_id, reference_tri_id))
            return;
        const Node2D &node = mesh.nodes[(size_t)node_id];
        const double dx = x - node.x;
        const double dy = y - node.y;
        nearest.push_back({dx * dx + dy * dy, node_id});
    };
    if (candidate_node_ids != nullptr)
    {
        for (const int node_id : *candidate_node_ids)
            append_nearest(node_id);
    }
    else
    {
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
            append_nearest(static_cast<int>(node_id));
    }
    std::sort(nearest.begin(), nearest.end(), [](const auto &lhs, const auto &rhs)
              { return lhs.first < rhs.first; });

    for (const auto &[dist2, node_id] : nearest)
    {
        (void)dist2;
        if (std::find(out.node_ids.begin(), out.node_ids.end(), node_id) != out.node_ids.end())
            continue;

        const Node2D &node = mesh.nodes[(size_t)node_id];
        const double dx = x - node.x;
        const double dy = y - node.y;
        const double h = std::max(1.0e-12, growth * support_scale * ctx.nodal_length_scale[(size_t)node_id]);
        double w = 0.0;
        double dwdx = 0.0;
        double dwdy = 0.0;
        compact_singular_weight(dx, dy, h, ctx.weight_epsilon, w, dwdx, dwdy);
        if (w <= 0.0)
        {
            // So usa peso residual como ultimo recurso, apos esgotar a
            // expansao do dominio de influencia real.
            w = 1.0 / (1.0e-12 + dist2);
            out.used_residual_fallback = true;
        }

        out.node_ids.push_back(node_id);
        out.w.push_back(w);
        out.dwdx.push_back(dwdx);
        out.dwdy.push_back(dwdy);
        if ((int)out.node_ids.size() >= min_support)
            break;
    }

    return out;
}

inline bool build_moment_system(
    const Mesh2D &mesh,
    const WeightedSupport &support,
    Mat3 &A,
    Mat3 &dA_dx,
    Mat3 &dA_dy,
    Mat3 &A_inv,
    double &cond_est)
{
    A = {};
    dA_dx = {};
    dA_dy = {};
    for (size_t k = 0; k < support.node_ids.size(); ++k)
    {
        const Node2D &node = mesh.nodes[(size_t)support.node_ids[k]];
        const std::array<double, 3> p = {1.0, node.x, node.y};
        add_outer_product_scaled(A, p, support.w[k]);
        add_outer_product_scaled(dA_dx, p, support.dwdx[k]);
        add_outer_product_scaled(dA_dy, p, support.dwdy[k]);
    }
    return inverse3(A, A_inv, cond_est) && cond_est < 1.0e12;
}

} // namespace detail

inline EFGMIContext2D make_context(
    const Mesh2D &mesh,
    double support_scale = 2.5,
    double weight_epsilon = 1.0e-3,
    int min_support = 6)
{
    EFGMIContext2D ctx;
    ctx.mesh = &mesh;
    ctx.nodal_length_scale = detail::compute_nodal_length_scale(mesh);
    ctx.support_scale =
        detail::parse_env_positive_double("TP3485_EFGMI_NODAL_SUPPORT_SCALE", support_scale);
    ctx.interface_support_scale =
        detail::parse_env_positive_double("TP3485_EFGMI_NODAL_INTERFACE_SUPPORT_SCALE",
                                          std::max(ctx.support_scale, 4.0));
    ctx.weight_epsilon = weight_epsilon;
    ctx.min_support =
        detail::parse_env_positive_int("TP3485_EFGMI_NODAL_MIN_SUPPORT", min_support);
    ctx.interface_min_support =
        detail::parse_env_positive_int("TP3485_EFGMI_NODAL_INTERFACE_MIN_SUPPORT",
                                       std::max(ctx.min_support, 9));
    return ctx;
}

inline EFGMIContext2D make_context(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    double support_scale = 2.5,
    double weight_epsilon = 1.0e-3,
    int min_support = 6)
{
    EFGMIContext2D ctx = make_context(mesh, support_scale, weight_epsilon, min_support);
    detail::build_material_component_metadata(ctx, eps_r_tri, mu_r_tri);
    return ctx;
}

inline EFGMIEdgeContext2D make_edge_context(
    const Mesh2D &mesh,
    const EdgeDofs &edge_dofs,
    double support_scale = 3.0,
    double weight_epsilon = 1.0e-3,
    int min_support = 9)
{
    EFGMIEdgeContext2D ctx;
    const int interior_min_support =
        detail::parse_env_positive_int("TP3485_EFGMI_EDGE_MIN_SUPPORT_INTERIOR",
                                       std::max(min_support, 12));
    const int boundary_min_support =
        detail::parse_env_positive_int("TP3485_EFGMI_EDGE_MIN_SUPPORT_BOUNDARY",
                                       std::max(interior_min_support, 15));
    const int near_boundary_bonus =
        detail::parse_env_positive_int("TP3485_EFGMI_EDGE_NEAR_BOUNDARY_BONUS", 2);
    const int patch_extra_edges =
        detail::parse_env_positive_int("TP3485_EFGMI_EDGE_PATCH_EXTRA_EDGES", 6);
    const int max_required_support =
        std::max(boundary_min_support, interior_min_support + near_boundary_bonus);

    ctx.nodal = make_context(mesh, support_scale, weight_epsilon, max_required_support);
    ctx.edge_dofs = &edge_dofs;
    ctx.node_to_edges.assign(mesh.nodes.size(), {});
    ctx.edge_support_nodes.assign(edge_dofs.edges.size(), {});
    ctx.edge_min_support.assign(edge_dofs.edges.size(), interior_min_support);
    ctx.edge_patch_nodes.assign(edge_dofs.edges.size(), {});
    ctx.edge_patch_edges.assign(edge_dofs.edges.size(), {});
    ctx.edge_patch_coeffs.assign(edge_dofs.edges.size(), {});
    ctx.edge_patch_cond_est.assign(edge_dofs.edges.size(), std::numeric_limits<double>::infinity());
    ctx.edge_patch_ok.assign(edge_dofs.edges.size(), 0);
    for (size_t edge_id = 0; edge_id < edge_dofs.edges.size(); ++edge_id)
    {
        const Edge &edge = edge_dofs.edges[edge_id];
        ctx.node_to_edges[(size_t)edge.n0].push_back(static_cast<int>(edge_id));
        ctx.node_to_edges[(size_t)edge.n1].push_back(static_cast<int>(edge_id));
    }

    std::vector<std::vector<int>> vertex_support_nodes(mesh.nodes.size());
    for (size_t center_id = 0; center_id < mesh.nodes.size(); ++center_id)
    {
        const Node2D &center = mesh.nodes[center_id];
        const double radius =
            std::max(1.0e-12, ctx.nodal.support_scale * ctx.nodal.nodal_length_scale[center_id]);
        const double radius2 = radius * radius;
        auto &support_nodes = vertex_support_nodes[center_id];
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            const Node2D &node = mesh.nodes[node_id];
            const double dx = node.x - center.x;
            const double dy = node.y - center.y;
            if (dx * dx + dy * dy <= radius2)
                support_nodes.push_back(static_cast<int>(node_id));
        }
        if (support_nodes.empty())
            support_nodes.push_back(static_cast<int>(center_id));
    }

    std::vector<unsigned char> mark(mesh.nodes.size(), 0);
    std::vector<int> touched;
    for (size_t edge_id = 0; edge_id < edge_dofs.edges.size(); ++edge_id)
    {
        const Edge &edge = edge_dofs.edges[edge_id];
        const Node2D &n0 = mesh.nodes[(size_t)edge.n0];
        const Node2D &n1 = mesh.nodes[(size_t)edge.n1];
        const bool near_boundary =
            edge.is_boundary || n0.is_boundary || n1.is_boundary;
        const int required_support = edge.is_boundary
                                         ? boundary_min_support
                                         : (near_boundary
                                                ? std::max(interior_min_support,
                                                           interior_min_support + near_boundary_bonus)
                                                : interior_min_support);
        const double mx = 0.5 * (n0.x + n1.x);
        const double my = 0.5 * (n0.y + n1.y);
        ctx.edge_min_support[edge_id] = required_support;

        auto &edge_nodes = ctx.edge_support_nodes[edge_id];
        edge_nodes.reserve(
            vertex_support_nodes[(size_t)edge.n0].size() +
            vertex_support_nodes[(size_t)edge.n1].size());
        touched.clear();

        const auto add_candidate = [&](int node_id)
        {
            if (mark[(size_t)node_id] != 0)
                return;
            mark[(size_t)node_id] = 1;
            touched.push_back(node_id);
            edge_nodes.push_back(node_id);
        };

        add_candidate(edge.n0);
        add_candidate(edge.n1);
        for (const int node_id : vertex_support_nodes[(size_t)edge.n0])
            add_candidate(node_id);
        for (const int node_id : vertex_support_nodes[(size_t)edge.n1])
            add_candidate(node_id);

        std::vector<std::pair<double, int>> nearest;
        nearest.reserve(mesh.nodes.size());
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            const Node2D &node = mesh.nodes[node_id];
            const double dx = node.x - mx;
            const double dy = node.y - my;
            nearest.push_back({dx * dx + dy * dy, static_cast<int>(node_id)});
        }
        std::sort(nearest.begin(), nearest.end(), [](const auto &lhs, const auto &rhs)
                  { return lhs.first < rhs.first; });
        for (const auto &[dist2, node_id] : nearest)
        {
            (void)dist2;
            if (detail::edge_support_has_coverage(
                    mesh,
                    edge,
                    edge_nodes,
                    ctx.nodal.nodal_length_scale,
                    required_support))
                break;
            add_candidate(node_id);
        }

        for (const int node_id : touched)
            mark[(size_t)node_id] = 0;
    }

    std::vector<unsigned char> edge_mark(edge_dofs.edges.size(), 0);
    std::vector<unsigned char> patch_node_mark(mesh.nodes.size(), 0);
    std::vector<int> edge_touched;
    std::vector<int> node_touched;
    for (size_t edge_id = 0; edge_id < edge_dofs.edges.size(); ++edge_id)
    {
        const Edge &center_edge = edge_dofs.edges[edge_id];
        const Node2D &n0 = mesh.nodes[(size_t)center_edge.n0];
        const Node2D &n1 = mesh.nodes[(size_t)center_edge.n1];
        const double mx = 0.5 * (n0.x + n1.x);
        const double my = 0.5 * (n0.y + n1.y);
        const int max_patch_edges = std::max(3, ctx.edge_min_support[edge_id] + patch_extra_edges);

        auto &patch_nodes = ctx.edge_patch_nodes[edge_id];
        auto &patch_edges = ctx.edge_patch_edges[edge_id];
        patch_nodes = ctx.edge_support_nodes[edge_id];
        node_touched.clear();
        for (const int node_id : patch_nodes)
        {
            if (patch_node_mark[(size_t)node_id] != 0)
                continue;
            patch_node_mark[(size_t)node_id] = 1;
            node_touched.push_back(node_id);
        }

        std::vector<std::pair<double, int>> edge_candidates;
        edge_candidates.reserve(ctx.edge_support_nodes[edge_id].size() * 6);
        edge_touched.clear();
        const auto add_patch_edge = [&](int candidate_edge_id)
        {
            if (edge_mark[(size_t)candidate_edge_id] != 0)
                return;
            edge_mark[(size_t)candidate_edge_id] = 1;
            edge_touched.push_back(candidate_edge_id);

            const Edge &candidate_edge = edge_dofs.edges[(size_t)candidate_edge_id];
            const Node2D &ea = mesh.nodes[(size_t)candidate_edge.n0];
            const Node2D &eb = mesh.nodes[(size_t)candidate_edge.n1];
            const double emx = 0.5 * (ea.x + eb.x);
            const double emy = 0.5 * (ea.y + eb.y);
            const double dist2 = (emx - mx) * (emx - mx) + (emy - my) * (emy - my);
            edge_candidates.push_back({dist2, candidate_edge_id});
        };

        add_patch_edge(static_cast<int>(edge_id));
        for (const int node_id : ctx.edge_support_nodes[edge_id])
            for (const int candidate_edge_id : ctx.node_to_edges[(size_t)node_id])
                add_patch_edge(candidate_edge_id);

        std::sort(edge_candidates.begin(), edge_candidates.end(), [&](const auto &lhs, const auto &rhs)
                  {
                      if (lhs.second == (int)edge_id)
                          return true;
                      if (rhs.second == (int)edge_id)
                          return false;
                      return lhs.first < rhs.first;
                  });
        for (const auto &[dist2, candidate_edge_id] : edge_candidates)
        {
            (void)dist2;
            patch_edges.push_back(candidate_edge_id);
            if ((int)patch_edges.size() >= max_patch_edges)
                break;
        }

        for (const int candidate_edge_id : patch_edges)
        {
            const Edge &candidate_edge = edge_dofs.edges[(size_t)candidate_edge_id];
            if (patch_node_mark[(size_t)candidate_edge.n0] == 0)
            {
                patch_node_mark[(size_t)candidate_edge.n0] = 1;
                node_touched.push_back(candidate_edge.n0);
                patch_nodes.push_back(candidate_edge.n0);
            }
            if (patch_node_mark[(size_t)candidate_edge.n1] == 0)
            {
                patch_node_mark[(size_t)candidate_edge.n1] = 1;
                node_touched.push_back(candidate_edge.n1);
                patch_nodes.push_back(candidate_edge.n1);
            }
        }

        for (const int edge_marked : edge_touched)
            edge_mark[(size_t)edge_marked] = 0;
        for (const int node_marked : node_touched)
            patch_node_mark[(size_t)node_marked] = 0;
    }

    for (size_t edge_id = 0; edge_id < edge_dofs.edges.size(); ++edge_id)
    {
        double cond_est = std::numeric_limits<double>::infinity();
        std::vector<double> coeffs;
        const bool ok = build_edge_patch_correction(ctx, static_cast<int>(edge_id), coeffs, cond_est);
        ctx.edge_patch_coeffs[edge_id] = std::move(coeffs);
        ctx.edge_patch_cond_est[edge_id] = cond_est;
        ctx.edge_patch_ok[edge_id] = ok ? 1u : 0u;
    }

    return ctx;
}

inline ShapeSample evaluate_shape_impl(
    const EFGMIContext2D &ctx,
    const std::vector<int> *candidate_node_ids,
    double x,
    double y,
    bool snap_to_node = false,
    int required_min_support = -1,
    int reference_tri_id = -1)
{
    if (ctx.mesh == nullptr)
        throw std::runtime_error("EFGMIContext2D sem mesh associado.");

    const Mesh2D &mesh = *ctx.mesh;
    if (snap_to_node)
    {
        const auto make_snap_sample = [&](int node_id) -> ShapeSample
        {
            ShapeSample sample;
            sample.node_ids = {node_id};
            sample.phi = {1.0};
            sample.dphidx = {0.0};
            sample.dphidy = {0.0};
            return sample;
        };

        if (candidate_node_ids != nullptr)
        {
            for (const int node_id : *candidate_node_ids)
            {
                if (!detail::node_belongs_to_material_component(ctx, node_id, reference_tri_id))
                    continue;
                const Node2D &node = mesh.nodes[(size_t)node_id];
                const double dx = x - node.x;
                const double dy = y - node.y;
                const double h = std::max(1.0e-12, ctx.nodal_length_scale[(size_t)node_id]);
                if (std::sqrt(dx * dx + dy * dy) <= ctx.snap_tolerance * std::max(1.0, h))
                    return make_snap_sample(node_id);
            }
        }
        else
        {
            for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
            {
                if (!detail::node_belongs_to_material_component(
                        ctx,
                        static_cast<int>(node_id),
                        reference_tri_id))
                    continue;
                const Node2D &node = mesh.nodes[node_id];
                const double dx = x - node.x;
                const double dy = y - node.y;
                const double h = std::max(1.0e-12, ctx.nodal_length_scale[node_id]);
                if (std::sqrt(dx * dx + dy * dy) <= ctx.snap_tolerance * std::max(1.0, h))
                    return make_snap_sample(static_cast<int>(node_id));
            }
        }
    }

    detail::WeightedSupport support;
    Mat3 A{}, dA_dx{}, dA_dy{}, A_inv{};
    double cond_est = std::numeric_limits<double>::infinity();
    bool ok = false;
    const int min_support =
        detail::required_min_support_for_triangle(ctx, reference_tri_id, required_min_support);
    double growth = 1.0;
    for (int step = 0; step < ctx.max_expand_steps; ++step)
    {
        const bool allow_residual_fallback = (step + 1 == ctx.max_expand_steps);
        support = detail::collect_support(
            ctx,
            x,
            y,
            growth,
            allow_residual_fallback,
            candidate_node_ids,
            min_support,
            reference_tri_id);
        if (support.compact_count < min_support && !allow_residual_fallback)
        {
            growth *= 1.5;
            continue;
        }
        ok = detail::build_moment_system(mesh, support, A, dA_dx, dA_dy, A_inv, cond_est);
        if (ok)
            break;
        growth *= 1.5;
    }
    if (!ok)
        throw std::runtime_error("EFGMI: matriz de momentos singular/instavel no ponto de integracao.");

    const std::array<double, 3> p = {1.0, x, y};
    const std::array<double, 3> dp_dx = {0.0, 1.0, 0.0};
    const std::array<double, 3> dp_dy = {0.0, 0.0, 1.0};
    const std::array<double, 3> g = detail::mat_vec(A_inv, p);

    ShapeSample sample;
    sample.node_ids = support.node_ids;
    sample.phi.resize(support.node_ids.size(), 0.0);
    sample.dphidx.resize(support.node_ids.size(), 0.0);
    sample.dphidy.resize(support.node_ids.size(), 0.0);

    for (size_t k = 0; k < support.node_ids.size(); ++k)
    {
        const Node2D &node = mesh.nodes[(size_t)support.node_ids[k]];
        const std::array<double, 3> pi = {1.0, node.x, node.y};
        const std::array<double, 3> Bi = {
            support.w[k] * pi[0],
            support.w[k] * pi[1],
            support.w[k] * pi[2]};
        const std::array<double, 3> dBi_dx = {
            support.dwdx[k] * pi[0],
            support.dwdx[k] * pi[1],
            support.dwdx[k] * pi[2]};
        const std::array<double, 3> dBi_dy = {
            support.dwdy[k] * pi[0],
            support.dwdy[k] * pi[1],
            support.dwdy[k] * pi[2]};

        const std::array<double, 3> v = detail::mat_vec(A_inv, Bi);
        sample.phi[k] = detail::dot3(p, v);

        const Mat3 dAinv_dx = [](
                                   const Mat3 &A_inv_local,
                                   const Mat3 &dA_local) {
            Mat3 tmp = detail::mat_mul(A_inv_local, dA_local);
            tmp = detail::mat_mul(tmp, A_inv_local);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    tmp.a[i][j] = -tmp.a[i][j];
            return tmp;
        }(A_inv, dA_dx);

        const Mat3 dAinv_dy = [](
                                   const Mat3 &A_inv_local,
                                   const Mat3 &dA_local) {
            Mat3 tmp = detail::mat_mul(A_inv_local, dA_local);
            tmp = detail::mat_mul(tmp, A_inv_local);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    tmp.a[i][j] = -tmp.a[i][j];
            return tmp;
        }(A_inv, dA_dy);

        const std::array<double, 3> dAinvBi_dx = detail::mat_vec(dAinv_dx, Bi);
        const std::array<double, 3> dAinvBi_dy = detail::mat_vec(dAinv_dy, Bi);
        const std::array<double, 3> AinvdBi_dx = detail::mat_vec(A_inv, dBi_dx);
        const std::array<double, 3> AinvdBi_dy = detail::mat_vec(A_inv, dBi_dy);

        sample.dphidx[k] =
            detail::dot3(dp_dx, v) +
            detail::dot3(p, dAinvBi_dx) +
            detail::dot3(p, AinvdBi_dx);
        sample.dphidy[k] =
            detail::dot3(dp_dy, v) +
            detail::dot3(p, dAinvBi_dy) +
            detail::dot3(p, AinvdBi_dy);
    }

    return sample;
}

inline ShapeSample evaluate_shape(
    const EFGMIContext2D &ctx,
    double x,
    double y,
    bool snap_to_node,
    int required_min_support,
    int reference_tri_id)
{
    return evaluate_shape_impl(ctx, nullptr, x, y, snap_to_node, required_min_support, reference_tri_id);
}

inline ShapeSample evaluate_shape_on_candidates(
    const EFGMIContext2D &ctx,
    const std::vector<int> &candidate_node_ids,
    double x,
    double y,
    bool snap_to_node,
    int required_min_support,
    int reference_tri_id)
{
    return evaluate_shape_impl(
        ctx,
        &candidate_node_ids,
        x,
        y,
        snap_to_node,
        required_min_support,
        reference_tri_id);
}

inline ScalarFieldSample evaluate_scalar_field(
    const EFGMIContext2D &ctx,
    const std::vector<double> &nodal_coeffs,
    double x,
    double y,
    bool snap_to_node = false)
{
    if (ctx.mesh == nullptr)
        throw std::runtime_error("EFGMIContext2D sem mesh associado.");
    if (nodal_coeffs.size() != ctx.mesh->nodes.size())
        throw std::runtime_error("EFGMI: nodal_coeffs.size() deve coincidir com mesh.nodes.size().");

    const ShapeSample sample = evaluate_shape(ctx, x, y, snap_to_node);
    ScalarFieldSample out;
    for (size_t k = 0; k < sample.node_ids.size(); ++k)
    {
        const double coeff = nodal_coeffs[(size_t)sample.node_ids[k]];
        out.value += coeff * sample.phi[k];
        out.dx += coeff * sample.dphidx[k];
        out.dy += coeff * sample.dphidy[k];
    }
    return out;
}

inline GlobalEdgeBasisValue evaluate_raw_global_edge_basis_from_sample(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    const ShapeSample &sample)
{
    if (ctx.edge_dofs == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D sem EdgeDofs associado.");
    if (edge_id < 0 || edge_id >= (int)ctx.edge_dofs->edges.size())
        throw std::runtime_error("EFGMI: edge_id invalido.");

    const Edge &edge = ctx.edge_dofs->edges[(size_t)edge_id];
    int idx_i = -1;
    int idx_j = -1;
    for (size_t k = 0; k < sample.node_ids.size(); ++k)
    {
        if (sample.node_ids[k] == edge.n0)
            idx_i = static_cast<int>(k);
        else if (sample.node_ids[k] == edge.n1)
            idx_j = static_cast<int>(k);
    }

    GlobalEdgeBasisValue out;
    if (idx_i < 0 || idx_j < 0)
        return out;

    const Node2D &ni = ctx.nodal.mesh->nodes[(size_t)edge.n0];
    const Node2D &nj = ctx.nodal.mesh->nodes[(size_t)edge.n1];
    const double dx = nj.x - ni.x;
    const double dy = nj.y - ni.y;
    const double L = std::sqrt(dx * dx + dy * dy);
    if (!(L > 0.0))
        return out;

    out.vx = L * (sample.phi[(size_t)idx_i] * sample.dphidx[(size_t)idx_j] -
                  sample.phi[(size_t)idx_j] * sample.dphidx[(size_t)idx_i]);
    out.vy = L * (sample.phi[(size_t)idx_i] * sample.dphidy[(size_t)idx_j] -
                  sample.phi[(size_t)idx_j] * sample.dphidy[(size_t)idx_i]);
    out.curl = 2.0 * L *
               (sample.dphidx[(size_t)idx_i] * sample.dphidy[(size_t)idx_j] -
                sample.dphidy[(size_t)idx_i] * sample.dphidx[(size_t)idx_j]);
    return out;
}

inline bool build_edge_patch_correction(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    std::vector<double> &coeffs,
    double &cond_est)
{
    coeffs.clear();
    cond_est = std::numeric_limits<double>::infinity();
    if (ctx.edge_dofs == nullptr || ctx.nodal.mesh == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D incompleto para correcao por patch.");
    if (edge_id < 0 || edge_id >= (int)ctx.edge_dofs->edges.size())
        throw std::runtime_error("EFGMI: edge_id invalido na correcao por patch.");

    const auto &patch_nodes = ctx.edge_patch_nodes[(size_t)edge_id];
    const auto &patch_edges = ctx.edge_patch_edges[(size_t)edge_id];
    const int n = static_cast<int>(patch_edges.size());
    if (patch_nodes.empty() || n <= 0)
        return false;

    int center_index = -1;
    for (int i = 0; i < n; ++i)
    {
        if (patch_edges[(size_t)i] == edge_id)
        {
            center_index = i;
            break;
        }
    }
    if (center_index < 0)
        return false;

    const Mesh2D &mesh = *ctx.nodal.mesh;
    std::vector<double> M((size_t)n * n, 0.0);
    for (int row = 0; row < n; ++row)
    {
        const Edge &test_edge = ctx.edge_dofs->edges[(size_t)patch_edges[(size_t)row]];
        const Node2D &na = mesh.nodes[(size_t)test_edge.n0];
        const Node2D &nb = mesh.nodes[(size_t)test_edge.n1];
        const double dx = nb.x - na.x;
        const double dy = nb.y - na.y;

        for (const auto &qp : kLineQuadP4)
        {
            const double x = (1.0 - qp.s) * na.x + qp.s * nb.x;
            const double y = (1.0 - qp.s) * na.y + qp.s * nb.y;
            const ShapeSample sample = evaluate_shape_on_candidates(
                ctx.nodal,
                patch_nodes,
                x,
                y,
                false,
                ctx.edge_min_support[(size_t)edge_id]);
            for (int col = 0; col < n; ++col)
            {
                const auto raw_basis =
                    evaluate_raw_global_edge_basis_from_sample(ctx, patch_edges[(size_t)col], sample);
                M[(size_t)row * n + col] += qp.weight * (raw_basis.vx * dx + raw_basis.vy * dy);
            }
        }
    }

    std::vector<double> rhs((size_t)n, 0.0);
    const Edge &center_edge = ctx.edge_dofs->edges[(size_t)edge_id];
    const Node2D &n0 = mesh.nodes[(size_t)center_edge.n0];
    const Node2D &n1 = mesh.nodes[(size_t)center_edge.n1];
    const double lx = n1.x - n0.x;
    const double ly = n1.y - n0.y;
    rhs[(size_t)center_index] = std::sqrt(lx * lx + ly * ly);

    if (!detail::solve_dense_linear_system(M, rhs, coeffs, cond_est))
        return false;
    return cond_est < 1.0e8;
}

inline GlobalEdgeBasisValue evaluate_global_edge_basis(
    const EFGMIEdgeContext2D &ctx,
    int edge_id,
    double x,
    double y)
{
    if (ctx.edge_dofs == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D sem EdgeDofs associado.");
    if (edge_id < 0 || edge_id >= (int)ctx.edge_dofs->edges.size())
        throw std::runtime_error("EFGMI: edge_id invalido.");

    const bool patch_ok =
        ((size_t)edge_id < ctx.edge_patch_ok.size()) && (ctx.edge_patch_ok[(size_t)edge_id] != 0);
    const auto &candidate_nodes =
        (patch_ok && (size_t)edge_id < ctx.edge_patch_nodes.size() &&
         !ctx.edge_patch_nodes[(size_t)edge_id].empty())
            ? ctx.edge_patch_nodes[(size_t)edge_id]
            : ctx.edge_support_nodes[(size_t)edge_id];
    const int required_min_support =
        ((size_t)edge_id < ctx.edge_min_support.size()) ? ctx.edge_min_support[(size_t)edge_id]
                                                        : ctx.nodal.min_support;
    const ShapeSample sample =
        evaluate_shape_on_candidates(ctx.nodal, candidate_nodes, x, y, false, required_min_support);

    if (!patch_ok ||
        (size_t)edge_id >= ctx.edge_patch_edges.size() ||
        (size_t)edge_id >= ctx.edge_patch_coeffs.size() ||
        ctx.edge_patch_edges[(size_t)edge_id].size() != ctx.edge_patch_coeffs[(size_t)edge_id].size())
    {
        return evaluate_raw_global_edge_basis_from_sample(ctx, edge_id, sample);
    }

    GlobalEdgeBasisValue out;
    const auto &patch_edges = ctx.edge_patch_edges[(size_t)edge_id];
    const auto &coeffs = ctx.edge_patch_coeffs[(size_t)edge_id];
    for (size_t k = 0; k < patch_edges.size(); ++k)
    {
        const auto raw_basis = evaluate_raw_global_edge_basis_from_sample(ctx, patch_edges[k], sample);
        out.vx += coeffs[k] * raw_basis.vx;
        out.vy += coeffs[k] * raw_basis.vy;
        out.curl += coeffs[k] * raw_basis.curl;
    }
    return out;
}

inline bool build_triangle_edge_transform(
    const EFGMIEdgeContext2D &ctx,
    int tri_id,
    Mat3 &transform,
    double &cond_est)
{
    transform = detail::identity3();
    cond_est = std::numeric_limits<double>::infinity();
    if (ctx.edge_dofs == nullptr || ctx.nodal.mesh == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D incompleto para transformacao de aresta.");
    if (tri_id < 0 || tri_id >= (int)ctx.nodal.mesh->tris.size())
        throw std::runtime_error("EFGMI: tri_id invalido na transformacao de aresta.");

    const Mesh2D &mesh = *ctx.nodal.mesh;
    const Tri &tri = mesh.tris[(size_t)tri_id];
    const TriEdges &te = ctx.edge_dofs->tri_edges[(size_t)tri_id];
    const TriGeomEdge tg = tri_geom_edge(mesh, tri);
    static constexpr int edge_i[3] = {0, 1, 2};
    static constexpr int edge_j[3] = {1, 2, 0};

    Mat3 H{};
    for (int row = 0; row < 3; ++row)
    {
        const Node2D &na = mesh.nodes[(size_t)tri.v[edge_i[row]]];
        const Node2D &nb = mesh.nodes[(size_t)tri.v[edge_j[row]]];
        const double dx = nb.x - na.x;
        const double dy = nb.y - na.y;

        for (const auto &qp : kLineQuadP4)
        {
            const double x = (1.0 - qp.s) * na.x + qp.s * nb.x;
            const double y = (1.0 - qp.s) * na.y + qp.s * nb.y;
            for (int col = 0; col < 3; ++col)
            {
                const auto raw_basis = evaluate_global_edge_basis(ctx, te.e[col], x, y);
                const double sign = static_cast<double>(te.sgn[col]);
                H.a[row][col] += qp.weight *
                                 ((sign * raw_basis.vx) * dx + (sign * raw_basis.vy) * dy);
            }
        }
    }

    Mat3 H_inv{};
    if (!detail::inverse3(H, H_inv, cond_est) || cond_est >= detail::kConsistencyCondLimit)
        return false;

    Mat3 target{};
    for (int i = 0; i < 3; ++i)
        target.a[i][i] = tg.L[(size_t)i];
    transform = detail::mat_mul(H_inv, target);
    return true;
}

inline LocalTriangleEdgeBasisSample evaluate_triangle_edge_basis(
    const EFGMIEdgeContext2D &ctx,
    int tri_id,
    double x,
    double y)
{
    if (ctx.edge_dofs == nullptr || ctx.nodal.mesh == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D sem dados suficientes para base vetorial local.");
    if (tri_id < 0 || tri_id >= (int)ctx.edge_dofs->tri_edges.size())
        throw std::runtime_error("EFGMI: tri_id invalido na avaliacao da base vetorial.");

    const TriEdges &te = ctx.edge_dofs->tri_edges[(size_t)tri_id];
    LocalTriangleEdgeBasisSample raw_local;
    for (int m = 0; m < 3; ++m)
    {
        const auto global_basis = evaluate_global_edge_basis(ctx, te.e[m], x, y);
        const double sign = static_cast<double>(te.sgn[m]);
        raw_local.vx[(size_t)m] = sign * global_basis.vx;
        raw_local.vy[(size_t)m] = sign * global_basis.vy;
        raw_local.curl[(size_t)m] = sign * global_basis.curl;
    }

    if ((size_t)tri_id >= ctx.tri_edge_transform.size())
        return raw_local;
    return apply_local_triangle_edge_transform(
        raw_local,
        ctx.tri_edge_transform[(size_t)tri_id]);
}

inline LocalTriangleShapeSample extract_local_triangle_shape(
    const ShapeSample &sample,
    const Tri &tri)
{
    LocalTriangleShapeSample out;
    for (int lv = 0; lv < 3; ++lv)
    {
        const int node_id = tri.v[lv];
        for (size_t k = 0; k < sample.node_ids.size(); ++k)
        {
            if (sample.node_ids[k] != node_id)
                continue;
            out.phi[(size_t)lv] = sample.phi[k];
            out.dphidx[(size_t)lv] = sample.dphidx[k];
            out.dphidy[(size_t)lv] = sample.dphidy[k];
            break;
        }
    }
    return out;
}

inline LocalTriangleShapeSample apply_local_triangle_transform(
    const LocalTriangleShapeSample &raw,
    const Mat3 &transform)
{
    LocalTriangleShapeSample out;
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            out.phi[(size_t)j] += raw.phi[(size_t)i] * transform.a[i][j];
            out.dphidx[(size_t)j] += raw.dphidx[(size_t)i] * transform.a[i][j];
            out.dphidy[(size_t)j] += raw.dphidy[(size_t)i] * transform.a[i][j];
        }
    }
    return out;
}

inline LocalTriangleShapeSample normalize_local_triangle_shape(
    const LocalTriangleShapeSample &raw)
{
    LocalTriangleShapeSample out = raw;
    const double phi_sum = raw.phi[0] + raw.phi[1] + raw.phi[2];
    if (std::abs(phi_sum) <= 1.0e-14)
        return out;

    const double dphi_sum_x = raw.dphidx[0] + raw.dphidx[1] + raw.dphidx[2];
    const double dphi_sum_y = raw.dphidy[0] + raw.dphidy[1] + raw.dphidy[2];
    for (int lv = 0; lv < 3; ++lv)
    {
        const double phi_raw = raw.phi[(size_t)lv];
        const double dphi_raw_x = raw.dphidx[(size_t)lv];
        const double dphi_raw_y = raw.dphidy[(size_t)lv];
        out.phi[(size_t)lv] = phi_raw / phi_sum;
        out.dphidx[(size_t)lv] =
            (dphi_raw_x * phi_sum - phi_raw * dphi_sum_x) / (phi_sum * phi_sum);
        out.dphidy[(size_t)lv] =
            (dphi_raw_y * phi_sum - phi_raw * dphi_sum_y) / (phi_sum * phi_sum);
    }
    return out;
}

inline LocalTriangleShapeSample make_consistent_local_triangle_shape(
    const ShapeSample &sample,
    const Tri &tri,
    const LocalTriangleConsistency &consistency)
{
    auto local = normalize_local_triangle_shape(extract_local_triangle_shape(sample, tri));
    if (consistency.nodal_ok)
        local = apply_local_triangle_transform(local, consistency.nodal_transform);
    return local;
}

inline LocalTriangleEdgeBasisSample build_local_triangle_edge_basis_raw(
    const TriGeomEdge &tg,
    const LocalTriangleShapeSample &local_shape)
{
    LocalTriangleEdgeBasisSample out;
    static constexpr int edge_i[3] = {0, 1, 2};
    static constexpr int edge_j[3] = {1, 2, 0};

    for (int m = 0; m < 3; ++m)
    {
        const int i = edge_i[m];
        const int j = edge_j[m];
        const double L = tg.L[(size_t)m];
        out.vx[(size_t)m] =
            L * (local_shape.phi[(size_t)i] * local_shape.dphidx[(size_t)j] -
                 local_shape.phi[(size_t)j] * local_shape.dphidx[(size_t)i]);
        out.vy[(size_t)m] =
            L * (local_shape.phi[(size_t)i] * local_shape.dphidy[(size_t)j] -
                 local_shape.phi[(size_t)j] * local_shape.dphidy[(size_t)i]);
        out.curl[(size_t)m] =
            2.0 * L * (local_shape.dphidx[(size_t)i] * local_shape.dphidy[(size_t)j] -
                       local_shape.dphidy[(size_t)i] * local_shape.dphidx[(size_t)j]);
    }
    return out;
}

inline LocalTriangleEdgeBasisSample apply_local_triangle_edge_transform(
    const LocalTriangleEdgeBasisSample &raw,
    const Mat3 &transform)
{
    LocalTriangleEdgeBasisSample out;
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            out.vx[(size_t)j] += raw.vx[(size_t)i] * transform.a[i][j];
            out.vy[(size_t)j] += raw.vy[(size_t)i] * transform.a[i][j];
            out.curl[(size_t)j] += raw.curl[(size_t)i] * transform.a[i][j];
        }
    }
    return out;
}

inline LocalTriangleEdgeBasisSample make_consistent_local_triangle_edge_basis(
    const ShapeSample &sample,
    const Tri &tri,
    const TriGeomEdge &tg,
    const LocalTriangleConsistency &consistency)
{
    const auto local_shape = make_consistent_local_triangle_shape(sample, tri, consistency);
    const auto raw_basis = build_local_triangle_edge_basis_raw(tg, local_shape);
    return apply_local_triangle_edge_transform(raw_basis, consistency.edge_transform);
}

inline bool build_local_triangle_nodal_transform(
    const EFGMIContext2D &ctx,
    const Tri &tri,
    Mat3 &transform,
    double &cond_est)
{
    if (ctx.mesh == nullptr)
        throw std::runtime_error("EFGMIContext2D sem mesh associado.");

    const Mesh2D &mesh = *ctx.mesh;
    Mat3 R{};
    for (int row = 0; row < 3; ++row)
    {
        const Node2D &node = mesh.nodes[(size_t)tri.v[row]];
        const auto sample = evaluate_shape(ctx, node.x, node.y, false);
        const auto local_shape = extract_local_triangle_shape(sample, tri);
        for (int col = 0; col < 3; ++col)
            R.a[row][col] = local_shape.phi[(size_t)col];
    }

    Mat3 R_inv{};
    if (!detail::inverse3(R, R_inv, cond_est) || cond_est >= detail::kConsistencyCondLimit)
    {
        transform = detail::identity3();
        return false;
    }

    transform = R_inv;
    return true;
}

inline bool build_local_triangle_edge_transform(
    const EFGMIContext2D &ctx,
    const Tri &tri,
    const TriGeomEdge &tg,
    const Mat3 &nodal_transform,
    Mat3 &transform,
    double &cond_est)
{
    if (ctx.mesh == nullptr)
        throw std::runtime_error("EFGMIContext2D sem mesh associado.");

    const Mesh2D &mesh = *ctx.mesh;
    Mat3 H{};
    static constexpr int edge_i[3] = {0, 1, 2};
    static constexpr int edge_j[3] = {1, 2, 0};

    for (int edge_id = 0; edge_id < 3; ++edge_id)
    {
        const Node2D &na = mesh.nodes[(size_t)tri.v[edge_i[edge_id]]];
        const Node2D &nb = mesh.nodes[(size_t)tri.v[edge_j[edge_id]]];
        const double dx = nb.x - na.x;
        const double dy = nb.y - na.y;

        for (const auto &qp : kLineQuadP4)
        {
            const double x = (1.0 - qp.s) * na.x + qp.s * nb.x;
            const double y = (1.0 - qp.s) * na.y + qp.s * nb.y;
            const auto sample = evaluate_shape(ctx, x, y, false);
            auto local_shape = normalize_local_triangle_shape(
                extract_local_triangle_shape(sample, tri));
            local_shape = apply_local_triangle_transform(local_shape, nodal_transform);
            const auto raw_basis = build_local_triangle_edge_basis_raw(tg, local_shape);
            for (int m = 0; m < 3; ++m)
            {
                H.a[edge_id][m] += qp.weight *
                                   (raw_basis.vx[(size_t)m] * dx +
                                    raw_basis.vy[(size_t)m] * dy);
            }
        }
    }

    Mat3 H_inv{};
    if (!detail::inverse3(H, H_inv, cond_est) || cond_est >= detail::kConsistencyCondLimit)
    {
        transform = detail::identity3();
        return false;
    }

    Mat3 target = {};
    for (int i = 0; i < 3; ++i)
        target.a[i][i] = tg.L[(size_t)i];
    transform = detail::mat_mul(H_inv, target);
    return true;
}

inline LocalTriangleConsistency build_local_triangle_consistency(
    const EFGMIContext2D &ctx,
    const Tri &tri,
    const TriGeomEdge &tg)
{
    LocalTriangleConsistency out;
    out.nodal_transform = detail::identity3();
    out.edge_transform = detail::identity3();
    out.nodal_ok = false;
    out.edge_ok = build_local_triangle_edge_transform(
        ctx,
        tri,
        tg,
        out.nodal_transform,
        out.edge_transform,
        out.edge_cond_est);
    return out;
}

inline EdgeBasisSample evaluate_edge_basis(
    const EFGMIEdgeContext2D &ctx,
    const ShapeSample &shape)
{
    if (ctx.edge_dofs == nullptr)
        throw std::runtime_error("EFGMIEdgeContext2D sem EdgeDofs associado.");

    std::unordered_map<int, int> local_index;
    local_index.reserve(shape.node_ids.size());
    for (size_t k = 0; k < shape.node_ids.size(); ++k)
        local_index.emplace(shape.node_ids[k], static_cast<int>(k));

    std::unordered_map<int, char> visited;
    visited.reserve(shape.node_ids.size() * 4);

    EdgeBasisSample out;
    for (const int node_id : shape.node_ids)
    {
        for (const int edge_id : ctx.node_to_edges[(size_t)node_id])
        {
            if (visited.find(edge_id) != visited.end())
                continue;
            visited.emplace(edge_id, 1);

            const Edge &edge = ctx.edge_dofs->edges[(size_t)edge_id];
            const int dof = ctx.edge_dofs->edge_to_dof[(size_t)edge_id];
            if (dof < 0)
                continue;

            const auto it0 = local_index.find(edge.n0);
            const auto it1 = local_index.find(edge.n1);
            if (it0 == local_index.end() || it1 == local_index.end())
                continue;

            const int i = it0->second;
            const int j = it1->second;
            const Node2D &ni = ctx.nodal.mesh->nodes[(size_t)edge.n0];
            const Node2D &nj = ctx.nodal.mesh->nodes[(size_t)edge.n1];
            const double dx = nj.x - ni.x;
            const double dy = nj.y - ni.y;
            const double L = std::sqrt(dx * dx + dy * dy);
            if (!(L > 0.0))
                continue;

            out.dof_ids.push_back(dof);
            out.vx.push_back(L * (shape.phi[i] * shape.dphidx[j] - shape.phi[j] * shape.dphidx[i]));
            out.vy.push_back(L * (shape.phi[i] * shape.dphidy[j] - shape.phi[j] * shape.dphidy[i]));
            out.curl.push_back(
                2.0 * L * (shape.dphidx[i] * shape.dphidy[j] - shape.dphidy[i] * shape.dphidx[j]));
        }
    }

    return out;
}

} // namespace efgmi2d
