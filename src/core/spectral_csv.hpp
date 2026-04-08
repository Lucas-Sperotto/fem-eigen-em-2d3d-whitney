/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/spectral_csv.hpp                                         */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exportacao didatica de matrizes e espectros em CSV.             */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: apoio operacional   */
/* aos executaveis das secoes 2.1 a 2.2.4.                                    */
/*****************************************************************************/

#pragma once

#include "dense.hpp"
#include "lapack_eig.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace spectral_csv
{

struct GeneralEigenColumnBinding
{
    int original_index = 0;
    double lambda_real = 0.0;
    double lambda_imag = 0.0;
    double beta_denominator = 0.0;
    int real_column = 0;
    int imag_column = -1;
    int imag_sign = 0;
};

template <typename Accessor>
inline bool write_crs_csv(
    const std::filesystem::path &path,
    int nr,
    int nc,
    Accessor accessor,
    double abs_tol = 1.0e-14)
{
    if (nr < 0 || nc < 0)
        return false;

    std::vector<int> row_ptr((size_t)nr + 1, 0);
    std::vector<int> col_idx;
    std::vector<double> values;
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            const double value = accessor(i, j);
            if (std::abs(value) <= abs_tol)
                continue;
            col_idx.push_back(j);
            values.push_back(value);
        }
        row_ptr[(size_t)i + 1] = static_cast<int>(col_idx.size());
    }

    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "section,index,value\n";
    out << "nrows,0," << nr << "\n";
    out << "ncols,0," << nc << "\n";
    out << "nnz,0," << values.size() << "\n";
    for (size_t i = 0; i < row_ptr.size(); ++i)
        out << "row_ptr," << i << "," << row_ptr[i] << "\n";
    for (size_t i = 0; i < col_idx.size(); ++i)
        out << "col_idx," << i << "," << col_idx[i] << "\n";
    for (size_t i = 0; i < values.size(); ++i)
        out << "values," << i << "," << values[i] << "\n";
    return static_cast<bool>(out);
}

inline bool write_dense_crs_csv(
    const std::filesystem::path &path,
    const DenseMat &mat,
    double abs_tol = 1.0e-14)
{
    return write_crs_csv(
        path,
        mat.n,
        mat.n,
        [&](int i, int j)
        { return mat(i, j); },
        abs_tol);
}

template <typename Accessor>
inline bool write_rect_like_crs_csv(
    const std::filesystem::path &path,
    int nr,
    int nc,
    Accessor accessor,
    double abs_tol = 1.0e-14)
{
    return write_crs_csv(path, nr, nc, accessor, abs_tol);
}

inline std::vector<int> symmetric_eigen_order(const GenEigResult &res)
{
    std::vector<int> order((size_t)res.n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b)
              {
                  if (res.w[(size_t)a] != res.w[(size_t)b])
                      return res.w[(size_t)a] < res.w[(size_t)b];
                  return a < b;
              });
    return order;
}

inline bool write_symmetric_eigenvalues_csv(
    const std::filesystem::path &path,
    const GenEigResult &res)
{
    const auto order = symmetric_eigen_order(res);
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "ordered_rank,solver_index,lambda\n";
    for (size_t rank = 0; rank < order.size(); ++rank)
    {
        const int solver_index = order[rank];
        out << (rank + 1) << ","
            << solver_index << ","
            << res.w[(size_t)solver_index] << "\n";
    }
    return static_cast<bool>(out);
}

inline bool write_symmetric_eigenvectors_csv(
    const std::filesystem::path &path,
    const GenEigResult &res)
{
    const auto order = symmetric_eigen_order(res);
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "ordered_rank,solver_index,dof_index,component_real,component_imag\n";
    for (size_t rank = 0; rank < order.size(); ++rank)
    {
        const int solver_index = order[rank];
        for (int dof = 0; dof < res.n; ++dof)
        {
            const double value = res.Zcol[(size_t)solver_index * res.n + dof];
            out << (rank + 1) << ","
                << solver_index << ","
                << dof << ","
                << value << ",0\n";
        }
    }
    return static_cast<bool>(out);
}

inline std::vector<GeneralEigenColumnBinding> general_eigen_bindings(
    const GenEigGeneralResult &res,
    double imag_tol = 1.0e-12)
{
    std::vector<GeneralEigenColumnBinding> out;
    out.reserve((size_t)res.n);

    for (int i = 0; i < res.n; ++i)
    {
        const double imag = res.lambda_im[(size_t)i];
        if (std::abs(imag) <= imag_tol)
        {
            out.push_back(
                {
                    i,
                    res.lambda_re[(size_t)i],
                    0.0,
                    res.beta[(size_t)i],
                    i,
                    -1,
                    0,
                });
            continue;
        }

        if (imag > 0.0 && (i + 1) < res.n)
        {
            out.push_back(
                {
                    i,
                    res.lambda_re[(size_t)i],
                    res.lambda_im[(size_t)i],
                    res.beta[(size_t)i],
                    i,
                    i + 1,
                    +1,
                });
            out.push_back(
                {
                    i + 1,
                    res.lambda_re[(size_t)i + 1],
                    res.lambda_im[(size_t)i + 1],
                    res.beta[(size_t)i + 1],
                    i,
                    i + 1,
                    -1,
                });
            ++i;
            continue;
        }

        out.push_back(
            {
                i,
                res.lambda_re[(size_t)i],
                res.lambda_im[(size_t)i],
                res.beta[(size_t)i],
                i,
                -1,
                0,
            });
    }

    std::sort(out.begin(), out.end(), [](const GeneralEigenColumnBinding &a, const GeneralEigenColumnBinding &b)
              {
                  if (a.lambda_real != b.lambda_real)
                      return a.lambda_real < b.lambda_real;
                  if (a.lambda_imag != b.lambda_imag)
                      return a.lambda_imag < b.lambda_imag;
                  return a.original_index < b.original_index;
              });
    return out;
}

inline bool write_general_eigenvalues_csv(
    const std::filesystem::path &path,
    const GenEigGeneralResult &res,
    double imag_tol = 1.0e-12)
{
    const auto bindings = general_eigen_bindings(res, imag_tol);
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "ordered_rank,solver_index,lambda_real,lambda_imag,beta_denominator\n";
    for (size_t rank = 0; rank < bindings.size(); ++rank)
    {
        const auto &item = bindings[rank];
        out << (rank + 1) << ","
            << item.original_index << ","
            << item.lambda_real << ","
            << item.lambda_imag << ","
            << item.beta_denominator << "\n";
    }
    return static_cast<bool>(out);
}

inline bool write_general_eigenvectors_csv(
    const std::filesystem::path &path,
    const GenEigGeneralResult &res,
    double imag_tol = 1.0e-12)
{
    const auto bindings = general_eigen_bindings(res, imag_tol);
    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(16);
    out << "ordered_rank,solver_index,dof_index,component_real,component_imag\n";
    for (size_t rank = 0; rank < bindings.size(); ++rank)
    {
        const auto &item = bindings[rank];
        for (int dof = 0; dof < res.n; ++dof)
        {
            const double real_part = res.VRcol[(size_t)item.real_column * res.n + dof];
            double imag_part = 0.0;
            if (item.imag_column >= 0 && item.imag_sign != 0)
                imag_part = item.imag_sign * res.VRcol[(size_t)item.imag_column * res.n + dof];
            out << (rank + 1) << ","
                << item.original_index << ","
                << dof << ","
                << real_part << ","
                << imag_part << "\n";
        }
    }
    return static_cast<bool>(out);
}

inline bool write_symmetric_problem_exports(
    const std::filesystem::path &dir,
    const std::string &prefix,
    const DenseMat &S,
    const DenseMat &T,
    const GenEigResult &res,
    double abs_tol = 1.0e-14)
{
    return write_dense_crs_csv(dir / (prefix + "_S_crs.csv"), S, abs_tol) &&
           write_dense_crs_csv(dir / (prefix + "_T_crs.csv"), T, abs_tol) &&
           write_symmetric_eigenvalues_csv(dir / (prefix + "_eigenvalues.csv"), res) &&
           write_symmetric_eigenvectors_csv(dir / (prefix + "_eigenvectors.csv"), res);
}

inline bool write_general_problem_exports_named(
    const std::filesystem::path &dir,
    const std::string &prefix,
    const std::string &left_label,
    const std::string &right_label,
    const DenseMat &left_mat,
    const DenseMat &right_mat,
    const GenEigGeneralResult &res,
    double abs_tol = 1.0e-14,
    double imag_tol = 1.0e-12)
{
    return write_dense_crs_csv(dir / (prefix + "_" + left_label + "_crs.csv"), left_mat, abs_tol) &&
           write_dense_crs_csv(dir / (prefix + "_" + right_label + "_crs.csv"), right_mat, abs_tol) &&
           write_general_eigenvalues_csv(dir / (prefix + "_eigenvalues.csv"), res, imag_tol) &&
           write_general_eigenvectors_csv(dir / (prefix + "_eigenvectors.csv"), res, imag_tol);
}

inline bool write_general_problem_exports(
    const std::filesystem::path &dir,
    const std::string &prefix,
    const DenseMat &A,
    const DenseMat &B,
    const GenEigGeneralResult &res,
    double abs_tol = 1.0e-14,
    double imag_tol = 1.0e-12)
{
    return write_general_problem_exports_named(
        dir,
        prefix,
        "A",
        "B",
        A,
        B,
        res,
        abs_tol,
        imag_tol);
}

template <typename Accessor>
inline void write_crs_with_warning(
    const std::filesystem::path &path,
    int nr,
    int nc,
    Accessor accessor,
    const std::string &label,
    double abs_tol = 1.0e-14)
{
    if (!write_crs_csv(path, nr, nc, accessor, abs_tol))
        throw std::runtime_error("Falha ao escrever matriz CRS: " + label);
}

} // namespace spectral_csv
