#include <RcppEigen.h>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> hrbf_atoms_rcpp(
    const Eigen::Map<Eigen::MatrixXd> mask_xyz_world,
    const Eigen::Map<Eigen::MatrixXd> centres_xyz_world,
    const Eigen::Map<Eigen::VectorXd> sigma_vec_mm,
    std::string kernel_type,
    double value_threshold = 1e-8)
{
    const int N = mask_xyz_world.rows();
    const int K = centres_xyz_world.rows();
    if (sigma_vec_mm.size() != K) {
        Rcpp::stop("sigma_vec_mm length must match centres rows");
    }
    if (mask_xyz_world.cols() != 3) {
        Rcpp::stop("mask_xyz_world must have exactly 3 columns");
    }
    if (centres_xyz_world.cols() != 3) {
        Rcpp::stop("centres_xyz_world must have exactly 3 columns");
    }
    if (kernel_type != "gaussian" &&
        kernel_type != "wendland_c4" &&
        kernel_type != "wendland_c6") {
        Rcpp::stop("kernel_type must be one of: gaussian, wendland_c4, wendland_c6");
    }
    if (!std::isfinite(value_threshold) || value_threshold < 0.0) {
        Rcpp::stop("value_threshold must be finite and non-negative");
    }
    for (int k = 0; k < K; ++k) {
        if (sigma_vec_mm[k] <= 0.0) {
            Rcpp::stop("All sigma values must be positive; sigma_vec_mm[%d] = %g", k, sigma_vec_mm[k]);
        }
    }
    const bool gaussian = (kernel_type == "gaussian");

    std::vector< Eigen::Triplet<double> > triplets;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    std::vector< std::vector<Eigen::Triplet<double>> > thread_triplets(nthreads);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<Eigen::Triplet<double>>& local = thread_triplets[tid];
#pragma omp for schedule(static)
        for (int k = 0; k < K; ++k) {
            double sigma = sigma_vec_mm[k];
            double inv_two_sigma2 = 1.0 / (2.0 * sigma * sigma);
            double inv_sigma = 1.0 / sigma;
            bool wendland_c4 = (kernel_type == "wendland_c4");
            double gaussian_max_dist2 = std::numeric_limits<double>::infinity();
            if (gaussian && value_threshold > 0.0) {
                gaussian_max_dist2 = -2.0 * sigma * sigma * std::log(value_threshold);
            }
            for (int n = 0; n < N; ++n) {
                double dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
                double dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
                double dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
                double dist2 = dx*dx + dy*dy + dz*dz;
                double value;
                if (gaussian) {
                    if (dist2 > gaussian_max_dist2) continue;
                    value = std::exp(-dist2 * inv_two_sigma2);
                } else {
                    double dist = std::sqrt(dist2) * inv_sigma;
                    if (dist >= 1.0) continue;
                    double base = 1.0 - dist;
                    if (wendland_c4) {
                        value = std::pow(base,6.0) *
                            (35.0*dist*dist + 18.0*dist + 3.0) / 3.0;
                    } else {
                        value = std::pow(base,8.0) * (32.0*dist*dist*dist +
                                                      25.0*dist*dist + 8.0*dist + 1.0);
                    }
                }
                if (std::fabs(value) > value_threshold) {
                    local.emplace_back(k, n, value);
                }
            }
        }
    }
    for (auto &vec : thread_triplets) {
        triplets.insert(triplets.end(), vec.begin(), vec.end());
    }
#else
    triplets.reserve(static_cast<size_t>(N) * 4);
    for (int k = 0; k < K; ++k) {
        double sigma = sigma_vec_mm[k];
        double inv_two_sigma2 = 1.0 / (2.0 * sigma * sigma);
        double inv_sigma = 1.0 / sigma;
        bool wendland_c4 = (kernel_type == "wendland_c4");
        double gaussian_max_dist2 = std::numeric_limits<double>::infinity();
        if (gaussian && value_threshold > 0.0) {
            gaussian_max_dist2 = -2.0 * sigma * sigma * std::log(value_threshold);
        }
        for (int n = 0; n < N; ++n) {
            double dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
            double dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
            double dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
            double dist2 = dx*dx + dy*dy + dz*dz;
            double value;
            if (gaussian) {
                if (dist2 > gaussian_max_dist2) continue;
                value = std::exp(-dist2 * inv_two_sigma2);
            } else {
                double dist = std::sqrt(dist2) * inv_sigma;
                if (dist >= 1.0) continue;
                double base = 1.0 - dist;
                if (wendland_c4) {
                    value = std::pow(base,6.0) *
                        (35.0*dist*dist + 18.0*dist + 3.0) / 3.0;
                } else {
                    value = std::pow(base,8.0) * (32.0*dist*dist*dist +
                                                  25.0*dist*dist + 8.0*dist + 1.0);
                }
            }
            if (std::fabs(value) > value_threshold) {
                triplets.emplace_back(k, n, value);
            }
        }
    }
#endif
    Eigen::SparseMatrix<double> result(K, N);
    if (!triplets.empty()) {
        result.setFromTriplets(triplets.begin(), triplets.end());
    }
    return result;
}
