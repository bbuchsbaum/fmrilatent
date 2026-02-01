#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <algorithm>
#include <numeric>

// [[Rcpp::depends(RcppEigen)]]

//' Generate DPSS (Slepian) basis via dense prolate matrix eigen-decomposition
//'
//' @param n Integer length of the time series.
//' @param NW Time-bandwidth product (N * W), where W is normalized half-bandwidth (cycles per sample).
//' @param k Number of tapers to return (columns).
//' @return Matrix (n x k) with columns ordered by decreasing eigenvalue (energy concentration).
//'
//' @note
//' This implementation builds the dense prolate matrix
//'   A[i,j] = 2W (i==j) else sin(2*pi*W*(i-j)) / (pi*(i-j))
//' and computes the top-k eigenvectors. It is O(n^3) and intended as a
//' simple, dependency-light starting point; consider replacing with a
//' tridiagonal solver for very long n.
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd generate_dpss_basis_rcpp(const int n, const double NW, int k) {
  if (n <= 0) Rcpp::stop("n must be positive");
  if (NW <= 0.0) Rcpp::stop("NW must be positive");
  if (k <= 0) Rcpp::stop("k must be positive");
  double W = NW / static_cast<double>(n);
  if (W <= 0.0 || W >= 0.5) Rcpp::stop("normalized half-bandwidth W must be in (0, 0.5)");
  if (k > n) k = n;

  Eigen::MatrixXd A(n, n);
  const double two_pi_W = 2.0 * M_PI * W;
  for (int i = 0; i < n; ++i) {
    A(i, i) = 2.0 * W;
    for (int j = i + 1; j < n; ++j) {
      double diff = static_cast<double>(i - j);
      double val = std::sin(two_pi_W * diff) / (M_PI * diff);
      A(i, j) = val;
      A(j, i) = val;
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  if (es.info() != Eigen::Success) {
    Rcpp::stop("Eigen decomposition failed in generate_dpss_basis_rcpp");
  }

  Eigen::VectorXd evals = es.eigenvalues();
  Eigen::MatrixXd evecs = es.eigenvectors(); // columns correspond to evals (ascending)

  // sort indices by descending eigenvalue
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&evals](int a, int b) {
    return evals[a] > evals[b];
  });

  Eigen::MatrixXd out(n, k);
  for (int col = 0; col < k; ++col) {
    out.col(col) = evecs.col(idx[col]);
  }
  return out;
}

//' Generate DPSS via tridiagonal prolate matrix (O(n^2))
//'
//' Uses the classic tridiagonal formulation (Percival & Walden 1993; Thomson 1982)
//' and LAPACK dstev for efficiency on long series.
//'
//' @param n Integer length of the time series.
//' @param NW Time-bandwidth product (N * W), where W is normalized half-bandwidth.
//' @param k Number of tapers to return.
//' @return Matrix (n x k) with columns ordered by decreasing eigenvalue.
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd generate_dpss_tridiag_rcpp(const int n, const double NW, int k) {
  if (n <= 0) Rcpp::stop("n must be positive");
  if (NW <= 0.0) Rcpp::stop("NW must be positive");
  if (k <= 0) Rcpp::stop("k must be positive");
  double W = NW / static_cast<double>(n);
  if (W <= 0.0 || W >= 0.5) Rcpp::stop("normalized half-bandwidth W must be in (0, 0.5)");
  if (k > n) k = n;

  // Tridiagonal coefficients
  Eigen::VectorXd d(n);         // main diagonal
  Eigen::VectorXd e(n - 1);     // off-diagonal

  for (int i = 0; i < n; ++i) {
    double term = (static_cast<double>(n - 1 - 2 * i) / 2.0);
    d[i] = std::cos(2.0 * M_PI * W) * term * term;
  }
  for (int i = 0; i < n - 1; ++i) {
    double ii = static_cast<double>(i + 1);          // 1-based in formula
    e[i] = (ii * (static_cast<double>(n) - ii)) / 2.0;
  }

  // LAPACK dstev expects raw arrays; computes all eigenpairs of symmetric tridiagonal
  Eigen::VectorXd evals = d;
  Eigen::VectorXd off = e;
  Eigen::MatrixXd z(n, n);
  std::vector<double> work(std::max(1, 2 * n - 2));
  char jobz = 'V';
  int nn = n;
  int ldz = n;
  int info = 0;
  F77_CALL(dstev)(&jobz, &nn, evals.data(), off.data(),
                  z.data(), &ldz, work.data(), &info FCONE);
  if (info != 0) {
    Rcpp::stop("dstev failed in generate_dpss_tridiag_rcpp (info=%d)", info);
  }

  // sort descending eigenvalues
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&evals](int a, int b) {
    return evals[a] > evals[b];
  });

  Eigen::MatrixXd out(n, k);
  for (int col = 0; col < k; ++col) {
    out.col(col) = z.col(idx[col]);
  }
  return out;
}
