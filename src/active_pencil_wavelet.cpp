#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// --- 1D CDF 5/3 lifting ----------------------------------------------------
inline void lift_1d_cdf53(std::vector<double>& data, int levels, bool forward) {
  int n = data.size();
  if (n < 2 || levels <= 0) return;

  auto forward_once = [&](std::vector<double>& v) {
    int nloc = v.size();
    int half = (nloc + 1) / 2;
    std::vector<double> s(half);
    std::vector<double> d(nloc - half);

    for (int i = 0; i < half; ++i) s[i] = v[2 * i];
    for (int i = 0; i < (nloc - half); ++i) d[i] = v[2 * i + 1];

    for (int i = 0; i < (nloc - half); ++i) {
      double s1 = s[i];
      double s2 = (i + 1 < half) ? s[i + 1] : s[half - 1];
      d[i] -= 0.5 * (s1 + s2);
    }
    for (int i = 0; i < half; ++i) {
      double d1 = (i > 0) ? d[i - 1] : d[0];
      double d2 = (i < (nloc - half)) ? d[i] : d[(nloc - half) - 1];
      s[i] += 0.25 * (d1 + d2);
    }
    for (int i = 0; i < half; ++i) v[i] = s[i];
    for (int i = 0; i < (nloc - half); ++i) v[half + i] = d[i];
  };

  auto inverse_once = [&](std::vector<double>& v) {
    int nloc = v.size();
    int half = (nloc + 1) / 2;
    std::vector<double> s(half);
    std::vector<double> d(nloc - half);
    for (int i = 0; i < half; ++i) s[i] = v[i];
    for (int i = 0; i < (nloc - half); ++i) d[i] = v[half + i];
    for (int i = 0; i < half; ++i) {
      double d1 = (i > 0) ? d[i - 1] : d[0];
      double d2 = (i < (nloc - half)) ? d[i] : d[(nloc - half) - 1];
      s[i] -= 0.25 * (d1 + d2);
    }
    for (int i = 0; i < (nloc - half); ++i) {
      double s1 = s[i];
      double s2 = (i + 1 < half) ? s[i + 1] : s[half - 1];
      d[i] += 0.5 * (s1 + s2);
    }
    for (int i = 0; i < half; ++i) v[2 * i] = s[i];
    for (int i = 0; i < (nloc - half); ++i) v[2 * i + 1] = d[i];
  };

  if (forward) {
    int curr = n;
    for (int l = 0; l < levels && curr >= 2; ++l) {
      std::vector<double> sub(data.begin(), data.begin() + curr);
      forward_once(sub);
      for (int i = 0; i < curr; ++i) data[i] = sub[i];
      curr = (curr + 1) / 2;
    }
  } else {
    std::vector<int> sizes;
    int curr = n;
    for (int l = 0; l < levels && curr >= 2; ++l) {
      sizes.push_back(curr);
      curr = (curr + 1) / 2;
    }
    for (int idx = static_cast<int>(sizes.size()) - 1; idx >= 0; --idx) {
      int sz = sizes[idx];
      std::vector<double> sub(data.begin(), data.begin() + sz);
      inverse_once(sub);
      for (int i = 0; i < sz; ++i) data[i] = sub[i];
    }
  }
}

// [[Rcpp::export]]
NumericMatrix cdf53_time_lift(NumericMatrix X, int levels, bool forward) {
  NumericMatrix out = clone(X);
  int nrow = out.nrow();
  int ncol = out.ncol();
  if (levels <= 0) return out;
  std::vector<double> buf(nrow);
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) buf[i] = out(i, j);
    lift_1d_cdf53(buf, levels, forward);
    for (int i = 0; i < nrow; ++i) out(i, j) = buf[i];
  }
  return out;
}

inline int flat_idx(int x, int y, int z, int nx, int ny) {
  return x + nx * (y + ny * z);
}

// [[Rcpp::export]]
NumericVector active_pencil_wavelet(NumericVector data_voxels,
                                    IntegerMatrix coords,
                                    IntegerVector dims,
                                    int levels,
                                    bool forward) {
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int total_vol = nx * ny * nz;
  int n_vox = coords.nrow();
  if (data_voxels.size() != n_vox) stop("Data size mismatch");

  std::vector<int> grid(total_vol, -1);
  for (int i = 0; i < n_vox; ++i) {
    int x = coords(i, 0) - 1;
    int y = coords(i, 1) - 1;
    int z = coords(i, 2) - 1;
    grid[flat_idx(x, y, z, nx, ny)] = i;
  }

  std::vector<double> buffer = as<std::vector<double>>(data_voxels);
  int max_dim = std::max({nx, ny, nz});
  std::vector<double> pencil(max_dim);

  auto process_axis = [&](int axis) {
    int a = (axis == 0) ? nx : (axis == 1 ? ny : nz);
    int b = (axis == 0) ? ny : (axis == 1 ? nz : nx);
    int c = (axis == 0) ? nz : (axis == 1 ? nx : ny);
    for (int c1 = 0; c1 < c; ++c1) {
      for (int b1 = 0; b1 < b; ++b1) {
        bool active = false;
        for (int a1 = 0; a1 < a; ++a1) {
          int idx;
          if (axis == 0)
            idx = grid[flat_idx(a1, b1, c1, nx, ny)];
          else if (axis == 1)
            idx = grid[flat_idx(c1, a1, b1, nx, ny)];
          else
            idx = grid[flat_idx(b1, c1, a1, nx, ny)];
          if (idx != -1) { active = true; break; }
        }
        if (!active) continue;

        for (int a1 = 0; a1 < a; ++a1) {
          int idx;
          if (axis == 0)
            idx = grid[flat_idx(a1, b1, c1, nx, ny)];
          else if (axis == 1)
            idx = grid[flat_idx(c1, a1, b1, nx, ny)];
          else
            idx = grid[flat_idx(b1, c1, a1, nx, ny)];
          pencil[a1] = (idx != -1) ? buffer[idx] : 0.0;
        }

        // operate only on first a entries
        std::vector<double> sub(a);
        for (int k = 0; k < a; ++k) sub[k] = pencil[k];
        lift_1d_cdf53(sub, levels, forward);
        for (int k = 0; k < a; ++k) pencil[k] = sub[k];

        for (int a1 = 0; a1 < a; ++a1) {
          int idx;
          if (axis == 0)
            idx = grid[flat_idx(a1, b1, c1, nx, ny)];
          else if (axis == 1)
            idx = grid[flat_idx(c1, a1, b1, nx, ny)];
          else
            idx = grid[flat_idx(b1, c1, a1, nx, ny)];
          if (idx != -1) buffer[idx] = pencil[a1];
        }
      }
    }
  };

  if (forward) {
    process_axis(0);
    process_axis(1);
    process_axis(2);
  } else {
    process_axis(2);
    process_axis(1);
    process_axis(0);
  }

  return wrap(buffer);
}
