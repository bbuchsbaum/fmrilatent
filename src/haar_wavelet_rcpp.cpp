#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {

inline size_t lin3D(int i, int j, int k, int nx, int ny, int nz) {
  return static_cast<size_t>(i) + static_cast<size_t>(nx) *
         (static_cast<size_t>(j) + static_cast<size_t>(ny) * static_cast<size_t>(k));
}

inline int bitCeilLog2(int v) {
  if (v <= 1) {
    return 0;
  }
  int p = 0;
  --v;
  while (v > 0) {
    v >>= 1;
    ++p;
  }
  return p;
}

inline unsigned long long mortonEncode3D(unsigned int x, unsigned int y,
                                         unsigned int z, int bits) {
  unsigned long long code = 0ULL;
  for (int b = 0; b < bits; ++b) {
    const unsigned long long xb = (x >> b) & 1U;
    const unsigned long long yb = (y >> b) & 1U;
    const unsigned long long zb = (z >> b) & 1U;
    code |= (xb << (3 * b + 0)) | (yb << (3 * b + 1)) | (zb << (3 * b + 2));
  }
  return code;
}

struct MortonEntry {
  unsigned long long code;
  int x;
  int y;
  int z;
};

}  // namespace

// [[Rcpp::export]]
IntegerVector get_morton_ordered_indices_rcpp(LogicalVector mask) {
  SEXP dimAttr = mask.attr("dim");
  if (dimAttr == R_NilValue) {
    stop("mask must be a 3D array");
  }
  IntegerVector dims(dimAttr);
  if (dims.size() != 3) {
    stop("mask must be a 3D array");
  }
  const int nx = dims[0], ny = dims[1], nz = dims[2];

  const int bits = bitCeilLog2(std::max(nx, std::max(ny, nz)));
  std::vector<MortonEntry> entries;
  entries.reserve(mask.size());

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (mask[lin3D(i, j, k, nx, ny, nz)] == TRUE) {
          unsigned long long code = mortonEncode3D(static_cast<unsigned int>(i),
                                                   static_cast<unsigned int>(j),
                                                   static_cast<unsigned int>(k),
                                                   bits);
          entries.push_back({code, i, j, k});
        }
      }
    }
  }

  std::sort(entries.begin(), entries.end(), [](const MortonEntry& a,
                                                const MortonEntry& b) {
    if (a.code != b.code) return a.code < b.code;
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
  });

  IntegerVector out(entries.size());
  for (size_t idx = 0; idx < entries.size(); ++idx) {
    const MortonEntry& e = entries[idx];
    out[idx] = static_cast<int>(lin3D(e.x, e.y, e.z, nx, ny, nz)) + 1;
  }
  return out;
}

// [[Rcpp::export]]
List precompute_haar_scalings_rcpp(LogicalVector mask, int levels) {
  SEXP dimAttr = mask.attr("dim");
  if (dimAttr == R_NilValue) {
    stop("mask must be a 3D array");
  }
  IntegerVector dims(dimAttr);
  if (dims.size() != 3) {
    stop("mask must be a 3D array");
  }

  int nx = dims[0], ny = dims[1], nz = dims[2];
  std::vector<unsigned char> occ(mask.size());
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        occ[lin3D(i, j, k, nx, ny, nz)] =
            mask[lin3D(i, j, k, nx, ny, nz)] == TRUE ? 1U : 0U;
      }
    }
  }

  List scalings(levels);
  int cx = nx, cy = ny, cz = nz;
  for (int lvl = 0; lvl < levels; ++lvl) {
    const int nbx = (cx + 1) / 2;
    const int nby = (cy + 1) / 2;
    const int nbz = (cz + 1) / 2;
    const int nblocks = nbx * nby * nbz;

    NumericVector sqrt_nvalid(nblocks);
    NumericVector sqrt_nvalid_div8(nblocks);
    std::vector<unsigned char> next_occ(nblocks, 0U);

    int idx = 0;
    for (int bi = 0; bi < nbx; ++bi) {
      const int x0 = 2 * bi;
      const int x1 = std::min(x0 + 1, cx - 1);
      for (int bj = 0; bj < nby; ++bj) {
        const int y0 = 2 * bj;
        const int y1 = std::min(y0 + 1, cy - 1);
        for (int bk = 0; bk < nbz; ++bk, ++idx) {
          const int z0 = 2 * bk;
          const int z1 = std::min(z0 + 1, cz - 1);
          int cnt = 0;
          for (int z = z0; z <= z1; ++z) {
            for (int y = y0; y <= y1; ++y) {
              for (int x = x0; x <= x1; ++x) {
                cnt += occ[lin3D(x, y, z, cx, cy, cz)];
              }
            }
          }
          if (cnt > 0) {
            next_occ[idx] = 1U;
          }
          sqrt_nvalid[idx] = std::sqrt(static_cast<double>(cnt));
          sqrt_nvalid_div8[idx] = std::sqrt(static_cast<double>(cnt) / 8.0);
        }
      }
    }

    List entry;
    entry["sqrt_nvalid"] = sqrt_nvalid;
    entry["sqrt_nvalid_div_8"] = sqrt_nvalid_div8;
    scalings[lvl] = entry;

    occ.swap(next_occ);
    cx = nbx;
    cy = nby;
    cz = nbz;
  }

  return scalings;
}

// [[Rcpp::export]]
List forward_lift_rcpp(NumericVector data_morton,
                       LogicalVector mask_flat_morton,
                       IntegerVector mask_dims,
                       int levels,
                       List scalings) {
  (void)mask_flat_morton;
  (void)mask_dims;
  std::vector<NumericVector> sqrtN(levels), sqrtN_div8(levels);
  std::vector<std::vector<int> > counts(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    List sc = scalings[lvl];
    NumericVector s = sc["sqrt_nvalid"];
    NumericVector s8 = sc["sqrt_nvalid_div_8"];
    sqrtN[lvl] = s;
    sqrtN_div8[lvl] = s8;
    std::vector<int> c(s.size());
    for (int i = 0; i < s.size(); ++i) {
      const double v = s[i];
      c[i] = static_cast<int>(std::llround(v * v));
    }
    counts[lvl].swap(c);
  }

  NumericVector current = clone(data_morton);
  std::vector<NumericVector> details(levels);

  for (int lvl = 0; lvl < levels; ++lvl) {
    const std::vector<int>& cnt = counts[lvl];
    const int numBlocks = static_cast<int>(cnt.size());

    int present = 0;
    for (int nv : cnt) if (nv > 0) ++present;
    NumericVector next_data(present);

    size_t total = 0;
    for (int nv : cnt) total += static_cast<size_t>(nv);
    NumericVector dvec(total);

    int idx_in = 0;
    size_t idx_out = 0;
    int idx_lp = 0;

    for (int b = 0; b < numBlocks; ++b) {
      const int nv = cnt[b];
      if (nv > 0) {
        double sum = 0.0;
        for (int t = 0; t < nv; ++t) sum += current[idx_in + t];
        const double avg = sum / static_cast<double>(nv);
        next_data[idx_lp++] = avg * sqrtN[lvl][b];
        const double s8 = sqrtN_div8[lvl][b];
        for (int t = 0; t < nv; ++t) dvec[idx_out + t] = (current[idx_in + t] - avg) * s8;
        idx_in += nv;
        idx_out += static_cast<size_t>(nv);
      }
    }

    details[lvl] = dvec;
    current = next_data;
  }

  return List::create(_["root_coeff"] = current,
                      _["detail_coeffs_by_level"] = details);
}

// [[Rcpp::export]]
NumericVector inverse_lift_rcpp(double root_coeff,
                                List detail_vecs,
                                LogicalVector mask_flat_morton,
                                IntegerVector mask_dims,
                                int levels,
                                List scalings) {
  (void)mask_flat_morton;
  (void)mask_dims;
  std::vector<NumericVector> sqrtN(levels), sqrtN_div8(levels);
  std::vector<std::vector<int> > counts(levels);
  for (int lvl = 0; lvl < levels; ++lvl) {
    List sc = scalings[lvl];
    NumericVector s = sc["sqrt_nvalid"];
    NumericVector s8 = sc["sqrt_nvalid_div_8"];
    sqrtN[lvl] = s;
    sqrtN_div8[lvl] = s8;
    std::vector<int> c(s.size());
    for (int i = 0; i < s.size(); ++i) {
      const double v = s[i];
      c[i] = static_cast<int>(std::llround(v * v));
    }
    counts[lvl].swap(c);
  }

  NumericVector current(1);
  current[0] = root_coeff;

  for (int lvl = levels - 1; lvl >= 0; --lvl) {
    const std::vector<int>& cnt = counts[lvl];
    const int numBlocks = static_cast<int>(cnt.size());

    size_t total = 0; for (int nv : cnt) total += static_cast<size_t>(nv);
    NumericVector next_data(total);

    NumericVector dvec = detail_vecs[lvl];
    if (static_cast<size_t>(dvec.size()) != total) {
      stop("Inverse lift detail length mismatch.");
    }

    size_t idx_out = 0;
    size_t idx_det = 0;
    int idx_lp = 0;
    for (int b = 0; b < numBlocks; ++b) {
      const int nv = cnt[b];
      if (nv > 0) {
        if (idx_lp >= static_cast<int>(current.size())) stop("Inverse lowpass underflow.");
        const double avg = current[idx_lp++] / sqrtN[lvl][b];
        const double s8 = sqrtN_div8[lvl][b];
        for (int t = 0; t < nv; ++t) next_data[idx_out + t] = dvec[idx_det + t] / s8 + avg;
        idx_out += static_cast<size_t>(nv);
        idx_det += static_cast<size_t>(nv);
      }
    }
    current = next_data;
  }

  return current;
}

// [[Rcpp::export]]
IntegerVector get_valid_finest_blocks_rcpp(LogicalVector mask) {
  SEXP dimAttr = mask.attr("dim");
  if (dimAttr == R_NilValue) {
    stop("mask must be a 3D array");
  }
  IntegerVector dims(dimAttr);
  if (dims.size() != 3) {
    stop("mask must be a 3D array");
  }

  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const int nbx = (nx + 1) / 2;
  const int nby = (ny + 1) / 2;
  const int nbz = (nz + 1) / 2;
  const int bits = bitCeilLog2(std::max(nbx, std::max(nby, nbz)));

  std::vector<unsigned long long> codes;
  codes.reserve(static_cast<size_t>(nbx) * nby * nbz);

  for (int bi = 0; bi < nbx; ++bi) {
    const int x0 = 2 * bi;
    const int x1 = std::min(x0 + 1, nx - 1);
    for (int bj = 0; bj < nby; ++bj) {
      const int y0 = 2 * bj;
      const int y1 = std::min(y0 + 1, ny - 1);
      for (int bk = 0; bk < nbz; ++bk) {
        const int z0 = 2 * bk;
        const int z1 = std::min(z0 + 1, nz - 1);
        bool any = false;
        for (int z = z0; z <= z1 && !any; ++z) {
          for (int y = y0; y <= y1 && !any; ++y) {
            for (int x = x0; x <= x1; ++x) {
              if (mask[lin3D(x, y, z, nx, ny, nz)] == TRUE) {
                any = true;
                break;
              }
            }
          }
        }
        if (any) {
          codes.push_back(
              mortonEncode3D(static_cast<unsigned int>(bi),
                             static_cast<unsigned int>(bj),
                             static_cast<unsigned int>(bk), bits));
        }
      }
    }
  }

  std::sort(codes.begin(), codes.end());
  IntegerVector out(codes.size());
  for (size_t idx = 0; idx < codes.size(); ++idx) {
    out[idx] = static_cast<int>(codes[idx]);
  }
  return out;
}

