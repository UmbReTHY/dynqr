#ifndef DYNQR_H_
#define DYNQR_H_

#include <cassert>
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <numeric>
#include <vector>

namespace dynqr {

template <typename Float_>
class dynamic_qr {
 public:
  typedef Float_ Float;

  /**
    @brief will create a dynamic qr for at most num_rows columns
  */
  dynamic_qr(int num_rows)
      : _q(eigen_mat::Identity(num_rows, num_rows)),
        _r_raw(new Float[num_rows * num_rows]),
        _r_begin(new Float*[num_rows]),
        _r_end(_r_begin) {
    DLOG(INFO) << "**QR-CTOR " << this << std::endl;
    CHECK(num_rows > 0);
    std::vector<int> offsets(num_rows);
    std::iota(offsets.begin(), offsets.end(), 0);
    init_row_ptr(offsets.begin());
  }

  // copy-constructor
  dynamic_qr(dynamic_qr const& orig)
      // note: orig.num_rows() == maximum number of columns possible
      : _q(orig._q),
        _r_raw(nullptr),
        _r_begin(new Float*[orig.num_rows()]),
        _r_end(_r_begin + orig.num_cols()) {
    DLOG(INFO) << "**QR-COPY-CTOR " << this << std::endl;
    int const num_elements = orig.num_rows() * orig.num_rows();
    _r_raw = new Float[num_elements];

    // since some column pointers might have been swapped,
    // we need to copy the exact same offsets
    std::vector<int> offsets(num_rows());
    for (int i = 0; i < offsets.size(); ++i)
      offsets[i] = (orig._r_begin[i] - orig._r_raw) / orig.num_rows();
    init_row_ptr(offsets.begin());

    std::copy(orig._r_raw, orig._r_raw + num_elements, _r_raw);

    CHECK(num_rows() == orig.num_rows());
    CHECK(num_cols() == orig.num_cols());
    CHECK(_q == orig._q);
    DCHECK(std::equal(_r_raw, _r_raw + num_elements, orig._r_raw));
  }

  // move-constructor
  dynamic_qr(dynamic_qr&& tmp)
      : _q(), _r_raw(tmp._r_raw), _r_begin(tmp._r_begin), _r_end(tmp._r_end) {
    DLOG(INFO) << "**QR-MOVE-CTOR " << this << std::endl;
    _q.swap(tmp._q);
    // ownership has been transferred
    tmp._r_raw = nullptr;
    tmp._r_begin = nullptr;
  }

  // move-assignment
  dynamic_qr& operator=(dynamic_qr&& tmp) {
    DLOG(INFO) << "**QR-MOVE-ASSIGN " << this << std::endl;
    if (this != &tmp) {
      _q.swap(tmp._q);
      _r_raw = tmp._r_raw;
      _r_begin = tmp._r_begin;
      _r_end = tmp._r_end;
      // ownership has been transferred
      tmp._r_raw = nullptr;
      tmp._r_begin = nullptr;
    }
    return *this;
  }

  dynamic_qr& operator=(dynamic_qr const&) = delete;

  int num_rows() const noexcept {
    CHECK(_q.rows() >= 0);
    return _q.rows();
  }

  int num_cols() const noexcept {
    CHECK(_r_end >= _r_begin);
    return _r_end - _r_begin;
  }

  template <typename Vector>
  void append_column(const Vector& col) {
    CHECK(num_cols() < num_rows());
    eigen_map new_col(*_r_end, num_rows());
    new_col = _q.transpose() * col;
    Float c, s;
    for (int i = num_rows() - 1; i-- > num_cols();) {
      // i is the index of 'a' in the givens routine
      givens(new_col[i], new_col[i + 1], &c, &s);  // c and s
      update_q(c, s, i);
      apply_half_givens(c, s, *_r_end, i);  // update R
      // new_col[i + 1] is implicitely set to 0
    }
    ++_r_end;
  }

  void delete_column(int pos) {
    CHECK(pos >= 0);
    CHECK(pos < num_cols());
    if (pos != num_cols() - 1) {
      // move all columns right of pos one to the left
      Float* tmp = _r_begin[pos];
      for (int i = pos, max = num_cols() - 1; i < max; ++i)
        _r_begin[i] = _r_begin[i + 1];
      --_r_end;       // num_cols() is now up-to-date again
      *_r_end = tmp;  // restored the address of the deleted column
      // R is now upper Hessenberg starting at column pos
      hessenberg_update(pos, num_cols());
    } else {
      --_r_end;
    }
  }

  /**
    @brief updates the QR decomposition for A + u * v^T
  */
  template <typename Vector1, typename Vector2>
  void rank_one_update(const Vector1& u, const Vector2& v) {
    CHECK(int(u.rows()) == num_rows());
    CHECK(int(v.rows()) == num_cols());

    // we now need to actually zero the subdiagonal entries, since we
    // access those during the updates of R
    for (int i = 0; i < (num_cols() - 1); ++i) _r_begin[i][i + 1] = 0;
    bool const is_not_square = num_rows() > num_cols();
    if (is_not_square)  // only the the last column has a sub-diagonal element
      _r_begin[num_cols() - 1][num_cols()] = 0;

    eigen_vector w = _q.transpose() * u;
    Float c, s;
    for (int i = num_rows() - 1; i-- > 0;) {
      givens(w[i], w[i + 1], &c, &s);
      apply_half_givens(c, s, w.data(), i);  // clear element i + 1 of w
      if (i < num_cols()) apply_givens(c, s, i, i);
      update_q(c, s, i);
    }
    // compute H = R + w * v^T
    // Since w = +-||w|| * e_1, this amounts to adding w[0] * v^t to R.row(0)
    for (int i = 0; i < num_cols(); ++i) _r_begin[i][0] += w[0] * v[i];
    // R is now upper Hessenberg
    hessenberg_update(0, (is_not_square ? num_cols() : num_cols() - 1));
  }

  /**
   @brief finds the least-squares solution to QR * x = b
 */
  template <typename Vector1, typename Vector2>
  void solve(const Vector1& b, Vector2& x) const {
    CHECK(num_cols() > 0);
    CHECK(int(x.rows()) == num_cols());
    eigen_vector u = _q.transpose() * b;
    // now back substitution: R * x = u
    for (int i = num_cols(); i-- > 0;) {
      for (int j = i + 1; j < num_cols(); ++j) u[i] -= x[j] * _r_begin[j][i];
      x[i] = u[i] / _r_begin[i][i];
    }
  }

  ~dynamic_qr() {
    DLOG(INFO) << "**QR-DESTRUCT " << this << std::endl;
    if (_r_begin) {
      delete[] _r_begin;
      _r_begin = nullptr;
    }
    if (_r_raw) {
      delete[] _r_raw;
      _r_raw = nullptr;
    }
  }

 private:
  /**
    @brief helper function to avoid code duplication in constructors
    @param begin iterator to the offsets of the columns into the raw pointer
  */
  template <class Iterator>
  void init_row_ptr(Iterator begin) {
    int const max_cols = _q.cols();
    auto const& num_rows = max_cols;
    CHECK(nullptr != _r_raw);  // to make sure constructors initialized
                               // the storage
    for (int i = 0; i < max_cols; ++i)
      _r_begin[i] = _r_raw + *(begin++) * num_rows;
  }

  void givens(const Float& a, const Float& b, Float* c, Float* s) const {
    if (0.0 == b) {
      *c = 1.0;
      *s = 0.0;
    } else {
      if (abs(b) > abs(a)) {
        const Float tau = -a / b;
        *s = 1 / hypot(1, tau);
        *c = *s * tau;
      } else {
        const Float tau = -b / a;
        *c = 1 / hypot(1, tau);
        *s = *c * tau;
      }
    }
    CHECK(abs(*c) <= 1) << "c = " << c;
    CHECK(abs(*s) <= 1) << "s = " << s;
  }

  void update_q(const Float& c, const Float& s, int k) {
    CHECK(0 <= k);
    CHECK(typename eigen_mat::Index(k + 1) < _q.cols());

    Float tmp;
    for (typename eigen_mat::Index i = 0; i < _q.cols(); ++i) {
      tmp = c * _q(i, k) - s * _q(i, k + 1);
      _q(i, k + 1) = s * _q(i, k) + c * _q(i, k + 1);
      _q(i, k) = tmp;
    }
  }

  void apply_half_givens(Float const& c, Float const& s, Float* col,
                         int a_idx) {
    CHECK(col);
    CHECK(0 <= a_idx);
    CHECK(a_idx < num_rows() - 1);

    col[a_idx] = c * col[a_idx] - s * col[a_idx + 1];
  }

  void apply_givens(Float const& c, Float const& s, int r_col_begin,
                    int a_idx) {
    CHECK(0 <= r_col_begin);
    CHECK(r_col_begin <= num_cols());
    CHECK(0 <= a_idx);
    CHECK(a_idx < num_rows() - 1);

    Float tmp;
    for (int j = r_col_begin; j < num_cols(); ++j) {
      tmp = c * _r_begin[j][a_idx] - s * _r_begin[j][a_idx + 1];
      _r_begin[j][a_idx + 1] =
          s * _r_begin[j][a_idx] + c * _r_begin[j][a_idx + 1];
      _r_begin[j][a_idx] = tmp;
    };
  }

  /**
    @param col_end the index of the column past the last column to update
  */
  void hessenberg_update(int k, int col_end) {
    Float c, s;
    for (int i = k; i < col_end; ++i) {
      givens(_r_begin[i][i], _r_begin[i][i + 1], &c, &s);
      apply_half_givens(c, s, _r_begin[i], i);
      apply_givens(c, s, i + 1, i);
      update_q(c, s, i);
    }
  }

  eigen_mat _q;
  Float* _r_raw;
  Float** _r_begin;
  Float** _r_end;
};

}  // namespace dynqr

#endif  // DYNQR_H_
