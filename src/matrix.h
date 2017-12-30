#ifndef MATRIX_H
#define MATRIX_H

#include <memory>

namespace dynqr {
namespace internal {
// @brief A square matrix which only supports A' * v, and basic element access
template <typename Float_>
class Matrix {
 public:
  typedef Float_ Float;

 private:
  std::unique_ptr<Float[]> data_;
};
}  // namespace internal
}  // namespace dynqr

#endif  // MATRIX_H
