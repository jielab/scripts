// Exercise chr1-tail indices without reading genotypes or touching full buffers.
#include "Simd.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>

int main()
{
  const Eigen::Index sites = 2054102, states = 69, batch = 16;
  const Eigen::Index count = sites * states * batch;
  if (count <= std::numeric_limits<int>::max()) return 1;
  // Eigen leaves this virtual address space uninitialised. Only the final site's
  // pages are touched; the test therefore needs very little physical memory.
  Eigen::ArrayXf alpha(count), beta(count), scale(sites * batch);
  const Eigen::Index begin = (sites - 1) * states * batch;
  alpha.segment(begin, states * batch).setOnes();
  beta.segment(begin, states * batch).setOnes();
  asmc::normalizeAlphaWithBeta(alpha, beta, scale, batch, states, sites - 1, sites);
  for (Eigen::Index index = begin; index < count; ++index) {
    if (!std::isfinite(alpha[index]) || std::abs(alpha[index] - 1.0f / states) > 1e-7f)
      return 2;
  }
  std::cout << "PASS: SIMD posterior indices beyond INT32_MAX at chr1 tail\n";
}
