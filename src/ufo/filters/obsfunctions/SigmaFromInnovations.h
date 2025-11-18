#pragma once

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "ufo/filters/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief ObsFunction that returns a Gaussian sigma as a function of innovation.
///
///  options:
///    innovation_mode:  <float>  # e_mode, the mode of the *true* noise in innovation units
///    innovation_grid:  [ ... ]  # 1D grid values r_i (centered innovations = e - e_mode)
///    variance_table:   [ ... ]  # same length, sigma^2(r_i)
///    obs group:  "ObsValue"     # (optional) group for y (defaults shown)
///    hofx group: "HofX"         # (optional) group for H(x)
///
/// For each obs & variable:
///   e = y - Hx
///   r = e - innovation_mode
///   sigma^2(r) is obtained by clamped linear interpolation in variance_table
///   output = sqrt( sigma^2(r) )
///
/// Outside the grid range we **clamp** to the nearest end (no extrapolation).
class SigmaFromInnovations : public ObsFunctionBase<float> {
 public:
  explicit SigmaFromInnovations(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const override;

  /// We don’t require pre-fetched variables explicitly; we will read
  /// ObsValue/HofX on demand.
  const ufo::Variables & requiredVariables() const override;

 private:
  // Options
  std::string obsGroup_;   ///< Group name for observations (default "ObsValue")
  std::string hofxGroup_;  ///< Group name for HofX (default "HofX")
  double eMode_;           ///< Innovation mode (location of true noise mode)

  std::vector<double> grid_;      ///< innovation_grid (centered innovations)
  std::vector<double> variance_;  ///< variance_table (sigma^2 at each grid point)

  ufo::Variables requiredVars_;   ///< left empty

  /// Interpolate sigma^2(r) with clamping at the endpoints.
  double varianceFromCenteredInnovation(double r) const;
};

}  // namespace ufo
