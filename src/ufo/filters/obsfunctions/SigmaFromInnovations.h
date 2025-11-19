#pragma once

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "ufo/filters/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// ObsFunction that returns a Gaussian sigma as a function of innovation.
///
/// YAML example:
///   obs function:
///     name: SigmaFromInnovations
///     options:
///       innovation_mode: 0.0            # e_mode, mode of true noise
///       innovation_grid: [ -4, -2, 0, 2, 4 ]   # grid in *innovation* space
///       variance_table:  [ ... same length ... ]
///       obs group:  ObsValue            # optional
///       hofx group: HofX                # optional
///
/// For each obs & variable:
///   e = y - H(x)
///   r = e - innovation_mode
///   sigma^2(r) obtained by clamped linear interpolation in variance_table
///   output = sqrt( sigma^2(r) )
///
/// Outside the grid range, we **clamp** to the nearest endpoint
/// (no extrapolation, just use min/max variance).
class SigmaFromInnovations : public ObsFunctionBase<float> {
 public:
  explicit SigmaFromInnovations(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const override;

  /// We will re-read ObsValue/HofX ourselves; nothing extra is required.
  const ufo::Variables & requiredVariables() const override;

 private:
  // Options
  std::string obsGroup_;    ///< Group name for observations (default "ObsValue")
  std::string hofxGroup_;   ///< Group name for HofX (default "HofX")
  double eMode_;            ///< innovation_mode (true noise mode)

  // Lookup tables; we store them **centered** internally.
  std::vector<double> grid_;      ///< centered innovation grid r_i = e_i - eMode_
  std::vector<double> variance_;  ///< variance_table: sigma^2(r_i)

  ufo::Variables requiredVars_;   ///< left empty

  /// Linear interpolation in sigma^2 with clamping at endpoints.
  double varianceFromCenteredInnovation(double r) const;
};

}  // namespace ufo
