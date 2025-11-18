#include "ufo/filters/obsfunctions/SigmaFromInnovations.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionMaker.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Factory registration
// -----------------------------------------------------------------------------

static ObsFunctionMaker<SigmaFromInnovations>
  makerSigmaFromInnovations("SigmaFromInnovations");

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------

SigmaFromInnovations::SigmaFromInnovations(const eckit::LocalConfiguration & conf)
  : obsGroup_(conf.getString("obs group", "ObsValue")),
    hofxGroup_(conf.getString("hofx group", "HofX")),
    eMode_(conf.getDouble("innovation_mode")),
    grid_(conf.getDoubleVector("innovation_grid")),
    variance_(conf.getDoubleVector("variance_table")),
    requiredVars_()
{
  if (grid_.size() < 2 || grid_.size() != variance_.size()) {
    throw eckit::UserError(
      "SigmaFromInnovations: 'innovation_grid' and 'variance_table' must have the "
      "same length >= 2", Here());
  }

  // Ensure strictly increasing grid (user responsibility to provide it this way).
  for (std::size_t i = 1; i < grid_.size(); ++i) {
    if (!(grid_[i] > grid_[i-1])) {
      throw eckit::UserError(
        "SigmaFromInnovations: 'innovation_grid' must be strictly increasing.", Here());
    }
  }

  // (Optional) sanity check: variances should be non-negative.
  for (std::size_t i = 0; i < variance_.size(); ++i) {
    if (variance_[i] < 0.0) {
      throw eckit::UserError(
        "SigmaFromInnovations: 'variance_table' entries must be >= 0.", Here());
    }
  }

  oops::Log::info() << "SigmaFromInnovations constructed: obs group = "
                    << obsGroup_ << ", hofx group = " << hofxGroup_
                    << ", mode = " << eMode_
                    << ", ngrid = " << grid_.size() << std::endl;
}

// -----------------------------------------------------------------------------
// requiredVariables
// -----------------------------------------------------------------------------

const ufo::Variables & SigmaFromInnovations::requiredVariables() const {
  // We read ObsValue/HofX directly from ObsSpace using the variables in out.varnames().
  // So we leave this empty.
  return requiredVars_;
}

// -----------------------------------------------------------------------------
// Helper: piecewise-linear interpolation of variance with endpoint clamping
// -----------------------------------------------------------------------------

double SigmaFromInnovations::varianceFromCenteredInnovation(double r) const {
  const std::size_t n = grid_.size();

  // Clamp at ends (do NOT extrapolate)
  if (r <= grid_.front()) {
    return variance_.front();
  }
  if (r >= grid_.back()) {
    return variance_.back();
  }

  // Find first grid point > r, then interpolate between that and previous.
  auto it = std::upper_bound(grid_.begin(), grid_.end(), r);
  std::size_t idx1 = static_cast<std::size_t>(std::distance(grid_.begin(), it));
  std::size_t idx0 = idx1 - 1;

  const double x0 = grid_[idx0];
  const double x1 = grid_[idx1];
  const double v0 = variance_[idx0];
  const double v1 = variance_[idx1];

  const double t = (r - x0) / (x1 - x0);  // in (0, 1)
  return v0 + t * (v1 - v0);
}

// -----------------------------------------------------------------------------
// compute
// -----------------------------------------------------------------------------

void SigmaFromInnovations::compute(const ObsFilterData & data,
                                   ioda::ObsDataVector<float> & out) const {
  const std::size_t nlocs = data.nlocs();
  const std::size_t nvars = out.nvars();
  const std::vector<std::string> & varNames = out.varnames();

  // Read y and H(x) for these variables from the requested groups.
  ioda::ObsDataVector<float> obs(data.obsdb(), varNames, obsGroup_);
  ioda::ObsDataVector<float> hofx(data.obsdb(), varNames, hofxGroup_);

  const float missing = util::missingValue<float>();

  for (std::size_t jv = 0; jv < nvars; ++jv) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      const float y   = obs[jv][jl];
      const float hx  = hofx[jv][jl];

      // If either is missing, propagate missing.
      if (!std::isfinite(static_cast<double>(y)) ||
          !std::isfinite(static_cast<double>(hx))) {
        out[jv][jl] = missing;
        continue;
      }

      // Innovation and centering by true noise mode.
      const double e = static_cast<double>(y - hx);
      const double r = e - eMode_;

      double var = varianceFromCenteredInnovation(r);
      if (var <= 0.0) {
        // Avoid NaN from sqrt; enforce small positive floor.
        var = std::numeric_limits<double>::min();
      }

      const double sigma = std::sqrt(var);
      out[jv][jl] = static_cast<float>(sigma);
    }
  }
}

}  // namespace ufo
