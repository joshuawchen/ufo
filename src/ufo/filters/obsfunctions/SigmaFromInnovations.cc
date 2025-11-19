#include "ufo/filters/obsfunctions/SigmaFromInnovations.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

// Register with the ObsFunction factory
static ObsFunctionMaker<SigmaFromInnovations> makerSigmaFromInnovations_(
    "SigmaFromInnovations");

// -----------------------------------------------------------------------------

SigmaFromInnovations::SigmaFromInnovations(const eckit::LocalConfiguration & conf)
  : obsGroup_("ObsValue"),
    hofxGroup_("HofX"),
    eMode_(0.0),
    invars_()
{
  oops::Log::trace() << "SigmaFromInnovations constructor" << std::endl;
  oops::Log::debug() << "SigmaFromInnovations config = " << conf << std::endl;

  // Optional groups
  conf.get("obs group",  obsGroup_);
  conf.get("hofx group", hofxGroup_);

  // Required scalar: innovation_mode
  conf.get("innovation_mode", eMode_);

  // Required arrays: innovation_grid and variance_table
  conf.get("innovation_grid", grid_);
  conf.get("variance_table", variance_);

  if (grid_.empty() || variance_.empty() || grid_.size() != variance_.size()) {
    throw eckit::UserError(
      "SigmaFromInnovations: innovation_grid and variance_table must be non-empty "
      "and of equal length", Here());
  }

  // Ensure innovation_grid is strictly increasing
  for (size_t i = 1; i < grid_.size(); ++i) {
    if (!(grid_[i] > grid_[i - 1])) {
      std::ostringstream oss;
      oss << "SigmaFromInnovations: innovation_grid must be strictly increasing; "
          << "grid[" << i-1 << "]=" << grid_[i-1]
          << ", grid[" << i   << "]=" << grid_[i];
      throw eckit::UserError(oss.str(), Here());
    }
  }

  // Ensure variances are non-negative
  for (double & v : variance_) {
    if (v < 0.0) {
      throw eckit::UserError(
        "SigmaFromInnovations: variance_table entries must be >= 0", Here());
    }
  }

  // We deliberately leave invars_ empty: we will pull ObsValue/HofX on demand
}

// -----------------------------------------------------------------------------

double SigmaFromInnovations::varianceFromCenteredInnovation(double r) const {
  const size_t n = grid_.size();
  if (n == 1) {
    return variance_[0];
  }

  // Clamp to endpoints (no extrapolation)
  if (r <= grid_.front()) {
    return variance_.front();
  }
  if (r >= grid_.back()) {
    return variance_.back();
  }

  // Find first grid point > r
  auto upper = std::upper_bound(grid_.begin(), grid_.end(), r);
  const size_t j = static_cast<size_t>(upper - grid_.begin());
  const size_t i = j - 1;

  const double x0 = grid_[i];
  const double x1 = grid_[j];
  const double y0 = variance_[i];
  const double y1 = variance_[j];

  const double t = (r - x0) / (x1 - x0);
  return (1.0 - t) * y0 + t * y1;
}

// -----------------------------------------------------------------------------

void SigmaFromInnovations::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "SigmaFromInnovations compute start" << std::endl;

  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  // We expect exactly one output variable for this ObsFunction instance
  ASSERT(out.nvars() == 1);
  const std::string varName = out.varnames()[0];  // e.g. "airTemperature"
