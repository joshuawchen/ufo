#include "ufo/filters/obsfunctions/SigmaFromInnovations.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

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
  : options_(),
    obsGroup_("ObsValue"),
    hofxGroup_("HofX"),
    eMode_(0.0),
    grid_(),
    variance_(),
    variable_(),
    requiredVars_()
{
  oops::Log::trace() << "SigmaFromInnovations constructor" << std::endl;
  oops::Log::debug() << "SigmaFromInnovations config = " << conf << std::endl;

  // Parse YAML options into options_
  options_.deserialize(conf);

  // Copy convenient values
  obsGroup_ = options_.obsGroup.value();
  hofxGroup_ = options_.hofxGroup.value();
  eMode_     = options_.innovationMode.value();
  grid_      = options_.innovationGrid.value();
  variance_  = options_.varianceTable.value();

  // Basic checks
  if (grid_.empty() || variance_.empty() || grid_.size() != variance_.size()) {
    throw eckit::UserError(
      "SigmaFromInnovations: innovation_grid and variance_table must be non-empty "
      "and have the same length.",
      Here());
  }

  // innovation_grid must be strictly increasing
  for (size_t i = 1; i < grid_.size(); ++i) {
    if (!(grid_[i] > grid_[i - 1])) {
      std::ostringstream oss;
      oss << "SigmaFromInnovations: innovation_grid must be strictly increasing; "
          << "grid[" << i-1 << "]=" << grid_[i-1]
          << ", grid[" << i   << "]=" << grid_[i];
      throw eckit::UserError(oss.str(), Here());
    }
  }

  // variances must be non-negative
  for (double & v : variance_) {
    if (v < 0.0) {
      throw eckit::UserError(
        "SigmaFromInnovations: variance_table entries must be >= 0.", Here());
    }
  }

  // Optional explicit target variable. When provided, declare <obs group>/<var>
  // and <hofx group>/<var> as required inputs, so the framework knows this
  // function depends on H(x) and schedules it in the post stage automatically.
  // When omitted, requiredVars_ stays empty and the target variable is taken
  // from the output variable in compute() ('defer to post: true' is then needed).
  if (options_.variable.value() != boost::none) {
    variable_ = *options_.variable.value();
    requiredVars_ += Variable(obsGroup_  + "/" + variable_);
    requiredVars_ += Variable(hofxGroup_ + "/" + variable_);
  }
}

// -----------------------------------------------------------------------------

double SigmaFromInnovations::varianceFromCenteredInnovation(double r) const {
  // clamp + linear interpolation in sigma^2, as requested

  const size_t n = grid_.size();
  if (n == 1) {
    return variance_[0];
  }

  // Clamp to endpoints: use min/max sigma^2 outside grid range
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
  const double y0 = variance_[i];  // sigma^2 at x0
  const double y1 = variance_[j];  // sigma^2 at x1

  const double t = (r - x0) / (x1 - x0);  // 0 <= t <= 1
  return (1.0 - t) * y0 + t * y1;
}

// -----------------------------------------------------------------------------

void SigmaFromInnovations::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "SigmaFromInnovations compute start" << std::endl;

  const size_t nlocs  = in.nlocs();
  const float missing = util::missingValue<float>();

  // One output variable per ObsFunction instance
  ASSERT(out.nvars() == 1);
  // Target variable: explicit 'variable' option if given, else the output slot.
  const std::string varName = variable_.empty() ? out.varnames()[0] : variable_;

  // Read y and H(x) for this variable
  std::vector<float> y, hofx;
  in.get(Variable(obsGroup_ + "/"  + varName), y);
  in.get(Variable(hofxGroup_ + "/" + varName), hofx);

  if (y.size() != nlocs || hofx.size() != nlocs) {
    throw eckit::UserError(
      "SigmaFromInnovations: size mismatch between ObsValue/HofX and nlocs.",
      Here());
  }

  // Compute sigma(e) for each location
  for (size_t j = 0; j < nlocs; ++j) {
    if (y[j] == missing || hofx[j] == missing) {
      out[0][j] = missing;
      continue;
    }

    const double e = static_cast<double>(y[j]) - static_cast<double>(hofx[j]);  // innovation
    const double r = e - eMode_;  // centered innovation, relative to noise mode

    double var = varianceFromCenteredInnovation(r);
    if (var < 0.0) {
      // Guard against tiny negatives from interpolation / roundoff
      var = 0.0;
    }

    out[0][j] = (var > 0.0) ? static_cast<float>(std::sqrt(var))
                            : 0.0f;
  }

  oops::Log::trace() << "SigmaFromInnovations compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & SigmaFromInnovations::requiredVariables() const {
  // Empty unless the optional 'variable' was set (see constructor). When empty,
  // use 'defer to post: true' in the filter to guarantee current-iterate H(x).
  return requiredVars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
