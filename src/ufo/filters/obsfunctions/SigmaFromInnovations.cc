#include "ufo/filters/obsfunctions/SigmaFromInnovations.h"

#include <algorithm>
#include <cmath>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/ObsFunctionMaker.h"

namespace ufo {

// Registration: this string is what you put under "name:" in YAML.
static ObsFunctionMaker<SigmaFromInnovations> makerSigmaFromInnovations(
    "SigmaFromInnovations");

// -----------------------------------------------------------------------------
// Constructor: parse options, center & sort the grid
// -----------------------------------------------------------------------------

SigmaFromInnovations::SigmaFromInnovations(const eckit::LocalConfiguration &conf)
  : obsGroup_("ObsValue"), hofxGroup_("HofX"), eMode_(0.0), requiredVars_() {
  // Optional group names
  conf.getIfExists("obs group",  obsGroup_);
  conf.getIfExists("hofx group", hofxGroup_);

  // Required: innovation_mode, innovation_grid, variance_table
  if (!conf.has("innovation_mode")) {
    throw eckit::BadValue("SigmaFromInnovations: 'innovation_mode' option is required", Here());
  }
  conf.get("innovation_mode", eMode_);

  if (!conf.has("innovation_grid") || !conf.has("variance_table")) {
    throw eckit::BadValue("SigmaFromInnovations: 'innovation_grid' and 'variance_table' "
                          "options are required", Here());
  }

  std::vector<double> rawGrid;
  std::vector<double> rawVar;
  conf.get("innovation_grid", rawGrid);
  conf.get("variance_table", rawVar);

  if (rawGrid.empty()) {
    throw eckit::BadValue("SigmaFromInnovations: innovation_grid cannot be empty", Here());
  }
  if (rawGrid.size() != rawVar.size()) {
    throw eckit::BadValue("SigmaFromInnovations: innovation_grid and variance_table "
                          "must have the same length", Here());
  }

  // Center the grid around innovation_mode and sort by r
  std::vector<std::pair<double, double>> pairs;
  pairs.reserve(rawGrid.size());
  for (std::size_t i = 0; i < rawGrid.size(); ++i) {
    const double r = rawGrid[i] - eMode_;  // centered innovation coordinate
    pairs.emplace_back(r, rawVar[i]);
  }

  std::sort(pairs.begin(), pairs.end(),
            [](const std::pair<double, double> &a,
               const std::pair<double, double> &b) {
              return a.first < b.first;
            });

  grid_.resize(pairs.size());
  variance_.resize(pairs.size());
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    grid_[i]     = pairs[i].first;
    variance_[i] = pairs[i].second;
  }
}

// -----------------------------------------------------------------------------
// Interpolation: sigma^2(r) with **clamping** at endpoints
// -----------------------------------------------------------------------------

double SigmaFromInnovations::varianceFromCenteredInnovation(const double r) const {
  // Single point: constant variance.
  if (grid_.size() == 1) {
    return variance_[0];
  }

  // Clamp to ends (no extrapolation).
  if (r <= grid_.front()) {
    return variance_.front();
  }
  if (r >= grid_.back()) {
    return variance_.back();
  }

  // Find first grid point strictly greater than r.
  auto it = std::upper_bound(grid_.begin(), grid_.end(), r);
  const std::size_t idx = static_cast<std::size_t>(std::distance(grid_.begin(), it));

  const double x0 = grid_[idx - 1];
  const double x1 = grid_[idx];
  const double v0 = variance_[idx - 1];
  const double v1 = variance_[idx];

  const double t = (r - x0) / (x1 - x0);
  return (1.0 - t) * v0 + t * v1;
}

// -----------------------------------------------------------------------------
// Compute: sigma(e) for every obs/variable
// -----------------------------------------------------------------------------

void SigmaFromInnovations::compute(const ObsFilterData &data,
                                   ioda::ObsDataVector<float> &out) const {
  const std::size_t nlocs = data.nlocs();
  const std::size_t nvars = out.nvars();

  // Read obs and HofX from requested groups for the *same* variables as "out".
  ioda::ObsDataVector<float> y(data.obsspace(), out.varnames(), obsGroup_);
  ioda::ObsDataVector<float> h(data.obsspace(), out.varnames(), hofxGroup_);

  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      const double e = static_cast<double>(y[jvar][jl]) -
                       static_cast<double>(h[jvar][jl]);
      const double r = e - eMode_;  // centered innovation

      double var = varianceFromCenteredInnovation(r);
      if (var < 0.0) var = 0.0;  // safety
      out[jvar][jl] = static_cast<float>(std::sqrt(var));
    }
  }
}

// -----------------------------------------------------------------------------
// requiredVariables
// -----------------------------------------------------------------------------

const ufo::Variables & SigmaFromInnovations::requiredVariables() const {
  // We are reading ObsValue/HofX directly via ObsDataVector, so nothing special
  // needs to be pre-fetched by the filter framework.
  return requiredVars_;
}

}  // namespace ufo
