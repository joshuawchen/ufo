#ifndef UFO_FILTERS_OBSFUNCTIONS_SIGMAFROMINNOVATIONS_H_
#define UFO_FILTERS_OBSFUNCTIONS_SIGMAFROMINNOVATIONS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// Options for SigmaFromInnovations ObsFunction.
///
/// YAML example:
///   test variables:
///   - name: ObsFunction/SigmaFromInnovations
///     options:
///       innovation_mode: 0.0
///       innovation_grid:   [-5.0, -4.0, ..., 4.0, 5.0]
///       variance_table:    [ ... same length as innovation_grid ... ]
///       obs group:         "ObsValue"       # optional
///       hofx group:        "HofX"           # optional
///       variable:          "airTemperature" # optional; see note below
class SigmaFromInnovationsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SigmaFromInnovationsParameters, Parameters)

 public:
  /// Mode of the true noise in innovation units (e_mode).
  oops::Parameter<double> innovationMode{
    "innovation_mode", 0.0, this};

  /// Grid values r_i (will be used on centered innovations r = e - e_mode).
  oops::Parameter<std::vector<double>> innovationGrid{
    "innovation_grid", {}, this};

  /// sigma^2(r_i) values, same length as innovationGrid.
  oops::Parameter<std::vector<double>> varianceTable{
    "variance_table", {}, this};

  /// Group name for observations (usually "ObsValue").
  oops::Parameter<std::string> obsGroup{
    "obs group", "ObsValue", this};

  /// Group name for HofX (background/analysis equivalent).
  oops::Parameter<std::string> hofxGroup{
    "hofx group", "HofX", this};

  /// Optional name of the observed variable this function targets.
  ///
  /// If set, requiredVariables() declares <obs group>/<variable> and
  /// <hofx group>/<variable>, so the framework knows this function depends on
  /// H(x) and schedules it after the obs operator (the post stage) on its own.
  /// If left unset, the target variable is inferred from the output variable at
  /// compute() time, in which case the filter must use 'defer to post: true'.
  oops::OptionalParameter<std::string> variable{
    "variable", this};
};

/// \brief ObsFunction that returns a Gaussian sigma as a function of innovation.
///
/// For each observation j:
///   e_j = y_j - Hx_j
///   r_j = e_j - innovation_mode
///   sigma^2(r_j) is obtained by clamped linear interpolation in variance_table
///   output_j = sqrt( sigma^2(r_j) )
///
/// Outside the grid range we **clamp** to the nearest end (no extrapolation).
class SigmaFromInnovations : public ObsFunctionBase<float> {
 public:
  static const std::string classname() { return "SigmaFromInnovations"; }

  explicit SigmaFromInnovations(const eckit::LocalConfiguration &);
  ~SigmaFromInnovations() override = default;

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const override;

  const ufo::Variables & requiredVariables() const override;

 private:
  // Parsed options
  SigmaFromInnovationsParameters options_;

  // Convenience copies
  std::string obsGroup_;
  std::string hofxGroup_;
  double eMode_;
  std::vector<double> grid_;
  std::vector<double> variance_;

  // Explicit target variable (empty if not provided; then inferred from output).
  std::string variable_;

  // Required variables. Populated with <obs group>/<var> and <hofx group>/<var>
  // only when the optional 'variable' is set, which lets the framework schedule
  // this function after H(x) automatically. Otherwise left empty and the
  // function relies on 'defer to post: true' for correct scheduling.
  ufo::Variables requiredVars_;

  /// Interpolate sigma^2(r) with clamping at the endpoints.
  double varianceFromCenteredInnovation(double r) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SIGMAFROMINNOVATIONS_H_
