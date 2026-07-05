#!/usr/bin/env python3
# (C) Copyright 2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Generates the obs file + YAML for the SigmaFromInnovations ObsFunction ctest.
# Run at build time by CMake. Deterministic (fixed seed) so the test is
# reproducible; regenerates against whatever IODA/HDF5 is present so it can't
# rot. An independent oracle computes the expected sigma -> TestReference.
import os, math
from datetime import datetime, timezone
import numpy as np

OUTDIR = os.environ.get("OUTDIR", ".")
SEED   = int(os.environ.get("SEED", "12345"))
VAR    = "airTemperature"
OUT_NC   = os.path.join(OUTDIR, "sigmafrominnovations_obs.nc4")
OUT_YAML = os.path.join(OUTDIR, "function_sigmafrominnovations.yaml")
rng = np.random.default_rng(SEED)

# ---- config (fixed via seed) ----
ngrid = 5
grid = np.sort(rng.uniform(-4.0, 4.0, ngrid)).astype(float)
for i in range(1, ngrid):
    if grid[i] <= grid[i-1]:
        grid[i] = grid[i-1] + 0.5
grid = grid.tolist()
variance = rng.uniform(0.25, 9.0, ngrid).round(6).tolist()
mode = float(round(rng.uniform(-1.0, 1.0), 6))

# ---- toy obs exercising every branch ----
r_targets = [grid[0]-2.0, grid[0], 0.5*(grid[0]+grid[1]),
             0.5*(grid[-2]+grid[-1]), grid[-1], grid[-1]+3.0]
r_targets = [float(r) for r in r_targets]
N = len(r_targets)
hofx = np.zeros(N)
y = np.array([r + mode for r in r_targets])

def var_from_r(r):
    n = len(grid)
    if n == 1:        return variance[0]
    if r <= grid[0]:  return variance[0]
    if r >= grid[-1]: return variance[-1]
    j = 0
    while j < n and grid[j] <= r:
        j += 1
    i = j - 1
    t = (r - grid[i]) / (grid[j] - grid[i])
    return (1.0 - t) * variance[i] + t * variance[j]

expected = np.array([math.sqrt(max(var_from_r(r), 0.0)) for r in r_targets])

lat = rng.uniform(-89, 89, N).astype("float32")
lon = rng.uniform(0, 359, N).astype("float32")
dt_val = int(datetime(2020, 1, 1, 3, 0, 0, tzinfo=timezone.utc).timestamp())
dtime = np.full(N, dt_val, dtype="int64")
EPOCH_UNITS = "seconds since 1970-01-01T00:00:00Z"

def write_netcdf4():
    from netCDF4 import Dataset
    ds = Dataset(OUT_NC, "w", format="NETCDF4")
    ds.createDimension("Location", N)
    ds.createVariable("Location", "i4", ("Location",))[:] = np.arange(N, dtype="int32")
    def gvar(path, data, dtype, units=None):
        grp, name = path.rsplit("/", 1)
        g = ds.createGroup(grp) if grp not in ds.groups else ds.groups[grp]
        v = g.createVariable(name, dtype, ("Location",)); v[:] = data
        if units: v.units = units
    gvar("MetaData/latitude",  lat,   "f4", "degrees_north")
    gvar("MetaData/longitude", lon,   "f4", "degrees_east")
    gvar("MetaData/dateTime",  dtime, "i8", EPOCH_UNITS)
    gvar("ObsValue/" + VAR,      y.astype("float32"),        "f4", "K")
    gvar("HofX/" + VAR,          hofx.astype("float32"),     "f4", "K")
    gvar("TestReference/" + VAR, expected.astype("float32"), "f4", "K")
    ds.setncattr("_ioda_layout", "ObsGroup")
    ds.setncattr("_ioda_layout_version", np.int32(0))
    ds.close()

def write_h5py():
    import h5py
    f = h5py.File(OUT_NC, "w")
    loc = f.create_dataset("Location", data=np.arange(N, dtype="int32")); loc.make_scale("Location")
    def gvar(path, data, units=None):
        d = f.create_dataset(path, data=data); d.dims[0].attach_scale(f["Location"])
        if units: d.attrs["units"] = units
    gvar("MetaData/latitude",  lat,   "degrees_north")
    gvar("MetaData/longitude", lon,   "degrees_east")
    gvar("MetaData/dateTime",  dtime, EPOCH_UNITS)
    gvar("ObsValue/" + VAR,      y.astype("float32"),        "K")
    gvar("HofX/" + VAR,          hofx.astype("float32"),     "K")
    gvar("TestReference/" + VAR, expected.astype("float32"), "K")
    f.attrs["_ioda_layout"] = "ObsGroup"; f.attrs["_ioda_layout_version"] = np.int32(0)
    f.close()

try:
    write_netcdf4()
except Exception:
    write_h5py()

def L(xs): return "[" + ", ".join("%.6g" % x for x in xs) + "]"
yaml = f"""time window:
  begin: 2020-01-01T00:00:00Z
  end:   2020-01-01T06:00:00Z

observations:
- obs space:
    name: SigmaFromInnovationsTest
    obsdatain:
      engine:
        type: H5File
        obsfile: {OUT_NC}
    simulated variables: [{VAR}]
  obs function:
    name: ObsFunction/SigmaFromInnovations
    options:
      obs group: ObsValue
      hofx group: HofX
      innovation_mode: {mode:.6g}
      innovation_grid: {L(grid)}
      variance_table: {L(variance)}
    variables: [{VAR}]
    tolerance: 1.0e-4
"""
open(OUT_YAML, "w").write(yaml)
print("wrote", OUT_NC, "and", OUT_YAML)
