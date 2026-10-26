# SPDX-License-Identifier: Apache-2.0
"""netCDF writer for the demography example's runs.

The demography example drives ED-MORE's law-free apply-primitives from Python with its own
vital-rate laws, so nothing in the Fortran output layer is involved: this module writes the
stand itself, in the ragged cohort/patch layout the plotting scripts beside it read.

**Why the example writes its own file at all.** The laws live in the example, not in Fortran, so
`edmore` cannot produce this run. The figures need a file, and this is the only thing that can
make one.

The schema is deliberately the one the existing plotters expect: one record per sampled year,
a `cohort` axis padded to a fixed capacity, a `patch` axis likewise, the CSR map
(`cohort_offset` / `cohort_count`) that says which cohorts belong to which patch, and the site
totals. Unused slots are left unwritten so netCDF fills them, which is what `_FillValue` is for
and what a reader masking on it expects (the same convention the Fortran output layer uses for
its own padding).
"""
from __future__ import annotations

import numpy as np
from netCDF4 import Dataset

#----- Per-cohort fields copied straight out of the engine, plus `basal_area`, which is derived
#      from dbh here rather than carried: it is pi*(dbh/2)^2 by definition and there is no reason
#      for the C-API to export a second copy of the same number.
_COHORT_REAL = ("dbh", "height", "nplant", "agb", "leaf_area", "growth_avg")
_COHORT_INT = ("pft", "owner_patch")
_PATCH_REAL = ("area", "age")
_PATCH_INT = ("dist_type", "cohort_offset", "cohort_count")


class StandWriter:
    """Accumulates yearly stand snapshots, then writes them as one ragged netCDF."""

    def __init__(self, path, cohort_cap=2048, patch_cap=64):
        self.path = str(path)
        self.cohort_cap = int(cohort_cap)
        self.patch_cap = int(patch_cap)
        self.records = []

    def sample(self, site, year):
        """Copy the whole stand out of the engine for one record."""
        r = {"year": int(year), "n_cohort": int(site.n_cohort), "n_patch": int(site.n_patch),
             "total_agb": float(site.total_agb), "total_lai": float(site.total_lai),
             "total_nplant": float(site.total_nplant),
             "total_basal_area": float(site.total_basal_area)}
        for f in _COHORT_REAL:
            r[f] = site.get(f)
        for f in _COHORT_INT:
            r[f] = site.get(f)
        r["global_cohort_id"] = site.get("global_id")
        for f in _PATCH_REAL:
            r["patch_" + f] = site.get_patch(f)
        for f in _PATCH_INT:
            r[f if f != "dist_type" else "dist_type"] = site.get_patch(f)
        r["global_patch_id"] = site.get_patch("global_id")
        #----- pi*(dbh/2)^2 in cm2, matching the Fortran `basal_area` column. -------------------
        r["basal_area"] = np.pi * (r["dbh"] / 2.0) ** 2
        #----- mean dbh, nplant-weighted -- the one site scalar not exposed as a total. ---------
        w = r["nplant"]
        r["mean_dbh"] = float((r["dbh"] * w).sum() / w.sum()) if w.size and w.sum() > 0 else 0.0
        self.records.append(r)

    def write(self):
        nt = len(self.records)
        if nt == 0:
            raise RuntimeError("StandWriter.write: nothing sampled")
        need_c = max(r["n_cohort"] for r in self.records)
        need_p = max(r["n_patch"] for r in self.records)
        if need_c > self.cohort_cap or need_p > self.patch_cap:
            raise RuntimeError(
                f"stand outgrew the writer's capacity: {need_c} cohorts / {need_p} patches "
                f"against {self.cohort_cap} / {self.patch_cap}. Raise the caps rather than "
                f"truncating -- a silently clipped stand is worse than a failed write.")

        with Dataset(self.path, "w", format="NETCDF4") as d:
            d.title = "ED-MORE demography with census-trained vital rates (examples/example03_demography)"
            d.Conventions = "CF-1.8"
            d.createDimension("time", None)
            d.createDimension("cohort", self.cohort_cap)
            d.createDimension("patch", self.patch_cap)

            def var(name, dims, kind="f8", units="", long_name=""):
                v = d.createVariable(name, kind, dims, zlib=True, complevel=1,
                                     fill_value=(-2147483647 if kind == "i4" else 9.969209968386869e+36))
                v.units = units
                v.long_name = long_name or name
                return v

            var("time", ("time",), units="year", long_name="decimal calendar year (period start)")
            var("year", ("time",), "i4", "1", "calendar year")
            var("month", ("time",), "i4", "1", "calendar month")
            var("day", ("time",), "i4", "1", "calendar day")
            var("n_cohort", ("time",), "i4", "-", "live cohorts this record")
            var("n_patch", ("time",), "i4", "-", "live patches this record")

            spec_c = {"dbh": ("cm", "cohort mean diameter at breast height"),
                      "height": ("m", "cohort height"),
                      "nplant": ("plant/m2", "stem density (per m2 of its own patch)"),
                      "agb": ("kgC/plant", "aboveground biomass per plant"),
                      "leaf_area": ("m2/plant", "leaf area per plant"),
                      "growth_avg": ("cm/yr", "moving-average diameter growth"),
                      "basal_area": ("cm2/plant", "basal area per plant")}
            for k, (u, ln) in spec_c.items():
                var(k, ("time", "cohort"), "f8", u, ln)
            var("pft", ("time", "cohort"), "i4", "-", "plant functional type")
            var("owner_patch", ("time", "cohort"), "i4", "-", "1-based owning patch index")
            var("global_cohort_id", ("time", "cohort"), "i4", "-",
                "persistent cohort id, never reused")

            var("patch_area", ("time", "patch"), "f8", "-", "patch area fraction of the site")
            var("patch_age", ("time", "patch"), "f8", "yr", "time since last disturbance")
            var("dist_type", ("time", "patch"), "i4", "-", "disturbance class")
            var("cohort_offset", ("time", "patch"), "i4", "-", "1-based first cohort of the patch")
            var("cohort_count", ("time", "patch"), "i4", "-", "cohorts in the patch")
            var("global_patch_id", ("time", "patch"), "i4", "-", "persistent patch id, never reused")

            var("total_agb", ("time",), "f8", "kgC/m2", "site aboveground biomass")
            var("total_lai", ("time",), "f8", "m2/m2", "site leaf area index")
            var("total_nplant", ("time",), "f8", "plant/m2", "site stem density")
            var("total_basal_area", ("time",), "f8", "m2/m2", "site basal area")
            var("mean_dbh", ("time",), "f8", "cm", "stem-weighted mean diameter")

            for i, r in enumerate(self.records):
                nc_, np_ = r["n_cohort"], r["n_patch"]
                d["time"][i] = float(r["year"])
                d["year"][i] = r["year"]
                d["month"][i] = 1
                d["day"][i] = 1
                d["n_cohort"][i] = nc_
                d["n_patch"][i] = np_
                for k in list(spec_c) + ["pft", "owner_patch"]:
                    if nc_:
                        d[k][i, :nc_] = r[k][:nc_]
                if nc_:
                    d["global_cohort_id"][i, :nc_] = r["global_cohort_id"][:nc_]
                if np_:
                    d["patch_area"][i, :np_] = r["patch_area"][:np_]
                    d["patch_age"][i, :np_] = r["patch_age"][:np_]
                    d["dist_type"][i, :np_] = r["dist_type"][:np_]
                    d["cohort_offset"][i, :np_] = r["cohort_offset"][:np_]
                    d["cohort_count"][i, :np_] = r["cohort_count"][:np_]
                    d["global_patch_id"][i, :np_] = r["global_patch_id"][:np_]
                for k in ("total_agb", "total_lai", "total_nplant", "total_basal_area", "mean_dbh"):
                    d[k][i] = r[k]
        return self.path
