# SPDX-License-Identifier: Apache-2.0
"""edmore.demography — Python companion for the Fortran demographic engine.

Loads ``libedmore.so`` (built with ``-DED_BUILD_PYLIB=ON``) and drives the
standalone carbon slow loop through the C-API, reading cohort/patch state back
into numpy. Importing this submodule loads the shared library; ``import edmore``
alone does not.
"""
from ._site import Config, Site

__all__ = ["Config", "Site"]
