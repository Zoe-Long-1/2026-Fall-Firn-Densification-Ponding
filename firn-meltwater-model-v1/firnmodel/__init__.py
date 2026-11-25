"""firnmodel: 1D firn column with heat conduction, densification and meltwater.

Standalone simulation package used by the web app (through Pyodide), the
verification notebook and the tests.

    import firnmodel
    results = firnmodel.run({"T_mean": -12.0, "years": 20})
"""

from .model import run, to_jsonable
from .params import DEFAULTS, make_params

__all__ = ["run", "to_jsonable", "DEFAULTS", "make_params"]
__version__ = "1.0.0"
