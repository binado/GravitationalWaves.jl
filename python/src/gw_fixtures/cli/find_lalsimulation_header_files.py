"""Print the directory holding liblal / liblalsimulation for installed lalsuite.

Keep in sync with PYTHON_LIBDIR_SNIPPET in packages/LALWaveforms/src/library.jl.
"""

import glob
import os
import sys

import lal
import lalsimulation  # noqa: F401


def main():
    site = os.path.dirname(os.path.dirname(lal.__file__))
    for d in (
        os.path.join(site, "lal", ".dylibs"),  # macOS wheel
        os.path.join(site, "lalsuite.libs"),  # Linux (auditwheel) wheel
        os.path.join(sys.prefix, "lib"),  # conda-style install
    ):
        if glob.glob(os.path.join(d, "liblalsimulation*")):
            print(d)
            return
    sys.exit("could not find liblalsimulation next to the lal Python package")


if __name__ == "__main__":
    main()
