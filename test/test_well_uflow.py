import contextlib
import math
import pathlib
import os
import subprocess

import pytest


@contextlib.contextmanager
def cd(path: pathlib.Path):
    """
    Change directory, and change it back after the with block.

    Examples
    --------
    >>> with imod.util.context.cd("docs"):
            do_something_in_docs()

    """
    curdir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(curdir)
        

class TestRun:
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = pathlib.Path(tmpdir)

        # Get the root path and add the bin folder to PATH so we can find the
        # recently built executable.
        project_root = pathlib.Path(__file__).parent.parent.resolve()
        os.environ["PATH"] += os.pathsep + str(project_root / "bin")

        # Copy the file
        source_file = project_root / 'test/data/well_uflow.dat'
        target_file = self.tmpdir / 'WELL_UFLOW.DAT'
        target_file.write_bytes(source_file.read_bytes())

        # Run gflow
        with cd(self.tmpdir):
            self.result = subprocess.run(
                ["gflow2", "WELL_UFLOW.DAT"],
                timeout=2,  # two seconds
                check=True,
                capture_output=True,
                text=True
            )
            
    def check_exitcode(self):
        assert self.result.returncode == 0
 
    def test_files(self):
        assert (self.tmpdir / "WELL_UFLOW.GRD").exists()
        assert (self.tmpdir / "WELL_UFLOW.XTR").exists()

    def test_heads_match_analytic_solution(self):
        """
        Compare the gridded heads with the analytic solution for a well in a
        uniform flow field, Haitjema (1995), Analytic Element Modeling of
        Groundwater Flow, section 3.1.10, equation (3.110):

            Phi = -Q0 * x + Q / (2 pi) * ln(r) + C

        For confined flow the head is Phi / (k H) plus a constant. The constant
        follows from the reference point in well_uflow.dat, which is located
        at the well (0, 0) with head 100; inside the well radius the potential
        is evaluated at the well radius.
        """
        # Parameters from test/data/well_uflow.dat. Units are feet and days;
        # the window coordinates (x, y) are in ft as well.
        k = 10.0  # hydraulic conductivity, perm [ft/day]
        H = 10.0  # aquifer thickness, thick [ft]
        Q0 = 1.0  # uniform flow rate in the x-direction, uniflow [ft^2/day]
        Q = 50.0  # well discharge [ft^3/day]
        rw = 1.0  # well radius [ft]
        href = 100.0  # head at the reference point (0, 0) [ft]
        T = k * H  # transmissivity [ft^2/day]

        # Read the Surfer ASCII grid (DSAA). Values are stored row by row,
        # starting at ymin, with nx values per row.
        tokens = (self.tmpdir / "WELL_UFLOW.GRD").read_text().split()
        assert tokens[0] == "DSAA"
        nx, ny = int(tokens[1]), int(tokens[2])
        xmin, xmax, ymin, ymax = (float(v) for v in tokens[3:7])
        heads = [float(v) for v in tokens[9:]]
        assert len(heads) == nx * ny

        dx = (xmax - xmin) / (nx - 1)
        dy = (ymax - ymin) / (ny - 1)
        max_error = 0.0
        for row in range(ny):
            y = ymin + row * dy
            for col in range(nx):
                x = xmin + col * dx
                r = max(math.hypot(x, y), rw)
                expected = (-Q0 * x + Q / (2.0 * math.pi) * math.log(r / rw)) / T + href
                max_error = max(max_error, abs(heads[row * nx + col] - expected))

        # GFLOW writes heads with 4 decimals, so allow for rounding [ft].
        assert max_error < 1.0e-3
