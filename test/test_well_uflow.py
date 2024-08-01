import contextlib
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
 