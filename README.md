# gflow1
GFLOW is the original analytic element groundwater solver written by Henk M.
Haitjema. GFLOW is an analytic element model, which solves steady state
groundwater flow in a single aquifer. GFLOW supports three-dimensional particle
tracking, but employs the Dupuit-Forchheimer approximation, thereby ignoring
resistance to vertical flow.

## Status
The initial commit of this repository consists of the Fortran code associated
with the public release version gflow1.exe (1/24/2018), which is the
computional core of GFLOW. The gflow1 code was originally only compiled with
the Lahey Fortran compiler on Windows. This compiler is discontinued as Lahey
Computer Systems has been permanently closed since December 31, 2022

The code has been successfully migrated to gfortran and supports multiple
platforms. It can now be compiled not just on Windows, but also on macOS and
Linux. Executables are available for Windows, macOS, and Linux under the
releases section.

## Goals
This repository serves two main goals:

-	To enable use of GFLOW models on Windows, macOS, and Linux now and in the
future.
- To publicly host the gflow1 code for archival purposes.

## Roadmap
We have the following (tentative) plans, primarily focused on improving
robustness and usability:

-	Support for mixed case paths: gflow1 is not case sensitive when dealing with
file paths and will internally transform file paths to upper case, since it was
a Windows only application. While Windows is not case sensitive for file paths,
macOS and Linux are. Ideally, the casing of paths can be preserved.
-	Support for long paths: gflow1 currently supports paths of no longer than 256
characters.
-	Fix compiler warnings: compilation with gflow1 currently emits many warnings.
Ideally, we address all of these.
-	Multi-compiler support: we currently only compile with gfortran. We intend to
add Intel Fortran as well (and address new warnings as they pop up).
-	The ` '-freal-4-real-8'` flag is required to force all reals to double
precision, since the executable otherwise crashes when run. Ideally, all type
declarations are updated to explicit double precision.
-	Add a basic set of integration tests.

## QGIS-plugin
A rudimentary QGIS plugin is available at
https://github.com/huite/gflow-plugin.

While it doesn't support all features of the GFLOW GUI, QGIS offers powerful
complementary capabilities including:

- Comprehensive geospatial data interaction
- Advanced visualization and data presentation
- Flexible vector geometry editing
- Access to various basemaps
