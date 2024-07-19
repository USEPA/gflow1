# Developing GFLOW

This document describes how to set up a development environment to modify, build and test GFLOW.

## Prerequisites

Before you can build and test GFLOW, you must install and configure the following on your development machine:

- git
- pixi
- a modern GNU Fortran compiler

### Git

[Git](https://git-scm.com) and/or the **GitHub app** (for [Mac](https://mac.github.com) or [Windows](https://windows.github.com)).
[GitHub's Guide to Installing Git](https://help.github.com/articles/set-up-git) is a good source of information.

### Pixi

Pixi is a package management tool for developers. It allows the developer to
install libraries and applications in a reproducible way and works
cross-platform on Windows, Mac, and Linux.

Pixi installation docs can be found [here](https://pixi.sh).

### Fortran compiler

GNU Fortran can be used to build GFLOW. It may be possible to build GFLOW with other compilers, but this cannot be guaranteed.

#### GNU Fortran

GNU Fortran can be installed on all three major platforms.

*Linux*

- Fedora-based: `dnf install gcc-gfortran`
- Debian-based: `apt install gfortran`

*macOS*

- [Homebrew](https://brew.sh/): `brew install gcc@13`
- [MacPorts](https://www.macports.org/): `sudo port install gcc13`

**Note:** Xcode 15 includes a new linker implementation which breaks GNU Fortran compatibility. A workaround is to set `LDFLAGS` to use the classic linker, for instance:

```shell
export LDFLAGS="$LDFLAGS -Wl,-ld_classic"
```

See [this ticket](https://github.com/mesonbuild/meson/issues/12282) on the Meson repository for more information.

*Windows*

[Minimalist GNU for Windows](https://www.mingw-w64.org/) is the recommended way to obtain the GCC toolchain on Windows. Several MinGW distributions are available.

To install with Chocolatey: `choco install mingw`

To install from SourceForge:

- Download the MinGW installer:
  https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe
- Run the installer. Make sure to change `Architecture` to `x86_64`. Leave the
  other settings on default.
- Find the `mingw64/bin` directory in the installation and add it
  to your PATH. Find `Edit the system environment variables` in your Windows
  Start Screen. Click the `Environmental Variables` button and double-click the
  `Path` variable in the User Variables (the top table). Click the `New` button
  and enter the location of the `mingw64/bin` directory.

Binaries may also be downloaded and installed from the [releases here](https://github.com/brechtsanders/winlibs_mingw/releases).

**Note:** the MinGW distribution [available on conda-forge](https://anaconda.org/conda-forge/m2w64-toolchain_win-64) includes an outdated version of GCC.


## Get the GFLOW repository

Fork and clone the GFLOW repository:

1. Login to your GitHub account or create one by following the instructions given [here](https://github.com/signup/free).
2. [Fork](http://help.github.com/forking) the [main GFLOW](https://github.com/USEPA/gflow1).
3. Clone your fork of the GFLOW repository and create an `upstream` remote pointing back to your fork.

After forking the GFLOW repository on GitHub.

1. Clone your fork of the GitHub repository to your computer.

```shell
  git clone git@github.com:<github username>/gflow1.git
```

2. Go to the GFLOW directory.

```shell
cd gflow1
```

3. Add the main GFLOW repository as an upstream remote to your repository.

```shell
git remote add upstream https://github.com/MODFLOW-USGS/modflow6.git
```

## Building

[Meson](https://mesonbuild.com/index.html) is the recommended build system for
GFLOW. It is included in the provided pixi environment.

Building GFLOW requires two steps:

- configure the build directory
- build the project

We can use pixi to do so:

```shell
pixi run setup
```

```shell
pixi run build
```
