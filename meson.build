project(
    'gflow1',
    'fortran',
    version: '0.0.1.dev0',
    license: 'MIT',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'optimization=2',
    'debug=false',
])
  
compiler = meson.get_compiler('fortran')
compiler_id = compiler.get_id()
compile_args = []
if compiler_id == 'gcc'
  # General options
  compile_args += ['-std=legacy', '-fno-align-commons', '-freal-4-real-8']
  # Define OS with gfortran for OS specific code
  # These are identical to pre-defined macros available with ifort
  system = build_machine.system()
  if system == 'linux'
    compile_args += '-D__linux__'
  elif system == 'darwin'
    compile_args += '-D__APPLE__'
  elif system == 'windows'
    compile_args += '-D_WIN32'
  endif
# TODO ifort, etc.
endif


main_file = 'src/gfmod.for'
src = [
    'src/lahey_utils.f90',
    'src/bufio.for',
    'src/dbcheck.for',
    'src/dbfun.for',
    'src/dbio.for',
    'src/dbmat.for',
    'src/dbmod.for',
    'src/dbnear.for',
    'src/dbnflow.for',
    'src/dbserv.for',
    'src/dicheck.for',
    'src/difun.for',
    'src/diio.for',
    'src/dimat.for',
    'src/dimod.for',
    'src/dinear.for',
    'src/dinflow.for',
    'src/gffun.for',
    'src/gfgrid.for',
    'src/gfio.for',
    'src/gfmat.for',
    'src/gfnflow.for',
    'src/gfserv.for',
    'src/gftrace.for',
    'src/gvcheck.for',
    'src/gvfun.for',
    'src/gvio.for',
    'src/gvmat.for',
    'src/gvmod.for',
    'src/gvnflow.for',
    'src/gvserv.for',
    'src/lkcheck.for',
    'src/lkfun.for',
    'src/lkio.for',
    'src/lkmat.for',
    'src/lkmod.for',
    'src/lknflow.for',
    'src/lkserv.for',
    'src/lscheck.for',
    'src/lsfun.for',
    'src/lsio.for',
    'src/lsmat.for',
    'src/lsmod.for',
    'src/lsnear.for',
    'src/lsnflow.for',
    'src/lsserv.for',
    'src/lsstream.for',
    'src/pdcheck.for',
    'src/pdfun.for',
    'src/pdio.for',
    'src/pdmat.for',
    'src/pdmod.for',
    'src/pdnear.for',
    'src/pdnflow.for',
    'src/pdserv.for',
    'src/twfun.for',
    'src/twio.for',
    'src/twmod.for',
    'src/twnear.for',
    'src/twserv.for',
    'src/w3check.for',
    'src/w3fun.for',
    'src/w3io.for',
    'src/w3mat.for',
    'src/w3mod.for',
    'src/w3near.for',
    'src/w3nflow.for',
    'src/w3serv.for',
    'src/wlcheck.for',
    'src/wlfun.for',
    'src/wlio.for',
    'src/wlmat.for',
    'src/wlmod.for',
    'src/wlnear.for',
    'src/wlnflow.for',
    'src/wlserv.for',
]
incdir = include_directories('./src')
project_bin_dir = join_paths(meson.current_source_dir(), 'bin')
executable(
    'gflow2',
    sources: [main_file] + src,
    include_directories: incdir,
    fortran_args : compile_args,
    install: true,
    install_dir: project_bin_dir,
)