incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(py3,
    [
      '-c',
      '''import os
import numpy as np
try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
  '''
    ],
    check: true
  ).stdout().strip()

  # We do need an absolute path to feed to `cc.find_library` below
  _incdir_numpy_abs = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check: true
  ).stdout().strip()
else
  _incdir_numpy_abs = incdir_numpy
endif
inc_np = include_directories(incdir_numpy)
numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']
np_dep = declare_dependency(include_directories: inc_np, compile_args: numpy_nodepr_api)

_cython_tree = [fs.copyfile('__init__.py')]
cython_args = ['-3', '--fast-fail', '--output-file', '@OUTPUT@', '--include-dir', '@BUILD_ROOT@', '@INPUT@']
cython_shared_module = []
cython_c_args = ['-DCYTHON_CCOMPLEX=0']
cython_cpp_args = cython_c_args
cython_cplus_args = ['--cplus'] + cython_args

cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, cython_shared_module]
)

cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, cython_shared_module]
)


python_sources = [
  'due.py',
  'exceptions.py',
  '__init__.py',
  'units.py',
  'version.py',
]

py3.install_sources(
  python_sources,
  subdir: 'MDAnalysis'
)

subdir('auxiliary')
subdir('guesser')
subdir('selections')
subdir('transformations')
subdir('visualization')
subdir('tests')
subdir('lib')
subdir('analysis')
subdir('converters')
subdir('coordinates')
subdir('core')
subdir('topology')
