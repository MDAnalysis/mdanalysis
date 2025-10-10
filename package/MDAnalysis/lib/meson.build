py3.extension_module('_cutil',
  [cython_gen_cpp.process('_cutil.pyx')],
  cpp_args: cython_cpp_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/lib'
)


py3.install_sources([
    'NeighborSearch.py',
    '__init__.py',
    'correlations.py',
    'log.py',
    'mdamath.py',
    'picklable_file_io.py',
    'pkdtree.py',
    'transformations.py',
    'util.py',
    '_distopia.py',
    'distances.py',
  ],
  subdir: 'MDAnalysis/lib'
)

