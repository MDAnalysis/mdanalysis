py3.extension_module('cutils',
  [cython_gen_cpp.process('cutils.pyx')],
  cpp_args: cython_cpp_args,
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/analysis/encore'
)

py3.install_sources([
'__init__.py',
'bootstrap.py',
'confdistmatrix.py',
'covariance.py',
'similarity.py',
'utils.py',
  ],
  subdir: 'MDAnalysis/analysis/encore'
)

subdir('clustering')
subdir('dimensionality_reduction')
