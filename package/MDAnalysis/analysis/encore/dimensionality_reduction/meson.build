affinityprop_src = ['src/spe.c']

py3.extension_module('stochasticproxembed',
  affinityprop_src + [cython_gen.process('stochasticproxembed.pyx')],
  c_args: cython_c_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/analysis/encore/dimensionality_reduction'
)

py3.install_sources([
'DimensionalityReductionMethod.py',
'__init__.py',
'reduce_dimensionality.py',
  ],
  subdir: 'MDAnalysis/analysis/encore/dimensionality_reduction'
)

