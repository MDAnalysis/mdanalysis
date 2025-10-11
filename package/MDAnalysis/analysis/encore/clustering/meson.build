affinityprop_src = ['src/ap.c']

py3.extension_module('affinityprop',
  affinityprop_src + [cython_gen.process('affinityprop.pyx')],
  cpp_args: cython_c_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/analysis/encore/clustering'
)

py3.install_sources([
    'ClusterCollection.py',
    'ClusteringMethod.py',
    '__init__.py',
    'cluster.py',
  ],
  subdir: 'MDAnalysis/analysis/encore/clustering'
)
