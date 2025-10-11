py3.extension_module('cython_util',
  [cython_gen_cpp.process('cython_util.pyx')],
  cpp_args: cython_cpp_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/lib/formats'
)

py3.extension_module('libdcd',
  [cython_gen_cpp.process('libdcd.pyx')],
  cpp_args: cython_cpp_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/lib/formats'
)

xdr_src = [
    'src/xdrfile_trr.c',
    'src/xdrfile_xtc.c',
    'src/xdrfile.c',
    'src/trr_seek.c',
    'src/xtc_seek.c',
]

py3.extension_module('libmdaxdr',
  xdr_src + [cython_gen_cpp.process('libmdaxdr.pyx')],
  cpp_args: cython_cpp_args,
  include_directories: [
    'include',
    ],
  dependencies: np_dep,
  install: true,
  subdir: 'MDAnalysis/lib/formats'
)

py3.install_sources([
    '__init__.py',
  ],
  subdir: 'MDAnalysis/lib/formats'
)
