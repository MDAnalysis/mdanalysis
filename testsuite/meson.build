project(
  'MDAnalysisTests',
  'c', 'cpp', 'cython',
  version: '2.10.0-dev0',
  license: 'LGPLv3+',
  meson_version: '>= 1.9.0',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c17',
  ],
)

py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()
subdir('MDAnalysisTests')
