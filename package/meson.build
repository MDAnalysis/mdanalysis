project(
  'MDAnalysis',
  'c', 'cython',
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
cc = meson.get_compiler('c')
cy = meson.get_compiler('cython')
cython = find_program(cy.cmd_array()[0])

if host_machine.system() == 'darwin' and cc.get_id() == 'clang'
  add_project_link_arguments('-stdlib=libc++,-mmacosx-version-min=10.9',
                             language : ['c'])
endif

subdir('MDAnalysis')
