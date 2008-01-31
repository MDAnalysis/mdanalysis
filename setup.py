import os
from distutils.core import setup, Extension
from distutils import sysconfig
from Pyrex.Distutils import build_ext

srcs = [['dcd/dcd.c'],['dcd/_dcdtest.pyx'],['numtools/_rms_matrix.pyx'],['delaunay/_delaunay.pyx', 'delaunay/blas.c', 'delaunay/tess.c']]
srcs = [[os.path.join('src', x) for x in s] for s in srcs]
include_dirs = ['src/dcd/include']

if __name__ == '__main__':
    DOC_FILES = ('COPYRIGHT', 'README', 'VERSION', 'CHANGELOG', 'KNOWN_BUGS', 'MAINTAINERS', 'TODO')
    LONG_DESCRIPTION = \
"""MDAnalysis is a tool for analyzing molecular dynamics trajectories.
"""
    CLASSIFIERS = ['Development Status :: 1 - Alpha',
                   'Environment :: Console',
                   'Intended Audience :: Scientists',
                   'License :: OSI Approved :: BSD License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Topic :: Scientific Software :: Biology',
                   'Topic :: Scientific Software :: Chemistry',]
    install_dir = sysconfig.get_python_lib() + os.sep + 'MDAnalysis'

    if 'DEBUG_CFLAGS' in os.environ:
        extra_compile_args = '\
            -std=c99 -pedantic -Wall -Wcast-align -Wcast-qual -Wpointer-arith \
            -Wchar-subscripts -Winline -Wnested-externs -Wbad-function-cast \
            -Wunreachable-code -Werror'
    else:
        extra_compile_args = ''

    define_macros = [('DEBUG', '1')]
    extensions = [Extension('_dcdmodule', srcs[0],
                            include_dirs = include_dirs,
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('_dcdtest', srcs[1],
                            include_dirs = include_dirs,
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('distances', ['src/numtools/distances.pyx'],
                            include_dirs = [],
                            libraries = ['m'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  #Extension('_rms_matrix', srcs[2],
                  #          include_dirs = include_dirs,
                  #          libraries = ['m'],
                  #          define_macros=define_macros,
                  #          extra_compile_args=extra_compile_args,
                  #          extra_link_args=["-framework","vecLib"]),
                  #Extension('_delaunay', srcs[3],
                  #          include_dirs = include_dirs,
                  #          libraries = ['m'],
                  #          define_macros=define_macros,
                  #          extra_compile_args=extra_compile_args)
                  ]

    setup(name              = 'MDAnalysis',
          version           = '0.1',
          description       = 'Python tools to support analysis of trajectories',
          author            = 'Naveen Michaud-Agrawal',
          author_email      = 'nmichaud@jhu.edu',
          url               = 'http://www.google.com',
          license           = 'BSD',
          packages          = [ 'MDAnalysis' ],
          package_dir       = {'MDAnalysis': 'python'},
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          #data_files        = [ (install_dir, DOC_FILES) ],
          #classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          cmdclass = {'build_ext': build_ext}
          )
