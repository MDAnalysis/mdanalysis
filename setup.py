import sys, os
from distutils import sysconfig
from numpy import get_numpy_include
from distutils.core import setup, Extension
from Pyrex.Distutils import build_ext

include_dirs = [get_numpy_include()]

if sys.platform == "darwin": # Mac OS X
    fast_numeric_include = ['/System/Library/Frameworks/vecLib.framework/Versions/A/Headers']
    fast_numeric_link = ["-framework","vecLib"]
elif sys.platform[:5] == "linux":
    fast_numeric_include = ['/opt/intel/cmkl/8.0/include']
    fast_numeric_link = ["-L/opt/intel/cmkl/8.0/lib/32", "-lmkl_lapack","-lmkl_lapack32","-lmkl_ia32","-lmkl","-lguide"]
else:
    fast_numeric_include = []
    fast_numeric_link = []

if __name__ == '__main__':
    DOC_FILES = ('README', 'LICENSE', 'CHANGELOG', 'TODO')
    LONG_DESCRIPTION = \
"""MDAnalysis is a tool for analyzing molecular dynamics trajectories.
"""
    CLASSIFIERS = ['Development Status :: 1 - Alpha',
                   'Environment :: Workstation',
                   'Intended Audience :: Scientists',
                   'License :: OSI Approved :: GPL License',
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
        define_macros = [('DEBUG', '1')]
    else:
        extra_compile_args = ''
        define_macros = []

    extensions = [Extension('_dcdmodule', ['src/dcd/dcd.c'],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('_dcdtest', ['src/dcd/_dcdtest.pyx'],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('distances', ['src/numtools/distances.pyx'],
                            include_dirs = include_dirs+['src/numtools'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  #Extension('rms_fitting', ['src/numtools/rms_fitting.pyx'],
                  #          libraries = ['m'],
                  #          define_macros=define_macros,
                  #          include_dirs = include_dirs+fast_numeric_include,
                  #          extra_link_args=fast_numeric_link,
                  #          extra_compile_args=extra_compile_args),
                  #Extension('delaunay', ['src/delaunay/delaunay.pyx', 'src/delaunay/blas.c', 'src/delaunay/tess.c'],
                  #          libraries = ['m'],
                  #          define_macros=define_macros,
                  #          include_dirs = include_dirs+fast_numeric_include+['src/delaunay'],
                  #          extra_link_args=fast_numeric_link,
                  #          extra_compile_args=extra_compile_args),
                  ]

    setup(name              = 'MDAnalysis',
          version           = '0.3',
          description       = 'Python tools to support analysis of trajectories',
          author            = 'Naveen Michaud-Agrawal',
          author_email      = 'nmichaud@jhu.edu',
          url               = '',
          license           = 'GPL',
          packages          = [ 'MDAnalysis' ],
          package_dir       = {'MDAnalysis': 'python'},
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          data_files        = [ (install_dir, DOC_FILES) ],
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          cmdclass = {'build_ext': build_ext}
          )
