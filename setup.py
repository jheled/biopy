from distutils.core import setup, Extension
import numpy, glob

module1 = Extension('biopy.cchelp',
                    include_dirs = [numpy.get_include()],
                    sources = ['biopy/cchelp.cc'])

module2 = Extension('biopy.cnexus',
                    sources = ['biopy/cnexus.c'])

module3 = Extension('biopy.treesset',
                    sources = ['biopy/treesset.cc'],
                    extra_compile_args=['-std=c++0x', '-Wno-invalid-offsetof'])

module4 = Extension('biopy.neutralsim',
                    sources = ['biopy/neutralsim.cc'],
                    extra_compile_args=['-std=c++0x'])

module5 = Extension('biopy.calign',
                    sources = ['biopy/calign.cc'],
                    extra_compile_args=['-std=c++0x'])

module7 = Extension('biopy.aalign',
                    sources = ['biopy/aalign.cc'],
                    extra_compile_args=['-std=c++0x'])


module6 = Extension('biopy.cclust',
                    sources = ['biopy/cclust.cc'],
                    extra_compile_args=['-std=c++0x'])

classifiers=[
  "Development Status :: 3 - Alpha",
  "Environment :: Console",
  "Intended Audience :: Developers",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: GNU Affero General Public License v3",
  "Programming Language :: Python",
  "Programming Language :: C",
  "Topic :: Scientific/Engineering :: Bio-Informatics"
  ]

from biopy import __version__

setup (name = 'biopy',
       version = __version__,
       description = 'Bioinformatics Python utilities',
       long_description = 'Bioinformatics Python utilities',
       author = 'Joseph Heled',
       author_email = 'jheled@gmail.com',
       url = 'http://http://code.google.com/p/biopy/',
       license = 'LGPL (V3)',
       classifiers = classifiers,
       platforms = ["Linux", "Mac OS-X"],
       packages = ['biopy'],
       package_dir={'biopy': 'biopy'},
       ext_modules = [module1,module2,module3,module4,module5,module6,module7],
       scripts = glob.glob('scripts/*'),
       data_files=[('doc',glob.glob('html/*.*')),
                   ('doc/_images', glob.glob('html/_images/*.*')),
                   ('doc/_images/math', glob.glob('html/_images/math/*.*')),
                   ('doc/_static', glob.glob('html/_static/*.*')),
                   ('biopy', ['biopy/readseq.h','biopy/seqslist.cc'])]
       )
