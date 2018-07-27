#! /usr/bin/env python

import os
from distutils.core import setup
from glob import glob
from distutils.extension import Extension

incdir_src = os.path.abspath("/direct/phenix+u/arbint/Jets/Analysis/simulation/LHAPDF-6.1.5/include")
incdir_build = os.path.abspath("../../include")
libdir = os.path.abspath("../../src/.libs")

ext = Extension("lhapdf",
                ["lhapdf.cpp"],
                include_dirs = [incdir_src, incdir_build],
                extra_compile_args= " -I/direct/phenix+u/arbint/Jets/Analysis/install/include".split(),
                library_dirs = [libdir],
                language = "C++",
                libraries = ["stdc++", "LHAPDF"])

setup(name = "LHAPDF",
      version = "6.1.5",
      ext_modules = [ext])
