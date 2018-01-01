from distutils.core import setup, Extension

# the c++ extension module
extension_mod = Extension("hilbert", ["Hilbert.cpp", "Hilbert.h"])

setup(name = "hilbert", ext_modules=[extension_mod])