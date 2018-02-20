from distutils.core import setup, Extension

name = "tools"      # name of the module
version = "1.0"  # the module's version number
sources=["tools_wrap.cxx", "WIsH/main.cpp", "WIsH/mm.cpp", "CAFE/main.cpp", "CAFE/kmer.cpp", "CAFE/dist_model.cpp", "CAFE/output.cpp", "CAFE/utils.cpp", "CAFE/seq_model.cpp"]

setup(name=name, version=version,
      # distutils detects .i files and compiles them automatically
    ext_modules=[Extension(name='_tools', # SWIG requires _ as a prefix for the module name
    sources=sources,
    extra_compile_args=['-std=c++11','-w'],
    include_dirs=['WIsH','CAFE'])
    ])
