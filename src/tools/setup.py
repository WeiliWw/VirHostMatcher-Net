import os  
import numpy
from distutils.sysconfig import get_config_vars
from distutils.core import setup, Extension  

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
cfg_vars = get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")




MOD = 'tools'  
setup(name=MOD, ext_modules=[Extension(MOD, 
    sources = ["tools.cpp", "WIsH.cpp", "mm.cpp"], 
    # extra_compile_args=['-fPIC','-pthread','-std=c++11'])], 
    extra_compile_args=['-w','-std=c++11'])],
    include_dirs=[numpy.get_include()]
    )  

