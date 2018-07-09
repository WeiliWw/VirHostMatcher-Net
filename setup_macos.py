import os  
import numpy
from distutils.sysconfig import get_config_vars
from distutils.core import setup, Extension  

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
cfg_vars = get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")




MOD = "tools"
version = "1.0"

setup(name=MOD, 
    version = version,
    ext_modules=[Extension(name = MOD, 
    sources = ["./src/tools/tools.cpp", "./src/tools/WIsH.cpp", "./src/tools/mm.cpp"], 
    # extra_compile_args=['-fPIC','-pthread','-std=c++11'])], 
    extra_compile_args=['-w','-stdlib=libc++','-std=c++11'])],
    include_dirs=[numpy.get_include(),"./src/tools/"]
    )  

