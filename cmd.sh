rm -R build/
#swig -python -c++ -Isrc tools.i
python setup.py build_ext   
python setup.py install --install-platlib=./src/

