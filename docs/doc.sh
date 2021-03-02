#run Doxygen for the C++ implementation
cd ../model/c++
doxygen

#run Sphinx for the Python implementation
cd ../python/sphinx
make html
cp -r _build/html/* ../../../docs/python
cd ../../..
