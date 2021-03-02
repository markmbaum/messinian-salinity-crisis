cd model/c++
python clean.py out
make
./bin/single.exe out
cd scripts
python plot_out.py

cd ../../python
python reference_model.py
cd ../..
