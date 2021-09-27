# Messinian Salinity Crisis [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5150871.svg)](https://doi.org/10.5281/zenodo.5150871) [![docs-python](https://img.shields.io/badge/docs-Python-green)](https://markmbaum.github.io/messinian-salinity-crisis/python/index.html) [![docs-c++](https://img.shields.io/badge/docs-C++-pink)](https://markmbaum.github.io/messinian-salinity-crisis/c++/index.html)

This repository contains code for simulating Mediterranean Sea level during the first stage of the [Messinian Salinity Crisis](https://en.wikipedia.org/wiki/Messinian_salinity_crisis), as detailed in
* [Baum, M. M. *Critical Analysis of Mediterranean Sea Level Limit Cycles During the Messinian Salinity Crisis*. Geologica Acta 11 (2021).
](https://revistes.ub.edu/index.php/GEOACTA/article/view/33997/35288)

My model is implemented in C++ and Python, with [documentation available here](https://markmbaum.github.io/messinian-salinity-crisis/). It is an adapted version of the model first described in
* [Garcia-Castellanos and Villasenor, *Messinian salinity crisis regulated by competing tectonics and erosion at the Gibraltar arc.* Nature 480.7377 (2011).](https://www.nature.com/articles/nature10651)

The original model code was recently made available [here](https://github.com/danigeos/asalted).

-----

My paper argues that the original model in the 2011 Nature paper is
1. **Irreproducible**. Starting from the same physical assumptions, the main results of the previous modeling cannot be produced.
2. **Incorrect**. It assumes a relationship between channel slope and Mediterranean sea-level that cannot be justified. Without this erroneous assumption, the main results of the paper (Mediterranean sea-level limit cycles) are probably impossible.
