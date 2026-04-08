# sXPS
solve XPS: wrapper for lmfit(XPS)

## About

sXPS is a wrapper for [lmfit](https://github.com/lmfit/lmfit-py) and [lmfitxps](https://github.com/Julian-Hochhaus/lmfitxps) that aims to simplify my workflow for X-ray Photoelectron Spectroscopy (XPS) analysis.

sXPS introduces a `Data` class that reads excel files containing data copied from the [Avantage Software](https://www.thermofisher.com/be/en/home/electron-microscopy/products/xps-instruments/features.html#avantage-software). The data model can be build with `Data.add_singlet` and `Data.add_doublet` functions that wrap the corresponding models from [lmfitxps](https://github.com/Julian-Hochhaus/lmfitxps). A Shirley background is included in the model by default.

Finally, the model can be fitted to the data using `Data.fit`, which is a wrapper for [lmfit](https://github.com/lmfit/lmfit-py).

For further details, please refer to the documentation of [lmfit](https://lmfit.github.io/lmfit-py/index.html) and [lmfitxps](https://lmfitxps.readthedocs.io/en/latest/index.html).