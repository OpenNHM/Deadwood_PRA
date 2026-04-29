# Deadwood_PRA: Decision support tool for assessing the protective effect of lying deadwood against snow avalanche release
This script computes avalanche release membership and assesses the protective effect of lying deadwood dominated forest areas against snow avalanche release. It is designed for disturbed mountain forests where lying deadwood strongly influences surface roughness and therefore the likelihood of avalanche release. The workflow combines terrain information, deadwood structure, and canopy coverage to derive spatially explicit fuzzy membership maps.

The code requires only a dense point cloud (ideally derived from UAV photogrammetry or ULS) and a reference DTM.

To run the provided R code, the paths to the input folders must be adjusted.

For testing the code, we provide following datasets:

- Test raw data (point clouds): [Zenodo Dataset](https://doi.org/10.5281/zenodo.19485391)  
- Results (calculated with default settings): [Zenodo Results](https://doi.org/10.5281/zenodo.19734443)

Due to data restrictions we cannot provide the required reference DTM.  
Users are referred to the official regional data provider.

For the study area Kals am Großglockner, the reference_DTM can be downloaded here:  
[TIRIS data portal](https://www.tirol.gv.at/sicherheit/geoinformation/geodaten-tiris/laserscandaten/)

For details, we refer to the corresponding publication:
Bührle, L. J., Baggio, T., Winiwarter, L., Adams, M., Lingua, E., Richter, P., Schulte, F., Holstein, K., Bebi, P., Marke, T., & Teich, M. (in prep.). A spatially explicit UAV-based decision support tool to assess the protective effect of lying deadwood on snow avalanche release.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
