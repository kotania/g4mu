# g4mu
### Muon simulations in Geant4 for indirect Cherenkov light yield studies



#### Prerequisites
* Geant4 11.0.0+
  * This project used Geant4 11.1.0 downloaded from [here](https://gitlab.cern.ch/geant4/geant4/-/archive/v11.1.0/geant4-v11.1.0.tar.gz) and built from source.
* [g4python](https://github.com/koichi-murakami/g4python.git)

To compile the geometry library, open the `CMakeLists.txt` file and set up the paths to the Geant4 and Geant4Py installations. Then, run:

```
mkdir build
cd build
cmake ..
make install
```
which will create the `geometry.so` object in the main directory.
