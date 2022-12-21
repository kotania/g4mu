# g4mu
Muon simulations in Geant4 for indirect Cherenkov light yield studies

To compile the geometry library, open the `CMakeLists.txt` file and set up the paths to the Geant4 and Geant4Py installations. Then, run:

```mkdir build
cd build
cmake ..
make install
```
which will create the `geometry.so` object in the main directory.
