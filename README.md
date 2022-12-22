# g4mu
### Muon simulations in Geant4 for indirect Cherenkov light yield studies



#### Prerequisites
* Geant4 11.0.0+
  * This project used Geant4 11.1.0 downloaded from [here](https://gitlab.cern.ch/geant4/geant4/-/archive/v11.1.0/geant4-v11.1.0.tar.gz) and built from source.
* [g4python](https://github.com/koichi-murakami/g4python.git)
  * Requires core Geant4 to be already built on the system
  * IMPORTANT: prior to building g4python, do the following:
    - uncomment [lines 34-36](https://github.com/koichi-murakami/g4python/blob/8fc4b88e6c6ff80a4993d29794feee19ec766d5f/source/particles/pyG4PrimaryVertex.cc#L34-L36) in `g4python/source/particles/pyG4PrimaryVertex.cc`;
    - uncomment [lines 36-42](https://github.com/koichi-murakami/g4python/blob/main/source/particles/pyG4PrimaryParticle.cc#L36-L42) in `g4python/source/particles/pyG4PrimaryParticle.cc`;
  
    Otherwise the `simulation/propagate_muons_root_output.py` script from this repository will crash with the following error: 

       ```
       geant4.G4particles.G4PrimaryVertex(G4PrimaryParticle): No constructor defined!
       ```
    -  append [lines 60-61](https://github.com/wkcwells/g4py/blob/master/source/particles/pyG4PrimaryParticle.cc#L60-L61) and [line 69](https://github.com/wkcwells/g4py/blob/master/source/particles/pyG4PrimaryParticle.cc#L69) from the external [g4py](https://github.com/wkcwells/g4py) repository to the end of the `g4python/source/particles/pyG4PrimaryParticle.cc` file. 

    Otherwise the `simulation/propagate_muons_root_output.py` script from this repository will crash with the following error:
       ```
       'geant4.G4particles.G4PrimaryParticle' object has no attribute 'SetMass' ('Set4Momentum').
       ```
      

#### How to run the simulation

_Step 1._ Compile the geometry library. 

Open the `CMakeLists.txt` file and set up the paths to the Geant4 and Geant4Py installations. Then, run:

```
mkdir build
cd build
cmake ..
make install
```
which will create the `geometry.so` object in the main directory.
