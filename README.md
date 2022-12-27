# g4mu
### Muon simulations in Geant4 for indirect Cherenkov light yield studies



#### Prerequisites
* Geant4 11.0.0+
  * This project used Geant4 11.1.0 downloaded from [here](https://gitlab.cern.ch/geant4/geant4/-/archive/v11.1.0/geant4-v11.1.0.tar.gz) and built from source.
  * Instructions to build Geant4: 
  
     ```
     wget https://gitlab.cern.ch/geant4/geant4/-/archive/v11.1.0/geant4-v11.1.0.tar.gz
     tar -xzvf geant4-v11.1.0
     mkdir geant4-v11.1.0-install geant4-v11.1.0-build
     cd geant4-v11.1.0-build
     cmake ../geant4-v11.1.0 -DCMAKE_INSTALL_PREFIX=../geant4-v11.1.0-install -DGEANT4_INSTALL_DATA=ON -DGEANT4_BUILD_TLS_MODEL=global-dynamic
     make -j8
     make install
     export GEANT_INSTALL_DIR=/path/to/geant4-v11.1.0-install
     export PATH=$PATH:$GEANT_INSTALL_DIR
     ```
  * To set the G4 environment variables (e.g. where the Geant4 data tables are stored), add `source /path/to/geant4-v11.1.0-install/bin/geant4.sh` to .bashrc or .bash_profile.
* [geant4_pybind](https://github.com/HaarigerHarald/geant4_pybind)
  * Requires core Geant4 to be already built on the system      
* ROOT
  * This project used ROOT 6.26/10 cloned from [here](https://github.com/root-project/root.git) and built from source.
      

#### How to run the simulation

1. Make an output folder in the g4mu directory:
```
cd g4mu
mkdir output
```
2. Run the simulation script from the g4mu directory as follows:
```
python simulation/propagate_muons_root_output_geant4_pybind.py -o output/icesim.root -n 1000 -e 20000 -rs 42
```
Arguments:
- `-n` is the number of events
- `-e` is the projectile energy in MeV
- `-rs` is the random seed (optional)
- `so` is the source particle name (optional; default is 'mu-')

The above command will run the Geant4 simulation and save the detailed information (including secondary energies, track lengths, and processes that created them) into ROOT trees (one tree entry per event). 
