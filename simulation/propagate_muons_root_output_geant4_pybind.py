"""
An adaptation of Daniel Dwyer's argon_box.py code (https://github.com/dadwyer/argon_box.git) 
for muon simulations in ice.

Author: Tania Kozynets
tetiana.kozynets@nbi.ku.dk
December 22, 2022

"""


# General tools
from copy import copy
import numpy as np
import sys

sys.path.append("../")

# Generators
from geant4_pybind import G4VUserPrimaryGeneratorAction, G4ParticleGun
from geant4_pybind import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex

# User-geant4 actions
from geant4_pybind import (
    G4UserRunAction,
    G4UserEventAction,
    G4UserSteppingAction,
    G4UserTrackingAction,
)

# Physics processes
from geant4_pybind import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QBBC

# Detector construction
from geant4_pybind import G4VUserDetectorConstruction

# Geant4 geometry
from geometry import icegeo

# Run Manager and UI commands
from geant4_pybind import G4RunManagerFactory, G4RunManagerType

from geant4_pybind import G4UImanager

# Random numbers
from geant4_pybind import HepRandom

# Physical Units
from geant4_pybind import MeV, cm, g, GeV, m, ns
from geant4_pybind import G4ThreeVector

# G4 Configuration
from geant4_pybind import G4UserLimits


# ROOT
from array import array
import ROOT
from ROOT import TTree, TFile

# Command line arguments
from argparse import ArgumentParser

SOURCE_POSITION = G4ThreeVector(0, 0, 0)
INIT_DIRECTION = G4ThreeVector(0, 0, 1)


parser = ArgumentParser()
parser.add_argument("-n", "--nevents", type=int, help="Number of events", default=10)
parser.add_argument("-so", "--source", help="Particle name or generator", default="mu-")
parser.add_argument(
    "-e", "--energy", type=float, help="Particle energy in MeV", default=20.0 * GeV
)
parser.add_argument("-o", "--output", help="Name of the ROOT output file")
parser.add_argument("-rs", "--seed", type=int, help="Random seed", default=312)
parser.add_argument(
    "--physlist",
    help="G4 Physics List: FTFP_BERT, QGSP_BERT, QGSP_BERT_HP",
    default="QGSP_BERT",
)

args = parser.parse_args()


class Simulation:
    def __init__(self):
        """Constructor"""
        self._ofilename = args.output
        self._random_seed = args.seed
        self._ofile = None
        self._otree = None
        self._treebuffer = None
        self._geom = None
        self._physlist_name = args.physlist
        self._physics_list = None
        self._source = None
        self._energies = None
        self._source_position = None
        self._generator = None
        self.run_manager = G4RunManagerFactory.CreateRunManager(G4RunManagerType.Serial)
        # return

    def initialize(self):
        """Initialize the simulation"""
        # Parse command line arguments
        self._parse_args()
        # Prepare output file
        self._prepare_output()
        # Initialize random number generator
        self._init_random()
        # Prepare geometry
        self._init_geometry()
        # Prepare physics list
        self._init_physics_list()
        # Prepare generator
        self._init_generator()
        # Prepare user actions
        self._init_user_actions()
        # Prepare run manager
        self._init_run_manager()
        # return

    def run(self):

        print("About to generate %s events..." % self._nevents)
        self.run_manager.BeamOn(self._nevents)

    def finalize(self):
        print("Finalizing")
        self._close_output()

    ######################

    def _parse_args(self):
        """Parse command line arguments"""

        self._args = args
        if self._args.nevents is not None:
            self._nevents = self._args.nevents
        if self._args.source is not None:
            self._source = self._args.source
        if self._args.energy is not None:
            self._energies = [
                self._args.energy * GeV,
            ]
        if self._args.output is not None:
            self._ofilename = self._args.output
        if self._args.seed is not None:
            self._random_seed = self._args.seed
        print("Configuration:")
        print("  nevents:", self._nevents)
        print("   source:", self._source)
        print("   energies:", self._energies)
        print("   output:", self._ofilename)
        print(" physlist:", self._physlist_name)
        print("     seed:", self._random_seed)
        # return

    def _prepare_output(self):
        print("Prepare output file...")
        """Prepare ROOT tree for output data"""
        ofile = TFile(self._ofilename, "RECREATE")
        otree = TTree("muons_in_ice", "Muon Simulation")
        tb = TreeBuffer()
        tb.maxInit = int(1e6)
        tb.maxTracks = int(1e7)
        tb.ev = array("i", [0])
        # Geant4 initial state particles
        tb.ni = array("i", [0])
        tb.pidi = array("i", [0] * tb.maxInit)
        tb.xi = array("d", [0] * tb.maxInit)
        tb.yi = array("d", [0] * tb.maxInit)
        tb.zi = array("d", [0] * tb.maxInit)
        tb.ti = array("d", [0] * tb.maxInit)
        tb.pxi = array("d", [0] * tb.maxInit)
        tb.pyi = array("d", [0] * tb.maxInit)
        tb.pzi = array("d", [0] * tb.maxInit)
        tb.ekini = array("d", [0] * tb.maxInit)
        tb.mi = array("d", [0] * tb.maxInit)
        # Geant4 track step data
        tb.nstep = array("i", [0])
        tb.tid = array("i", np.zeros(tb.maxTracks, dtype="int64"))
        tb.pid = array("i", np.zeros(tb.maxTracks, dtype="int64"))
        tb.parid = array("i", np.zeros(tb.maxTracks, dtype="int64"))
        proc_strings = ROOT.vector("string")()
        tb.proc = proc_strings
        tb.track_length = array("d", np.zeros(tb.maxTracks, dtype="float64"))
        tb.ekin = array("d", np.zeros(tb.maxTracks, dtype="float64"))
        # Hook for ancestor particle
        tb.ancestor = None

        otree.Branch("ev", tb.ev, "ev/I")
        otree.Branch("ni", tb.ni, "ni/I")
        otree.Branch("pidi", tb.pidi, "pidi[ni]/I")
        otree.Branch("xi", tb.xi, "xi[ni]/D")
        otree.Branch("yi", tb.yi, "yi[ni]/D")
        otree.Branch("zi", tb.zi, "zi[ni]/D")
        otree.Branch("ti", tb.ti, "ti[ni]/D")
        otree.Branch("pxi", tb.pxi, "pxi[ni]/D")
        otree.Branch("pyi", tb.pyi, "pyi[ni]/D")
        otree.Branch("pzi", tb.pzi, "pzi[ni]/D")
        otree.Branch("ekini", tb.ekini, "ekini[ni]/D")
        otree.Branch("mi", tb.mi, "mi[ni]/D")
        otree.Branch("nstep", tb.nstep, "nstep/I")
        otree.Branch("tid", tb.tid, "tid[nstep]/I")
        otree.Branch("pid", tb.pid, "pid[nstep]/I")
        otree.Branch("parid", tb.parid, "parid[nstep]/I")
        otree.Branch("proc", tb.proc)
        otree.Branch("track_length", tb.track_length, "track_length[nstep]/D")
        otree.Branch("ekin", tb.ekin, "ekin[nstep]/D")

        self._ofile = ofile
        self._otree = otree
        self._treebuffer = tb
        # return

    def _close_output(self):
        print("Closing output file")
        """Write ROOT data and close file"""
        self._otree.Write()
        self._ofile.Close()
        # return

    def _init_random(self):
        """Initialize random number generator"""
        HepRandom.setTheSeed(self._random_seed)
        # return

    def _init_geometry(self):
        print("Initializing geometry")
        self._geom = icegeo.IceCubeDetectorConstruction()
        self.run_manager.SetUserInitialization(self._geom)
        # return

    def _init_physics_list(self):
        print("Initializing physics list")
        if self._physlist_name == "FTFP_BERT":
            self._physics_list = FTFP_BERT()
        elif self._physlist_name == "QGSP_BERT":
            self._physics_list = QGSP_BERT()
        elif self._physlist_name == "QGSP_BERT_HP":
            self._physics_list = QGSP_BERT_HP()
        elif self._physlist_name == "QBBC":
            self._physics_list = QBBC()
        else:
            raise ValueError(
                "Invalid physics list: '%r'.  Please choose from:FTFP_BERT, QGSP_BERT, QGSP_BERT_HP, QBBC"
                % self._physlist_name
            )
        self.run_manager.SetUserInitialization(self._physics_list)
        # return

    def _init_generator(self):
        print("Initializing generator")
        self._source_pos = SOURCE_POSITION
        self._generator = PrimaryGeneratorAction(self._treebuffer)
        self.run_manager.SetUserAction(self._generator)
        # return

    def _init_user_actions(self):

        print("Initializing user actions")
        self._run_action = RunAction()
        self._event_action = EventAction(
            treebuffer=self._treebuffer, outputtree=self._otree
        )
        self._step_action = SteppingAction(
            treebuffer=self._treebuffer, geometry=self._geom
        )

        self.run_manager.SetUserAction(self._run_action)
        self.run_manager.SetUserAction(self._event_action)
        self.run_manager.SetUserAction(self._step_action)

    def _init_run_manager(self):
        print("Initializing run manager")
        self.run_manager.Initialize()
        # return


###################################################################

# Tree Buffer
class TreeBuffer:
    """Dummy class for collecting tree data buffers"""

    pass


###################################################################


class PrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
    def __init__(self, treebuffer):
        super().__init__()

        self.isInitialized = False
        self.fParticleGun = G4ParticleGun(1)
        self._tb = treebuffer
        # default particle kinematic
        particleTable = G4ParticleTable.GetParticleTable()
        self.particle_name = args.source
        particle = particleTable.FindParticle(self.particle_name)
        self.fParticleGun.SetParticleDefinition(particle)
        self.fParticleGun.SetParticleMomentumDirection(INIT_DIRECTION)
        self.fParticleGun.SetParticleEnergy(args.energy)

    def initialize(self):

        particleTable = G4ParticleTable.GetParticleTable()
        self.particleDef = particleTable.FindParticle(self.particle_name)
        self.isInitialized = True

    def GeneratePrimaries(self, anEvent):
        # this function is called at the begining of each event

        if not self.isInitialized:
            self.initialize()
        # No ancestor for this generator
        self._tb.ancestor = None

        position = SOURCE_POSITION
        time = 0.0
        self.fParticleGun.SetParticlePosition(position)
        self.fParticleGun.GeneratePrimaryVertex(anEvent)
        mass = self.particleDef.GetPDGMass()
        tb = self._tb
        tb.pidi[0] = self.particleDef.GetPDGEncoding()
        tb.xi[0] = position.x / cm
        tb.yi[0] = position.y / cm
        tb.zi[0] = position.z / cm
        tb.ti[0] = time / ns
        tb.pxi[0] = 0
        tb.pyi[0] = 0
        tb.pzi[0] = np.sqrt(args.energy**2 - mass**2) / GeV
        tb.ekini[0] = (args.energy - mass) / GeV
        tb.mi[0] = mass / GeV


class RunAction(G4UserRunAction):

    print("*** Beginning of the run ***")

    def EndOfRunAction(self, run):

        print("*** End of the run ***")
        print(
            "- Run summary : (id= %d, #events= %d)"
            % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
        )


# ------------------------------------------------------------------
class EventAction(G4UserEventAction):
    def __init__(self, treebuffer, outputtree):

        """Constructor"""
        G4UserEventAction.__init__(self)
        self._tb = treebuffer
        self._otree = outputtree

    def BeginOfEventAction(self, event):

        self._tb.ni[0] = 0
        self._tb.nstep[0] = 0
        # return

    def EndOfEventAction(self, event):
        """Record event"""

        tb = self._tb
        # Event ID
        tb.ev[0] = event.GetEventID()
        # Primary particles
        ni = 0
        for pv_idx in range(event.GetNumberOfPrimaryVertex()):
            pvtx = event.GetPrimaryVertex(pv_idx)
            for pp_idx in range(pvtx.GetNumberOfParticle()):
                ppart = pvtx.GetPrimary(pp_idx)
                tb.pidi[ni] = ppart.GetPDGcode()
                tb.xi[ni] = pvtx.GetX0() / cm
                tb.yi[ni] = pvtx.GetY0() / cm
                tb.zi[ni] = pvtx.GetZ0() / cm
                tb.ti[ni] = pvtx.GetT0() / ns
                tb.pxi[ni] = ppart.GetPx() / GeV
                tb.pyi[ni] = ppart.GetPy() / GeV
                tb.pzi[ni] = ppart.GetPz() / GeV
                # Ick, G4py interface doesn't provide particle energy func :/
                pmomSq = ppart.GetPx() ** 2 + ppart.GetPy() ** 2 + ppart.GetPz() ** 2
                pmass = ppart.GetMass()
                penergy = np.sqrt(pmomSq + pmass**2)
                tb.ekini[ni] = (penergy - pmass) / GeV
                tb.mi[ni] = pmass / GeV
                ni += 1
                if ni == tb.maxInit:
                    print("Warning: Reached max number of initial particles.")
                    break
        tb.ni[0] = ni
        self._otree.Fill()
        tb.proc.clear()


# ------------------------------------------------------------------
class SteppingAction(G4UserSteppingAction):
    def __init__(self, treebuffer, geometry):
        """Constructor"""
        G4UserSteppingAction.__init__(self)
        self._tb = treebuffer
        self._geom = geometry

    def UserSteppingAction(self, step):
        """Collect data for current simulation step"""
        tb = self._tb
        istp = tb.nstep[0]
        if istp >= tb.maxTracks:
            print("Reached maximum tracks:", istp)
            return

        track = step.GetTrack()
        prestep = step.GetPreStepPoint()
        poststep = step.GetPostStepPoint()
        tb.tid[istp] = track.GetTrackID()
        tb.pid[istp] = track.GetDefinition().GetPDGEncoding()
        tb.parid[istp] = track.GetParentID()
        tb.ekin[istp] = prestep.GetKineticEnergy() / GeV
        process = track.GetCreatorProcess()
        if process == None:
            process = "None"
        else:
            process = str(process.GetProcessName())
        tb.proc.push_back(process)
        tb.track_length[istp] = track.GetTrackLength()
        tb.nstep[0] += 1


###################################################################

if "__main__" == __name__:
    # Run the simulation

    sim = Simulation()
    sim.initialize()
    sim.run()
    sim.finalize()
