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
from geant4 import G4VUserPrimaryGeneratorAction, G4ParticleGun
from geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex

# User-geant4 actions
from geant4 import G4UserRunAction, G4UserEventAction, G4UserSteppingAction

# Physics processes
from geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP

# Detector construction
from geant4 import G4VUserDetectorConstruction

# Geant4 geometry
import geometry

# Run Manager and UI commands
from geant4 import gRunManager, gApplyUICommand, G4UImanager

# Random numbers
from geant4 import HepRandom

# Physical Units
from geant4 import MeV, cm, g, GeV, m, ns
from geant4 import G4ThreeVector

# G4 Configuration
from geant4 import G4UserLimits


# ROOT
from array import array
from ROOT import TTree, TFile

# Command-line arguments
from argparse import ArgumentParser


class MySimulation:
    "My Simple Simulation"

    def __init__(self):
        """Constructor"""
        self._ofilename = "ice_sim.root"
        self._random_seed = 23
        self._ofile = None
        self._otree = None
        self._treebuffer = None
        self._geom = None
        self._physlist_name = "QGSP_BERT"
        self._physics_list = None
        self._source = None
        self._energies = None
        self._source_position = None
        self._generator = None
        return

    def initialize(self):
        """Initialize the simulation"""
        # Parse command-line arguments
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
        return

    def run(self):
        """Run the simulation"""
        gRunManager.BeamOn(0)
        # gRunManager.BeamOn(self._nevents)
        return

    def finalize(self):
        print("finalize")
        """Finalize the simulation"""
        self._close_output()
        return

    ######################

    def _parse_args(self):
        """Parse command line arguments"""
        parser = ArgumentParser()
        parser.add_argument("--nevents", type=int, help="Number of events", default=1)
        parser.add_argument("--source", help="Particle name or generator", default="e-")
        parser.add_argument(
            "--energy", type=float, help="Particle energy", default=2.0 * GeV
        )
        parser.add_argument("--output", help="Output filename")
        parser.add_argument("--seed", type=int, help="Random seed")
        parser.add_argument(
            "--physlist",
            help="G4 Physics List: FTFP_BERT, QGSP_BERT, QGSP_BERT_HP",
            default="QGSP_BERT",
        )

        self._args = parser.parse_args()
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
        return

    def _prepare_output(self):
        print("prepare output")
        """Prepare ROOT tree for output data"""
        ofile = TFile(self._ofilename, "RECREATE")
        otree = TTree("muons_in_ice", "Muon Simulation")
        tb = TreeBuffer()

        tb.maxInit = 100
        tb.maxTracks = 100000
        tb.ev = array("i", [0])
        # Ancestor particles (e.g. muon)
        # Note: For simplicity, assume only one potential ancestor
        tb.pida = array("i", [0])
        tb.xa = array("d", [0])
        tb.ya = array("d", [0])
        tb.za = array("d", [0])
        tb.ta = array("d", [0])
        tb.pxa = array("d", [0])
        tb.pya = array("d", [0])
        tb.pza = array("d", [0])
        tb.ekina = array("d", [0])
        tb.ma = array("d", [0])
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
        tb.tid = array("i", [0] * tb.maxTracks)
        tb.pid = array("i", [0] * tb.maxTracks)
        tb.parid = array("i", [0] * tb.maxTracks)
        tb.xs = array("d", [0] * tb.maxTracks)
        tb.ys = array("d", [0] * tb.maxTracks)
        tb.zs = array("d", [0] * tb.maxTracks)
        tb.xe = array("d", [0] * tb.maxTracks)
        tb.ye = array("d", [0] * tb.maxTracks)
        tb.ze = array("d", [0] * tb.maxTracks)
        tb.ekin = array("d", [0] * tb.maxTracks)
        tb.edep = array("d", [0] * tb.maxTracks)
        # Hook for ancestor particle
        tb.ancestor = None

        otree.Branch("ev", tb.ev, "ev/I")
        otree.Branch("pida", tb.pida, "pida/I")
        otree.Branch("xa", tb.xa, "xa/D")
        otree.Branch("ya", tb.ya, "ya/D")
        otree.Branch("za", tb.za, "za/D")
        otree.Branch("ta", tb.ta, "ta/D")
        otree.Branch("pxa", tb.pxa, "pxa/D")
        otree.Branch("pya", tb.pya, "pya/D")
        otree.Branch("pza", tb.pza, "pza/D")
        otree.Branch("ekina", tb.ekina, "ekina/D")
        otree.Branch("ma", tb.ma, "ma/D")
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
        otree.Branch("ekin", tb.ekin, "ekin[nstep]/D")
        otree.Branch("edep", tb.edep, "edep[nstep]/D")
        otree.Branch("xs", tb.xs, "xs[nstep]/D")
        otree.Branch("ys", tb.ys, "ys[nstep]/D")
        otree.Branch("zs", tb.zs, "zs[nstep]/D")
        otree.Branch("xe", tb.xe, "xe[nstep]/D")
        otree.Branch("ye", tb.ye, "ye[nstep]/D")
        otree.Branch("ze", tb.ze, "ze[nstep]/D")

        self._ofile = ofile
        self._otree = otree
        self._treebuffer = tb
        return

    def _close_output(self):
        print("close output")
        """Write ROOT data, and close file"""
        self._otree.Write()
        self._ofile.Close()
        return

    def _init_random(self):
        print("init random")
        """Initialize random number generator"""
        HepRandom.setTheSeed(self._random_seed)
        return

    def _init_geometry(self):
        print("init geometry")
        """Initialize geant geometry"""
        self._geom = geometry.geometry()
        gRunManager.SetUserInitialization(self._geom)
        return

    def _init_physics_list(self):
        print("init phys list")
        """Initialize the physics list"""
        # Use standard physics list
        if self._physlist_name == "FTFP_BERT":
            self._physics_list = FTFP_BERT()
        elif self._physlist_name == "QGSP_BERT":
            self._physics_list = QGSP_BERT()
        elif self._physlist_name == "QGSP_BERT_HP":
            self._physics_list = QGSP_BERT_HP()
        else:
            raise ValueError(
                "Invalid physics list: '%r'.  Please choose from:FTFP_BERT, QGSP_BERT, QGSP_BERT_HP"
                % self._physlist_name
            )
        gRunManager.SetUserInitialization(self._physics_list)
        return

    def _init_generator(self):
        print("init generator")
        """Initialize particle generator"""
        self._source_pos = G4ThreeVector(0, 0, 0)
        self._generator = MyParticleGeneratorAction(
            particleName=self._source,
            energies=self._energies,
            position=self._source_pos,
            treebuffer=self._treebuffer,
        )
        gRunManager.SetUserAction(self._generator)
        return

    def _init_user_actions(self):

        print("init user actions")
        """Initialize user actions"""
        self._run_action = MyRunAction()
        self._event_action = MyEventAction(
            treebuffer=self._treebuffer, outputtree=self._otree
        )
        self._step_action = MySteppingAction(
            treebuffer=self._treebuffer, geometry=self._geom
        )
        gRunManager.SetUserAction(self._run_action)
        gRunManager.SetUserAction(self._event_action)
        gRunManager.SetUserAction(self._step_action)
        return

    def _init_run_manager(self):
        print("init run manager")
        """Initialize the Geant4 run manager"""
        gApplyUICommand("/run/initialize")
        return


###################################################################

# Tree Buffer
class TreeBuffer:
    """Dummy class for collecting tree data buffers"""

    pass


###################################################################

# Particle interaction generator
class MyParticleGeneratorAction(G4VUserPrimaryGeneratorAction):
    "Generator for single type of particles (e.g. e-, gammas, etc)"

    def __init__(
        self,
        treebuffer,
        particleName="mu-",
        energies=[
            1.0 * GeV,
        ],
        position=G4ThreeVector(0, 0, 0),
    ):
        print("init generator action (constructor)")
        G4VUserPrimaryGeneratorAction.__init__(self)
        self.isInitialized = False
        self.particleName = particleName
        self.energies = energies
        self.position = position
        self.particleDef = None
        self._tb = treebuffer
        pass

    def initialize(self):

        print("init generator action (method)")
        # Prepare generator
        particleTable = G4ParticleTable.GetParticleTable()
        print("generated particle table")
        self.particleDef = particleTable.FindParticle(self.particleName)
        print("found particle")
        self.isInitialized = True
        print("initialized the particle")
        return

    def GeneratePrimaries(self, event):

        print("generate primaries")
        # self.particleGun.GeneratePrimaryVertex(event)
        if not self.isInitialized:
            self.initialize()
        # No ancestor for this generator
        self._tb.ancestor = None
        # Create primaries
        position = self.GenerateVertexPosition()
        print("generated vertex position")
        time = 0.0
        vertex = G4PrimaryVertex(position, time)
        mass = self.particleDef.GetPDGMass()
        print("got pdg mass")
        # Record initial particle
        tb = self._tb
        tb.pidi[0] = self.particleDef.GetPDGEncoding()
        print("got pdg encoding")
        tb.xi[0] = position.x / cm
        tb.yi[0] = position.y / cm
        tb.zi[0] = position.z / cm
        tb.ti[0] = time / ns
        tb.pxi[0] = 0
        tb.pyi[0] = 0
        tb.pzi[0] = 0
        tb.ekini[0] = 0
        tb.mi[0] = mass / GeV
        print("set up initial particle properties")
        print(self.energies)
        for enr in self.energies:
            particle = G4PrimaryParticle(self.particleDef.GetPDGEncoding())
            # Particle emission is in +z direction
            particle.Set4Momentum(0, 0, np.sqrt(enr**2 - mass**2), enr)
            # Add to total initial momentum and kinetic energy
            tb.pzi[0] += np.sqrt((self.energies[0]) ** 2 - mass**2) / GeV
            tb.ekini[0] += (enr - mass) / GeV
            # Ensure mass is exact
            particle.SetMass(mass)
            # Set direction
            particleDirection = self.GenerateParticleDirection()
            particle.SetMomentumDirection(particleDirection)
            vertex.SetPrimary(particle)
        event.AddPrimaryVertex(vertex)
        print("added primary vertex")

    def GenerateVertexPosition(self):
        """Generate a vertex position"""
        return G4ThreeVector(self.position)

    def GenerateParticleDirection(self):
        """Generate a particle direction"""
        return G4ThreeVector(0, 0, 1)


class MyRunAction(G4UserRunAction):
    "My Run Action"

    def EndOfRunAction(self, run):

        print("end of run action")
        print("*** End of Run")
        print(
            "- Run sammary : (id= %d, #events= %d)"
            % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())
        )


# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
    "My Event Action"

    def __init__(self, treebuffer, outputtree):

        print("init event action")
        """Constructor"""
        G4UserEventAction.__init__(self)
        self._tb = treebuffer
        self._otree = outputtree

    def BeginOfEventAction(self, event):

        print("begin of event action")
        self._tb.ev[0] = -1
        self._tb.pida[0] = 0
        self._tb.xa[0] = 0
        self._tb.ya[0] = 0
        self._tb.za[0] = 0
        self._tb.ta[0] = 0
        self._tb.pxa[0] = 0
        self._tb.pya[0] = 0
        self._tb.pza[0] = 0
        self._tb.ekina[0] = 0
        self._tb.ma[0] = 0
        self._tb.edep[0] = 0
        self._tb.ni[0] = 0
        self._tb.nstep[0] = 0
        return

    def EndOfEventAction(self, event):
        """Record event"""

        print("end of event action")
        tb = self._tb
        # Event ID
        tb.ev[0] = event.GetEventID()
        # Ancestor particle (e.g. neutrino)
        if tb.ancestor != None:
            # Log info
            heppart = tb.ancestor
            tb.pida[0] = heppart["pdgid"]
            heppos = heppart["position"]
            hepmom = heppart["momentum"]
            tb.xa[0] = heppos[0] / cm
            tb.ya[0] = heppos[1] / cm
            tb.za[0] = heppos[2] / cm
            tb.ta[0] = heppart["time"] / ns
            tb.pxa[0] = hepmom[0] / GeV
            tb.pya[0] = hepmom[1] / GeV
            tb.pza[0] = hepmom[2] / GeV
            tb.ekina[0] = (heppart["energy"] - heppart["mass"]) / GeV
            tb.ma[0] = heppart["mass"] / GeV
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
        return


# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
    "My Stepping Action"

    def __init__(self, treebuffer, geometry):
        """Constructor"""
        G4UserSteppingAction.__init__(self)
        self._tb = treebuffer
        self._geom = geometry

        print("stepping action")

    def UserSteppingAction(self, step):
        """Collect data for current simulation step"""
        tb = self._tb
        istp = tb.nstep[0]
        if istp >= tb.maxTracks:
            print("Reached maximum tracks:", istp)
            return

        print("stepping action")
        track = step.GetTrack()
        prestep = step.GetPreStepPoint()
        poststep = step.GetPostStepPoint()
        tb.tid[istp] = track.GetTrackID()
        tb.pid[istp] = track.GetDefinition().GetPDGEncoding()
        tb.parid[istp] = track.GetParentID()
        tb.ekin[istp] = prestep.GetKineticEnergy() / GeV
        # Sum simulated energy deposition by volume
        tb.edep[istp] = step.GetTotalEnergyDeposit() / GeV
        # Capture step position
        prepos = prestep.GetPosition()
        postpos = poststep.GetPosition()
        tb.xs[istp] = prepos.x / cm
        tb.ys[istp] = prepos.y / cm
        tb.zs[istp] = prepos.z / cm
        tb.xe[istp] = postpos.x / cm
        tb.ye[istp] = postpos.y / cm
        tb.ze[istp] = postpos.z / cm
        tb.nstep[0] += 1
        return


###################################################################

if "__main__" == __name__:
    # Run the simulation

    mysim = MySimulation()
    mysim.initialize()
    mysim.run()
    mysim.finalize()
