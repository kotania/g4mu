import sys

sys.path.append("../")

from geant4 import *
from geant4.utils import ParticleGun
import geometry

# User Action Initialization
class AppBuilder(G4VUserActionInitialization):
    def Build(self):
        # setup particle gun
        particle_gun = ParticleGun()
        self.SetUserAction(particle_gun)

        # setup particle gun
        pg = particle_gun.GetGun()
        pg.SetParticleByName("mu-")
        pg.SetParticleEnergy(100.0 * GeV)
        pg.SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0))
        pg.SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0) * m)


# set detector construction
geom = geometry.geometry()
gRunManager.SetUserInitialization(geom)


# set physics list
physics_list = QBBC()
gRunManager.SetUserInitialization(physics_list)

# user action initialization
app_builder = AppBuilder()
gRunManager.SetUserInitialization(app_builder)

gApplyUICommand("/analysis/setFileName output_trial")

gRunManager.Initialize()
gRunManager.BeamOn(0)

gTrackingManager.SetVerboseLevel(1)
gRunManager.BeamOn(1)
