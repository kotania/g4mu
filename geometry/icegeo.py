from geant4_pybind import *


class IceCubeDetectorConstruction(G4VUserDetectorConstruction):
    def __init__(self):
        super().__init__()
        self.fScoringVolume = None

    def Construct(self):
        nist = G4NistManager.Instance()

        # Option to switch on/off checking of volumes overlaps
        checkOverlaps = True

        # World
        world_sizeXY = 2000 * m
        world_sizeZ = 10000 * m

        H = nist.FindOrBuildElement("H")
        O = nist.FindOrBuildElement("O")
        # D = G4Element("Deuterium", "D", Zeff=1.0, Aeff=2.013553 * g / mole)
        # O18 = G4Element("IsoOxygen", "O18", Zeff=8.0, Aeff=17.99916 * g / mole)

        normal_ice = G4Material("normal_ice", 0.9216 * g / cm3, nComponents=2)
        normal_ice.AddElement(H, nAtoms=2)
        normal_ice.AddElement(O, nAtoms=1)

        # semiheavy_ice = G4Material("semiheavy_ice", 0.9712 * g / cm3, nComponents=3)
        # semiheavy_ice.AddElement(H, nAtoms=1)
        # semiheavy_ice.AddElement(D, nAtoms=1)
        # semiheavy_ice.AddElement(O, nAtoms=1)

        # H2O18 = G4Material("H2O18", 1.0254 * g / cm3, nComponents=2)
        # H2O18.AddElement(H, nAtoms=2)
        # H2O18.AddElement(O18, nAtoms=1)

        # iso_ice = G4Material("iso_ice", 0.92 * g / cm3, nComponents=3)
        # iso_ice.AddMaterial(normal_ice, 99.7 * perCent)
        # iso_ice.AddMaterial(semiheavy_ice, 0.03 * perCent)
        # iso_ice.AddMaterial(H2O18, 0.2 * perCent)

        # # We got Air bubbles in the ice!
        # AirBubble = nist.FindOrBuildMaterial("G4_AIR")

        # # The more realistic ice
        # Ice = G4Material("Ice", 0.9216 * g / cm3, nComponents=2)
        # Ice.AddMaterial(iso_ice, 99.9892 * perCent)
        # Ice.AddMaterial(AirBubble, 0.0108 * perCent)
        # Ice.GetIonisation().SetMeanExcitationEnergy(75.0 * eV)

        solidWorld = G4Box(
            "World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ
        )

        logicWorld = G4LogicalVolume(solidWorld, normal_ice, "World")

        physWorld = G4PVPlacement(
            None,  # no rotation
            G4ThreeVector(),  # at (0,0,0)
            logicWorld,  # its logical volume
            "World",  # its name
            None,  # its mother  volume
            False,  # no boolean operation
            0,  # copy number
            checkOverlaps,
        )  # overlaps checking

        # Set logicWorld as scoring volume
        self.fScoringVolume = logicWorld

        # always return the physical World
        return physWorld
