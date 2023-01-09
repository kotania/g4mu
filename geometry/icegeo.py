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
        normal_ice.GetIonisation().SetMeanExcitationEnergy(75.0 * eV)

        mpt_ice = G4MaterialPropertiesTable()

        # Refractive index data from https://atmos.uw.edu/ice_optical_constants/Warren_and_Brandt_2008.pdf
        # (see also https://atmos.uw.edu/ice_optical_constants/)
        # The below energies correspond to the 265-675nm wavelength range
        
        photon_energies = [1.84239837, 1.85369213, 1.86512521, 1.8767002 , 1.88841975,
       1.90028659, 1.91230352, 1.92447339, 1.93679916, 1.94928384,
       1.96193051, 1.97474235, 1.98772262, 2.00087466, 2.01420191,
       2.02770788, 2.0413962 , 2.05527059, 2.06933486, 2.08359295,
       2.09804888, 2.1127068 , 2.12757097, 2.14264579, 2.15793575,
       2.1734455 , 2.18917981, 2.20514359, 2.22134191, 2.23777996,
       2.2544631 , 2.27139688, 2.28858696, 2.30603922, 2.32375969,
       2.34175462, 2.36003043, 2.37859373, 2.39745139, 2.41661044,
       2.43607817, 2.45586211, 2.47597001, 2.49640992, 2.5171901 ,
       2.53831914, 2.55980589, 2.58165952, 2.6038895 , 2.62650563,
       2.64951808, 2.67293734, 2.6967743 , 2.72104024, 2.74574684,
       2.77090621, 2.79653092, 2.822634  , 2.84922897, 2.87632986,
       2.90395125, 2.93210828, 2.96081668, 2.99009281, 3.01995368,
       3.05041699, 3.08150114, 3.11322532, 3.1456095 , 3.17867449,
       3.21244199, 3.24693462, 3.282176  , 3.31819078, 3.35500469,
       3.39264464, 3.43113874, 3.4705164 , 3.51080839, 3.55204694,
       3.59426579, 3.63750032, 3.68178763, 3.72716664, 3.77367824,
       3.82136534, 3.87027309, 3.92044896, 3.97194293, 4.02480761,
       4.07909848, 4.13487405, 4.19219606, 4.25112973, 4.31174399,
       4.37411178, 4.43831031, 4.50442138, 4.57253175, 4.64273351]

        n_phase = [1.3075115 , 1.307623  , 1.307705  , 1.307787  , 1.307869  ,
       1.307951  , 1.3080495 , 1.3081725 , 1.3082955 , 1.308379  ,
       1.308461  , 1.3085645 , 1.3086875 , 1.3088105 , 1.3089335 ,
       1.3090565 , 1.3091795 , 1.3093025 , 1.3094255 , 1.3095485 ,
       1.3096715 , 1.3097945 , 1.3099175 , 1.3100405 , 1.3101635 ,
       1.3102865 , 1.3104095 , 1.3105325 , 1.310674  , 1.310838  ,
       1.311002  , 1.311166  , 1.31133   , 1.3114705 , 1.3115935 ,
       1.311722  , 1.311886  , 1.31205   , 1.3122425 , 1.3124475 ,
       1.312642  , 1.312806  , 1.31297   , 1.3131675 , 1.3133725 ,
       1.3135775 , 1.3137825 , 1.3139875 , 1.3141925 , 1.3143975 ,
       1.314623  , 1.314869  , 1.315115  , 1.315361  , 1.315607  ,
       1.315853  , 1.316099  , 1.3163525 , 1.3166395 , 1.3169265 ,
       1.3172135 , 1.3175005 , 1.3178    , 1.318128  , 1.318456  ,
       1.3188195 , 1.3191885 , 1.3195575 , 1.3199265 , 1.3202955 ,
       1.32076575, 1.32123725, 1.32170875, 1.32218025, 1.32265175,
       1.32312325, 1.32359475, 1.32406625, 1.32453775, 1.325071  ,
       1.325809  , 1.326547  , 1.327285  , 1.328023  , 1.328761  ,
       1.329499  , 1.330237  , 1.330975  , 1.331713  , 1.332451  ,
       1.333189  , 1.333951  , 1.335345  , 1.336739  , 1.338133  ,
       1.339527  , 1.340921  , 1.342315  , 1.343709  , 1.345103]

        ice_ephot = G4doubleVector([energy * eV for energy in photon_energies])
        ice_refr = G4doubleVector(n_phase)
        ice_bins = len(ice_refr)
        mpt_ice.AddProperty("RINDEX", ice_ephot, ice_refr, ice_bins)
        normal_ice.SetMaterialPropertiesTable(mpt_ice)

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
