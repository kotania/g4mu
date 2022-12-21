//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace IceGeo
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int ncomponents, natoms;
    
  //***Elements
  G4Element* H = new G4Element("Hydrogen",   "H",  z=1.,  a=  1.00794*g/mole);
  G4Element* C = new G4Element("C",          "C",  z=6.,  a= 12.01   *g/mole);
  G4Element* N = new G4Element("Nitrogen",   "N",  z=7.,  a= 14.0067 *g/mole);
  G4Element* O = new G4Element("Oxygen",     "O",  z=8.,  a= 15.9994 *g/mole);
    
  G4Element* D = new G4Element("Deuterium","D",  z=1.,  a=  2.013553*g/mole);
  G4Element* O18 = new G4Element("IsoOxygen","O18",  z=8.,  a= 17.99916*g/mole);


  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();


  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 2000*m;
  G4double world_sizeZ  = 2000*m;

  // Ice material definition from https://github.com/claudiok/clsim/blob/master/private/geant4/TrkDetectorConstruction.cxx

  G4Material* normal_ice = new G4Material("normal_ice", 0.9216*g/cm3 , ncomponents=2);
  normal_ice->AddElement(H, natoms=2);
  normal_ice->AddElement(O, natoms=1);

  G4Material* semiheavy_ice = new G4Material("semiheavy_ice", 0.9712*g/cm3 , ncomponents=3);
  semiheavy_ice->AddElement(H, natoms=1);
  semiheavy_ice->AddElement(D, natoms=1);
  semiheavy_ice->AddElement(O, natoms=1);

  G4Material* H2O18 = new G4Material("H2O18", 1.0254*g/cm3 , ncomponents=2);
  H2O18->AddElement(H, natoms=2);
  H2O18->AddElement(O18, natoms=1);

  G4Material* iso_ice = new G4Material("iso_ice", 0.92*g/cm3, ncomponents=3);
  iso_ice->AddMaterial(normal_ice, 99.7*perCent);
  iso_ice->AddMaterial(semiheavy_ice, 0.03*perCent);
  iso_ice->AddMaterial(H2O18, 0.2*perCent);

  // We got Air bubbles in the ice!
  G4NistManager* man = G4NistManager::Instance();
  G4Material* AirBubble  = man->FindOrBuildMaterial("G4_AIR");

  // The more realistic ice 
  G4Material* Ice = new G4Material("Ice", 0.9216*g/cm3, ncomponents=2);
  Ice->AddMaterial(iso_ice,  99.9892*perCent);
  Ice->AddMaterial(AirBubble, 0.0108*perCent);
  Ice->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        Ice,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  
  // Set logicWorld as scoring volume
  //
  fScoringVolume = logicWorld;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
