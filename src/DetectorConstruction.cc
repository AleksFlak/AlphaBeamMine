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
/// \brief Implementation of the DetectorConstruction class


#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"

#include "CommandLineParser.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

#include "RunAction.hh"

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes invisGrey(false, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visYellow(true, G4Colour(1.0, 1.0, 0.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visPurple(true, G4Colour(1.0, 0.0, 1.0));

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  // R = {155 * micrometer, 175 * micrometer, 195 * micrometer, 215 * micrometer, 235 * micrometer, 255 * micrometer, 275 * micrometer, 295 * micrometer, 315 * micrometer, 335 * micrometer};

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *S_Steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

  // Define materials
  G4Material* voxelMaterial = new G4Material("voxelMaterial", 1.089 * g/cm3, 9);
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_H"), 10 * perCent); // 10% hydrogen
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_C"), 19.9 * perCent); // 19.9% carbon
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_N"), 4.2 * perCent);  // 4.2% nitrogen
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_O"), 65 * perCent);   // 65% oxygen
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Na"), 0.2 * perCent); // 0.2% sodium
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_P"), 0.1 * perCent);  // 0.1% phosphorus
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_S"), 0.2 * perCent);  // 0.2% sulfur
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Cl"), 0.3 * perCent); // 0.3% chlorine
  voxelMaterial->AddMaterial(man->FindOrBuildMaterial("G4_K"), 0.1 * perCent);  // 0.1% potassium

  G4Material* skullMaterial = new G4Material("skullMaterial", 1.165 * g/cm3, 11);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_H"), 8.8 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_C"), 39.5 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_N"), 2.6 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_O"), 39.5 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Na"), 0.1 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Mg"), 0.1 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_P"), 2.8 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_S"), 0.2 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Cl"), 0.1 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_K"), 0.1 * perCent);
  skullMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Ca"), 6.2 * perCent);

  G4Material* skinMaterial = new G4Material("skinMaterial", 1.088 * g/cm3, 9);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_H"), 10 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_C"), 19.9 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_N"), 4.2 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_O"), 65 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Na"), 0.2 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_P"), 0.1 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_S"), 0.2 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_Cl"), 0.3 * perCent);
  skinMaterial->AddMaterial(man->FindOrBuildMaterial("G4_K"), 0.1 * perCent);

  G4double nucleusSize = 300 * nm;
  G4double margin = 20 * nm;

  G4Box *solidWorld = new G4Box("world", 100 * mm, 100 * mm, 500 * mm);

  G4Box *solidWater = new G4Box("water", 10 * mm, 10 * mm, 400 * mm);

  G4Box *solidVoxel = new G4Box("voxel", nucleusSize / 2 + margin, nucleusSize / 2 + margin, nucleusSize / 2 + margin);

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    air,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicWater = new G4LogicalVolume(solidWater,
                                                    waterMaterial,
                                                    "water");
  G4PVPlacement *physiWater = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                logicWater,
                                                "water",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  // Define solids for the skull material in front and behind the voxel structure
  G4Box* solidSkullFront = new G4Box("skullFront", 10 * um, 10 * um, 2.5 * um);
  G4Box* solidSkullBack = new G4Box("skullBack", 10 * um, 10 * um, 2.5 * um);

  // Define logical volumes for the skull material in front and behind
  G4LogicalVolume* logicSkullFront = new G4LogicalVolume(solidSkullFront, skullMaterial, "skullFront");
  G4LogicalVolume* logicSkullBack = new G4LogicalVolume(solidSkullBack, skullMaterial, "skullBack");
 
  // Define skin solids for behind and in front of voxel structure 
  G4Box* solidSkinFront = new G4Box("skinFront", 10 * um, 10 * um, 6.9 * um);
  G4Box* solidSkinBack = new G4Box("skinBack", 10 * um, 10 * um, 6.9 * um);

  // Define logical volumes for the skin material in front and behind
  G4LogicalVolume* logicSkinFront = new G4LogicalVolume(solidSkinFront, skinMaterial, "skinFront");
  G4LogicalVolume* logicSkinBack = new G4LogicalVolume(solidSkinBack, skinMaterial, "skinBack");

  // Place the skin materials in front and behind the voxel structure
  G4PVPlacement* physiSkinFront = new G4PVPlacement(0, G4ThreeVector(0, 0, -93.1 * um), logicSkinFront, "skinFront", logicWater, false, 0);
  //G4PVPlacement* physiSkinBack = new G4PVPlacement(0, G4ThreeVector(0, 0, 143.1 * um), logicSkinBack, "skinBack", logicWater, false, 0);

  // Place the skull materials in front and behind the voxel structure
  G4PVPlacement* physiSkullFront = new G4PVPlacement(0, G4ThreeVector(0, 0, -49.1 * um), logicSkullFront, "skullFront", logicWater, false, 0);
  //G4PVPlacement* physiSkullBack = new G4PVPlacement(0, G4ThreeVector(0, 0, 99.1 * um), logicSkullBack, "skullBack", logicWater, false, 0);

  G4LogicalVolume *logicVoxel = new G4LogicalVolume(solidVoxel,
                                                    voxelMaterial, // Use the custom material
                                                    "voxel");
  G4UserLimits *userLimits = new G4UserLimits();
  userLimits->SetMaxAllowedStep(3 * nm);
  G4int noVoxels = 0;
  G4double spacing = 0.5;

  // displacent factor to set the origin just before the skinFront object
  G4ThreeVector displacement = -G4ThreeVector(0, 0, -110 * um);

  for (G4int i = -5; i < 5; i++)
  {
    for (G4int j = -5; j < 5; j++)
    {
      for (G4int k = 0; k < 15000; k++) //vary the second number to vary the voxelised region size
      {

        // G4RotationMatrix* rot = new G4RotationMatrix(-1*theta,
        //                                   0,
        //                                   0) ;

        G4PVPlacement *physiCell = new G4PVPlacement(0,
                                                     G4ThreeVector(0.5 * um + i * spacing * um, 0.5 * um + j * spacing * um, 5 * um + k * spacing * um),
                                                     logicVoxel,
                                                     "voxel",
                                                     logicWater,
                                                     0,
                                                     k,
                                                     0);
        noVoxels++;
        physiCell->SetTranslation(physiCell->GetTranslation() + displacement); //applying displacement to voxelised region
      }
    }
  }

  // applying displacement to all the objects
  physiSkinFront->SetTranslation(physiSkinFront->GetTranslation() + displacement);
  //physiSkinBack->SetTranslation(physiSkinBack->GetTranslation() + displacement);
  physiSkullFront->SetTranslation(physiSkullFront->GetTranslation() + displacement);
  //physiSkullBack->SetTranslation(physiSkullBack->GetTranslation() + displacement);

  //applying colour to objects to distuingish in gui
  G4cout << "placed " << noVoxels << " voxels. " << G4endl;
  logicWater->SetVisAttributes(&invisGrey);
  //logicSkullBack->SetVisAttributes(&visRed);
  logicSkullFront->SetVisAttributes(&visRed);
  //logicSkinBack->SetVisAttributes(&visYellow);
  logicSkinFront->SetVisAttributes(&visYellow);
  logicVoxel->SetVisAttributes(&visBlue);
  logicWorld->SetVisAttributes(&invisGrey);

  return physiWorld;
}