//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
//                                                          
// $Id: MGGeneratorCosmicMessenger.cc
//      
// CLASS IMPLEMENTATION:  MGGeneratorCosmicMessenger.cc
//
//---------------------------------------------------------------------------//
/* AUTHOR: Michael Gold
* CONTACT: mgold@unm.edu
* FIRST SUBMISSION: Mar 4 ,2019
* 
* REVISION:
* SPECIAL NOTES:
* 
*
*/
// 
//---------------------------------------------------------------------------//
/**
 MGGeneratorCosmicMessenger.cc 
  
  MGGeneratorCosmicMessenger.c
 
 c* AUTHOR: Neil McFadden
* CONTACT: nmcfadde@unm.edu
* FIRST SUBMISSION:
* 
* REVISION:
*
* 
*/
//---------------------------------------------------------------------------//
//

#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "generators/MGGeneratorCosmic.hh"
#include "generators/MGGeneratorCosmicMessenger.hh"

//
//---------------------------------------------------------------------------//

MGGeneratorCosmicMessenger::MGGeneratorCosmicMessenger(MGGeneratorCosmic *generator) : fCosmicGenerator(generator)
{
	// /MG/generator/LiquidArgon
	fCosmicDirectory = new G4UIdirectory("/MG/generator/Cosmic/");
	fCosmicDirectory->SetGuidance("Set to generate optical photons @128 nm in argon inside cryostat");

  fCosmicSetScintMean = new G4UIcmdWithADoubleAndUnit("/MG/generator/Cosmic/SetScintMean",this);
  fCosmicSetScintMean->SetGuidance("Default scintillation gaussian mean is 128 nm, can change it with this command");
  fCosmicSetScintMean->SetDefaultUnit("nm");
  fCosmicSetScintMean->SetUnitCategory("Length");
  fCosmicSetScintMean->SetUnitCandidates("nm micron mm cm m km");

  fCosmicSetScintSigma = new G4UIcmdWithADoubleAndUnit("/MG/generator/Cosmic/SetScintSigma",this);
  fCosmicSetScintSigma->SetGuidance("Default scintillation gaussian sigma is 2.929 nm, can change it with this command");
  fCosmicSetScintSigma->SetDefaultUnit("nm");
  fCosmicSetScintSigma->SetUnitCategory("Length");
  fCosmicSetScintSigma->SetUnitCandidates("nm micron mm cm m km");

  fCosmicSetRadius= new G4UIcmdWithADoubleAndUnit("/MG/generator/Cosmic/SetRadius",this);
  fCosmicSetRadius->SetGuidance("Define Max Radius to Generate points inside a Cylidrical Croystat that has diameter = height");
  fCosmicSetRadius->SetDefaultUnit("cm");
  fCosmicSetRadius->SetUnitCategory("Length");
  fCosmicSetRadius->SetUnitCandidates("micron mm cm m km");

  fCosmicSetHeight= new G4UIcmdWithADoubleAndUnit("/MG/generator/Cosmic/SetHeight",this);
  fCosmicSetHeight->SetGuidance("Define Max Height to Generate points inside a Cylidrical Croystat that has diameter = height");
  fCosmicSetHeight->SetDefaultUnit("cm");
  fCosmicSetHeight->SetUnitCategory("Length");
  fCosmicSetHeight->SetUnitCandidates("micron mm cm m km");

  //example
  // /MG/generator/Cosmic/SetCenterVector 0.0 0.0 100.0 cm
  fCosmicSetCenterVector= new G4UIcmdWith3VectorAndUnit("/MG/generator/Cosmic/SetCenterVector",this);
  fCosmicSetCenterVector->SetGuidance("Set Center of generator volume inside a Cylidrical Croystat that has diameter = height");
  fCosmicSetCenterVector->SetGuidance("Default value is (0,0,0) but arrays are not centered on (0,0,0)");
  fCosmicSetCenterVector->SetDefaultUnit("cm");
  fCosmicSetCenterVector->SetUnitCategory("Length");
  fCosmicSetCenterVector->SetUnitCandidates("micron mm cm m km");

  fCosmicSetParticle = new G4UIcmdWithAString("/MG/generator/Cosmic/SetParticle",this);
  fCosmicSetParticle->SetGuidance("Define particle type that will be generated in LAr");

  fCosmicSetEnergy = new G4UIcmdWithADoubleAndUnit("/MG/generator/Cosmic/SetEnergy",this);
  fCosmicSetEnergy->SetGuidance("Define energy to Generate points inside a Cylidrical Croystat that has diameter = height");
  fCosmicSetEnergy->SetDefaultUnit("MeV");
  fCosmicSetEnergy->SetUnitCategory("Energy");
  fCosmicSetEnergy->SetUnitCandidates("eV KeV MeV GeV TeV");
}

//---------------------------------------------------------------------------//

MGGeneratorCosmicMessenger::MGGeneratorCosmicMessenger(const MGGeneratorCosmicMessenger & other) : G4UImessenger(other)
{;}

//---------------------------------------------------------------------------//

MGGeneratorCosmicMessenger::~MGGeneratorCosmicMessenger()
{
	delete fCosmicDirectory;
  delete fCosmicSetRadius;
  delete fCosmicSetHeight;
  delete fCosmicSetParticle;
  delete fCosmicSetEnergy;
  delete fCosmicSetCenterVector;
  delete fCosmicSetScintMean;
  delete fCosmicSetScintSigma;
}

//---------------------------------------------------------------------------//

void MGGeneratorCosmicMessenger::SetNewValue(G4UIcommand *cmd, G4String str)
{
  if(cmd == fCosmicSetRadius){
    fCosmicGenerator->SetRadius(fCosmicSetRadius->GetNewDoubleValue(str));
  }
  else if(cmd == fCosmicSetScintMean){
    fCosmicGenerator->SetMean(fCosmicSetScintMean->GetNewDoubleValue(str));
  }
  else if(cmd == fCosmicSetScintSigma){
    fCosmicGenerator->SetSigma(fCosmicSetScintSigma->GetNewDoubleValue(str));
  }
  else if(cmd == fCosmicSetHeight){
    fCosmicGenerator->SetHeight(fCosmicSetHeight->GetNewDoubleValue(str));
  }
  else if(cmd == fCosmicSetCenterVector){
    fCosmicGenerator->SetCenterVector(fCosmicSetCenterVector->GetNew3VectorValue(str));
  }
  else if(cmd == fCosmicSetParticle){
    fCosmicGenerator->SetParticleType(str);
  }
  else if(cmd == fCosmicSetEnergy){
    fCosmicGenerator->SetParticleEnergy(fCosmicSetEnergy->GetNewDoubleValue(str));
  }
}
