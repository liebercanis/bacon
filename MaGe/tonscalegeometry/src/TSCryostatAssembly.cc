//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
//                                                          
// $Id: MGcodeTemplate.cc,v 1.1 2004-12-09 08:58:35 pandola Exp $ 
//      
// CLASS IMPLEMENTATION:  TSCryostatAssembly.cc
//
//---------------------------------------------------------------------------//
/**
 * SPECIAL NOTES:
 *Modified version of MJDemoCryostatAssembly 
 *Assembly origin:  Coaxial with barrel, coincident with crossarm axis
 */
// 
//---------------------------------------------------------------------------//
/**
 * AUTHOR: Matthew Green
 * CONTACT: mpgreen@physics.unc.edu
 * FIRST SUBMISSION: Sept 10, 2010
 * 
 * REVISION:
 * 
 * 09-10-2010, Created, M. Green
 * 03-02-2011, Modified by R. Martin to remove cross arm
 */
//---------------------------------------------------------------------------//
//

#include "G4LogicalVolume.hh"	
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4AffineTransform.hh"

//---------------------------------------------------------------------------//

#include "io/MGLogger.hh"
#include "tonscalegeometry/TSCryostatAssembly.hh"
#include "mjdemonstrator/MJDemoCryostatTopLid.hh"
#include "mjdemonstrator/MJDemoCryostatBottomLid.hh"
#include "tonscalegeometry/TSCryostatHoop.hh"
#include "mjdemonstrator/MJDemoCryostatCrossArmTube.hh"
#include "mjdemonstrator/MJDemoCryostatRailAssembly.hh"

using namespace CLHEP;


//---------------------------------------------------------------------------//

TSCryostatAssembly::TSCryostatAssembly(G4String fName, G4String serialNumber) :
  MJVDemoAssembly(fName, serialNumber, "LMFE")
{
  fTopLidPtr = new MJDemoCryostatTopLid(fName+"_CryostatTopLid_001", "867-5309");  
  fHoopPtr = new TSCryostatHoop(fName+"_CryostatHoop_001", "867-5309");  
  fBottomLidPtr = new MJDemoCryostatBottomLid(fName+"_CryostatBottomLid_001", "867-5309");  
  //fCrossArmTubePtr = new MJDemoCryostatCrossArmTube(fName+"_CrossArmTube_001", "867-5309");
  //fRailAssembly1Ptr = new MJDemoCryostatRailAssembly(fName+"_RailAssembly_001", "867-5309"); 
  //fRailAssembly2Ptr = new MJDemoCryostatRailAssembly(fName+"_RailAssembly_002", "867-5309"); 
}

// -----------------------------------------------------------------------

TSCryostatAssembly::TSCryostatAssembly(const TSCryostatAssembly & rhs) : 
  MJVDemoAssembly(rhs)
{
  fTopLidPtr = new MJDemoCryostatTopLid(fName+"_CryostatTopLid_001", "867-5309");  
  fHoopPtr = new TSCryostatHoop(fName+"_CryostatHoop_001", "867-5309");  
  fBottomLidPtr = new MJDemoCryostatBottomLid(fName+"_CryostatBottomLid_001", "867-5309");  
  // fCrossArmTubePtr = new MJDemoCryostatCrossArmTube(fName+"_CrossArmTube_001", "867-5309");
  //fRailAssembly1Ptr = new MJDemoCryostatRailAssembly(fName+"_RailAssembly_001", "867-5309"); 
  //fRailAssembly2Ptr = new MJDemoCryostatRailAssembly(fName+"_RailAssembly_002", "867-5309"); 
}

// -----------------------------------------------------------------------

TSCryostatAssembly::~TSCryostatAssembly()
{
  delete fTopLidPtr;
  delete fHoopPtr;
  delete fBottomLidPtr;
  //delete fCrossArmTubePtr;
  //delete fRailAssembly1Ptr;
  //delete fRailAssembly2Ptr;
}

// -----------------------------------------------------------------------

void TSCryostatAssembly::Place(G4ThreeVector* assemPosition, G4RotationMatrix* assemRotation, G4LogicalVolume* motherLogical)
{
  G4AffineTransform *assemAffine = new G4AffineTransform(assemRotation,*assemPosition);

  G4ThreeVector *hooplocalPos = new G4ThreeVector(0, 0, 3.0*25.4*mm);
  G4RotationMatrix *hooplocalRot = new G4RotationMatrix();
  G4AffineTransform *hoopaffine = new G4AffineTransform(hooplocalRot,*hooplocalPos);
  *hoopaffine *= *assemAffine;  
  G4ThreeVector *hoopglobalPos = new G4ThreeVector(hoopaffine->NetTranslation());
  G4RotationMatrix *hoopglobalRot= new G4RotationMatrix(hoopaffine->NetRotation());
  fHoopPtr->Place(hoopglobalPos, hoopglobalRot, motherLogical);
  delete hooplocalPos;
  delete hooplocalRot;
  delete hoopaffine;
  
  G4ThreeVector *topLidlocalPos = new G4ThreeVector(0, 0, 2.8483*25.4*mm);//2.8473
  G4RotationMatrix *topLidlocalRot = new G4RotationMatrix();
  G4AffineTransform *topLidaffine = new G4AffineTransform(topLidlocalRot,*topLidlocalPos);
  *topLidaffine *= *assemAffine;  
  G4ThreeVector *topLidglobalPos = new G4ThreeVector(topLidaffine->NetTranslation());
  G4RotationMatrix *topLidglobalRot= new G4RotationMatrix(topLidaffine->NetRotation());
  fTopLidPtr->Place(topLidglobalPos, topLidglobalRot, motherLogical);
  delete topLidlocalPos;
  delete topLidlocalRot;
  delete topLidaffine;
  
  G4ThreeVector *bottomLidlocalPos = new G4ThreeVector(0, 0, -2.8383*25.4*mm);//-2.8373
  G4RotationMatrix *bottomLidlocalRot = new G4RotationMatrix();
  G4AffineTransform *bottomLidaffine = new G4AffineTransform(bottomLidlocalRot,*bottomLidlocalPos);
  *bottomLidaffine *= *assemAffine;  
  G4ThreeVector *bottomLidglobalPos = new G4ThreeVector(bottomLidaffine->NetTranslation());
  G4RotationMatrix *bottomLidglobalRot= new G4RotationMatrix(bottomLidaffine->NetRotation());
  fBottomLidPtr->Place(bottomLidglobalPos, bottomLidglobalRot, motherLogical);
  delete bottomLidlocalPos;
  delete bottomLidlocalRot;
  delete bottomLidaffine;  

  /*
  G4ThreeVector *crossArmTubelocalPos = new G4ThreeVector(6.176*25.4*mm, 0, 0);
  G4RotationMatrix *crossArmTubelocalRot = new G4RotationMatrix();
  crossArmTubelocalRot->rotateY(-pi/2);
  G4AffineTransform *crossArmTubeaffine = new G4AffineTransform(crossArmTubelocalRot,*crossArmTubelocalPos);
  *crossArmTubeaffine*=*assemAffine;  
  G4ThreeVector *crossArmTubeglobalPos = new G4ThreeVector(crossArmTubeaffine->NetTranslation());
  G4RotationMatrix *crossArmTubeglobalRot= new G4RotationMatrix(crossArmTubeaffine->NetRotation());
  fCrossArmTubePtr->Place(crossArmTubeglobalPos, crossArmTubeglobalRot, motherLogical);
  delete crossArmTubelocalPos;
  delete crossArmTubelocalRot;
  delete crossArmTubeaffine;  */

	/*
  G4ThreeVector *RailAssembly1localPos = new G4ThreeVector(0, 0, 2.13*25.4*mm);
  G4RotationMatrix *RailAssembly1localRot = new G4RotationMatrix();
  G4AffineTransform *RailAssembly1affine = new G4AffineTransform(RailAssembly1localRot,*RailAssembly1localPos);
  *RailAssembly1affine*=*assemAffine;  
  G4ThreeVector *RailAssembly1globalPos = new G4ThreeVector(RailAssembly1affine->NetTranslation());
  G4RotationMatrix *RailAssembly1globalRot= new G4RotationMatrix(RailAssembly1affine->NetRotation());
  fRailAssembly1Ptr->Place(RailAssembly1globalPos, RailAssembly1globalRot, motherLogical);
  delete RailAssembly1localPos;
  delete RailAssembly1localRot;
  delete RailAssembly1affine;  

  G4ThreeVector *RailAssembly2localPos = new G4ThreeVector(0, 0, -2.3025*25.4*mm);
  G4RotationMatrix *RailAssembly2localRot = new G4RotationMatrix();
  RailAssembly2localRot->rotateX(-pi);
  G4AffineTransform *RailAssembly2affine = new G4AffineTransform(RailAssembly2localRot,*RailAssembly2localPos);
  *RailAssembly2affine*=*assemAffine;  
  G4ThreeVector *RailAssembly2globalPos = new G4ThreeVector(RailAssembly2affine->NetTranslation());
  G4RotationMatrix *RailAssembly2globalRot= new G4RotationMatrix(RailAssembly2affine->NetRotation());
  fRailAssembly2Ptr->Place(RailAssembly2globalPos, RailAssembly2globalRot, motherLogical);
  delete RailAssembly2localPos;
  delete RailAssembly2localRot;
  delete RailAssembly2affine;  
	*/

  delete assemAffine;
}

// -----------------------------------------------------------------------

