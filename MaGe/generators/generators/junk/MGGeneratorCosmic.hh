//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
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
/**                                                            
* $Id:
*      
* CLASS DECLARATION:  MGGeneratorPb210.hh
*
*---------------------------------------------------------------------------//
*
* DESCRIPTION: 
*
*/ 
// Begin description of class here
/** simple paramterized cosmic ray muon spectrun
*/
// End class description
//
/**  
* SPECIAL NOTES:
* Using this class requires GSS.  
*
*/
// 
// --------------------------------------------------------------------------//
/** 
* AUTHOR: Michael Gold
* CONTACT: mgold@unm.edu
* FIRST SUBMISSION: Mar 4 ,2019
* 
* REVISION:
* 
*/
// --------------------------------------------------------------------------//

#ifndef _MGGENERATORCOSMIC_HH
#define _MGGENERATORCOSMIC_HH

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

#include "TF1.h"
#include "TH1F.h"

#include "generators/MGVGenerator.hh"

//---------------------------------------------------------------------------//

class G4Event;
class G4Messenger;
class G4ParticleGun;
class G4Run;

class MGGeneratorCosmic : public MGVGenerator
{

  public:
    MGGeneratorCosmic();

    MGGeneratorCosmic(const MGGeneratorCosmic &);

    ~MGGeneratorCosmic();
    //public interface
    void GeneratePrimaryVertex(G4Event *event);
    void SetParticlePosition(G4ThreeVector pos) { fCurrentPosition = pos;} 

    void DirectionDecider();
    void EnergyDecider();
    void PositionDecider();
    void ParticleDecider();
    G4bool IsInArgon(G4ThreeVector rp);
    
    //Messenger Commands
    void SetRadius(G4double r){fRadius = r;}
    void SetHeight(G4double h){fInnerHeight = h;}
    void SetCenterVector(G4ThreeVector vec){fCenterVector = vec;}
    void SetParticleType(G4String str){fParticleType = str;}
    void SetParticleEnergy(G4double nrg){fEnergy = nrg;}
    void SetMean(G4double mean){fMean = mean;}
    void SetSigma(G4double sigma){fSigma = sigma;}



  private:
    static const G4double LambdaE;
	  G4ParticleGun*	fParticleGun;
    //particle properties
    G4double  fCurrentEnergy; // energy of current particle
    G4ThreeVector fCurrentPosition; // current position of particle
    G4ThreeVector fDirection; // direction of momentum
    G4double fRadius = 0;
    G4double fInnerHeight = 0;
    G4ThreeVector fCenterVector;
    G4String fParticleType = "opticalphoton";
    G4double fEnergy = 0;
    G4double fMean = 0;
    G4double fSigma = 0;

};
#endif
