// ------------------------------------------------------------- 
// Implementation of the TG4Sol_Event class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "TG4Sol_Event.hh"

/********************************************************************/
TG4Sol_Event::TG4Sol_Event(){
    Zero();
}

/********************************************************************/
TG4Sol_Event::~TG4Sol_Event(){
}

/********************************************************************/
void TG4Sol_Event::Zero(){
    // Primary deecay vertex
    WasDecay = 0;
    MotherName = "";
    MotherMass = 0.0;
    MotherTrackID = 0.;
    DecayTime = 0.0;
    //
    InteractionPoint_X = 0.0;
    InteractionPoint_Y = 0.0;
    InteractionPoint_Z = 0.0;
    //
    DecayVertex_X = 0.0;
    DecayVertex_Y = 0.0;
    DecayVertex_Z = 0.0;
    //
    MotherMomentumAtDecay_X = 0.0;
    MotherMomentumAtDecay_Y = 0.0;
    MotherMomentumAtDecay_Z = 0.0;
    //
    DaughterNames.resize(0);
    DaughterMasses.resize(0);
    DaughterCharges.resize(0);
    DaughterMomentums_X.resize(0);
    DaughterMomentums_Y.resize(0);
    DaughterMomentums_Z.resize(0);
    DaughterTrackID.resize(0);
    //
    // beam composition
    BeamNames.resize(0);
    BeamMasses.resize(0);
    BeamCharges.resize(0);
    BeamMomentums_X.resize(0);
    BeamMomentums_Y.resize(0);
    BeamMomentums_Z.resize(0);
    BeamTrackID.resize(0);
}

/********************************************************************/
void TG4Sol_Event::Display(){
    std::cout << "The TG4Sol_Event class\n";
    //---------------------------------------------------------------
    // Primary deecay vertex
    std::cout << "Primary Event composition: \n";
    for(UInt_t i=0;i<BeamNames.size();i++){
	std::cout << " Name: " << BeamNames[i] << "\n";
	std::cout << " Particle Mass: " <<
	    BeamMasses[i] << "\n";
	std::cout << " Particle Charge: " <<
	    BeamCharges[i] << "\n";
	std::cout << " Particle Momentum: (" <<
	    BeamMomentums_X[i] << ";";
	std::cout << BeamMomentums_Y[i] << ";";
	std::cout << BeamMomentums_Z[i] << ")\n";
    }
    if(WasDecay > 0){
	std::cout << "Mother: " << MotherName << "\n";
	std::cout << "Mother mass: " << MotherMass << "\n";
	std::cout << "Decaying time: " << DecayTime << "\n";
	std::cout << "Interaction point at (";
	std::cout << InteractionPoint_X << ";";
	std::cout << InteractionPoint_Y << ";";
	std::cout << InteractionPoint_Z << ")\n";
	std::cout << "Decay Vertex at (";
	std::cout << DecayVertex_X << ";";
	std::cout << DecayVertex_Y << ";";
	std::cout << DecayVertex_Z << ")\n";
	std::cout << "Mother Momentum at decay point (";
	std::cout << MotherMomentumAtDecay_X << ";";
	std::cout << MotherMomentumAtDecay_Y << ";";
	std::cout << MotherMomentumAtDecay_Z << ")\n";
	std::cout << "Daughters:\n";
	for(UInt_t i=0;i<DaughterNames.size();i++){
	    std::cout << " Name: " << DaughterNames[i] << "\n";
	    std::cout << " Daughter Mass: " <<
		DaughterMasses[i] << "\n";
	    std::cout << " Daughter Charge: " <<
		DaughterCharges[i] << "\n";
	    std::cout << " Daughter Momentum: (" <<
		DaughterMomentums_X[i] << ";";
	    std::cout << DaughterMomentums_Y[i] << ";";
	    std::cout << DaughterMomentums_Z[i] << ")\n";
	}
    }
}

/********************************************************************/
ClassImp(TG4Sol_Event)
