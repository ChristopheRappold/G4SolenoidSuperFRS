#include "TG4Sol_Event.hh"
#include "TG4Sol_Hit.hh"

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TH2F.h"
#include "TVector3.h"
#include "TCanvas.h"

#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <iostream>
#include <tuple>

#include "ROOT/TSeq.hxx"

//#include "ROOT/TThreadedObject.hxx"
//#include "ROOT/TProcessExecutor.hxx"
//#include "TTreeProcessorMT.hxx"

void runData(const std::string& nameList, const std::set<std::string>& ParticleList = std::set<std::string>(), const std::string& nameOut = "" )
{
  //int nthreads = 4;
  //ROOT::EnableImplicitMT(
  
  TChain* Chain = new TChain("G4Tree");
  TFile* f_first = nullptr;
  std::cout << "Files :" << nameList << std::endl;
  if(nameList.find(".root") != std::string::npos)
    {
      std::cout << "Load from single file " << nameList << std::endl;
      int temp_nb = Chain->AddFile(nameList.c_str());
      f_first = new TFile(nameList.c_str());
      std::cout << " Loaded " << temp_nb << " files " << std::endl;

    }
  else
    {
      std::cout << "Adding Chain from List files" << std::endl;
      std::ifstream List(nameList.c_str());
      std::string infiles;
      int nb_file = 0;
      while(std::getline(List, infiles))
        {	  
	  std::cout << infiles << std::endl;
	  int temp_nb = Chain->AddFile(infiles.c_str());
          if(nb_file==0)
	    f_first = new TFile(infiles.c_str());
	  ++nb_file;
	}
      std::cout << " Loaded " << nb_file << " files " << std::endl;
    }

  //TTree* tt = (TTree*) f->Get("G4Tree");
  std::vector<std::string>* nameDetInFile = (std::vector<std::string>*)(f_first->Get("nameDet"));
  std::vector<std::string> nameDet ;
  for(auto name : *nameDetInFile)
    nameDet.emplace_back(name);
  f_first->Close();
  //f_first->Delete();
  std::cout<<" load nameDet done !\n";
  for(auto name : nameDet)
    std::cout<<name<<"\n";

  size_t id_Si3 = 0;
  size_t id_FMF2 = 0;
  size_t id_TrckFwd = 0; 
  for(size_t id = 0; id < nameDet.size(); ++id)
    {
      if(nameDet[id]=="HypHI_InSi_log3")
	id_Si3 = id;
      if(nameDet[id]=="FMF2_log")
	id_FMF2 = id;
      if(nameDet[id] == "HypHI_TrackFwd_log")
	id_TrckFwd = id;
    }

    
  //ROOT::TTreeProcessorMT tp(Chain);
  TTreeReader reader(Chain);
  
  //ROOT::TThreadedObject<TH2F> h_HitPatternXAll("h_HitPatternX","h_HitPatternX",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_HitPatternYAll("h_HitPatternY","h_HitPatternY",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_ParticlePhiAll("h_Phi","h_Phi",360*5,-180,180,30,0,30);
 
  auto f_Process = [&](TTreeReader& reader)
		   {
		     
		     TTreeReaderValue<TG4Sol_Event> Revent(reader,"TG4Sol_Event");
		     std::vector<TTreeReaderArray<TG4Sol_Hit>*> AllHits;
		     for(auto name : nameDet)
		       AllHits.emplace_back(new TTreeReaderArray<TG4Sol_Hit>(reader,name.c_str()));
		     
		     TH2F* h_HitPatternX = new TH2F("h_HitPatternX","h_HitPatternX",2000,-20,20,20,0,20);
		     TH2F* h_HitPatternY = new TH2F("h_HitPatternY","h_HitPatternY",2000,-20,20,20,0,20);
		     TH2F* h_ParticlePhi = new TH2F("h_Phi","h_Phi",360*5,-180,180,30,0,30);
  
		     // auto h_HitPatternX = h_HitPatternX.Get();
		     // auto h_HitPatternY = h_HitPatternY.Get();
		     // auto h_ParticlePhi = h_ParticlePhiAll.Get();
		     
		     const auto Entries = Chain->GetEntries();
		     std::cout << " Entries :" << Entries << std::endl;
		     int timing = 0;
		     int first_event = 0;
		     while(reader.Next())
		       {
			 int nb = reader.GetCurrentEntry();
			 auto event = Revent.Get();
			 if(static_cast<int>(static_cast<double>(nb) / static_cast<double>(Entries) * 10) == timing)
			   {
			     std::cout <<" Processing :" << timing * 10 << "%  Event #"<<nb<<" \n";
			     ++timing;
			   }
			 
			 std::set<int> validTrack;
			 
			 for(size_t id = 0;id<event->BeamNames.size();++id)
			   {
			     int trackID = event->BeamTrackID[id];
			     //std::cout<<"beam : "<<event->BeamNames[id]<<" #"<<trackID<<"\n";
			     if(ParticleList.size()>0)
			       {
				 auto it_par = ParticleList.find(event->BeamNames[id]);
				 if(it_par != ParticleList.end())
				   validTrack.insert(trackID);
			       }
			     else
			       validTrack.insert(trackID);
			     
			   }
			 for(size_t id = 0;id<event->DaughterNames.size();++id)
			   {
			     int trackID = event->DaughterTrackID[id];
			     if(ParticleList.size()>0)
			       {
				 auto it_par = ParticleList.find(event->DaughterNames[id]);
				 if(it_par != ParticleList.end())
				   validTrack.insert(trackID);
			       }
			     else
			       validTrack.insert(trackID);
			   }
			 
			 for(size_t it_br = 0 ;it_br < AllHits.size(); ++it_br)
			   {
			     if(first_event==0)
			       {
				 TString nameBranch = AllHits[it_br]->GetBranchName();
				 if(it_br == id_FMF2 || it_br == id_TrckFwd)
				   {
				     for(auto index : ROOT::TSeqI(3)) 
				       {
					 TString nameBranch2 = nameBranch;
					 nameBranch2 += "_";
					 nameBranch2 += index;
					 std::cout<<"name "<<it_br<<" : "<<nameBranch2<<"\n";
					 h_HitPatternX->Fill(h_HitPatternX->GetXaxis()->GetXmin()-10,nameBranch2,1.);
					 h_HitPatternY->Fill(h_HitPatternY->GetXaxis()->GetXmin()-10,nameBranch2,1.);
				       }
				   }
				 else
				   {
				     std::cout<<"name "<<it_br<<" : "<<nameBranch<<"\n";
				     h_HitPatternX->Fill(h_HitPatternX->GetXaxis()->GetXmin()-10,nameBranch,1.);
				     h_HitPatternY->Fill(h_HitPatternY->GetXaxis()->GetXmin()-10,nameBranch,1.);
				   }
			       }
			     
			     for(auto hit : *AllHits[it_br])
			       {
				 //hit.Print();
				 auto it_find = validTrack.find(hit.TrackID);
				 
				 if(it_find != validTrack.end())
				   {
				     //std::cout<<"Branch: "<<Hits->GetBranchName()<<" "<<hit.HitPosX<<" "<<hit.HitPosY<<" "<<hit.HitPosZ<<" "<<" "<<hit.LayerID<<" "<<hit.Pname<<"\n";
				     TString nameBranch = AllHits[it_br]->GetBranchName();
				     if(it_br == id_FMF2 || it_br == id_TrckFwd)
				       {
					 nameBranch += "_";
					 nameBranch += hit.LayerID;
				       }
				     h_HitPatternX->Fill(hit.HitPosX,nameBranch,1.);
				     h_HitPatternY->Fill(hit.HitPosY,nameBranch,1.);
				   }
			       }
			   }
			 if(first_event==0)
			   ++first_event;

			 std::unordered_map<int,double> PhiPerTrack;
			 
			 for(auto hit : *AllHits[id_Si3])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     
			     if(it_find != validTrack.end())
			       {
				 TVector3 TempMom(hit.MomX,hit.MomY,hit.MomZ);
				 PhiPerTrack.insert(std::make_pair(hit.TrackID,TempMom.Phi()));
			       }
			   }
			 for(auto hit : *AllHits[id_FMF2])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     
			     if(it_find != validTrack.end())
			       {
				 if(hit.LayerID == 0)
				   {
				     TVector3 TempMom(hit.MomX,hit.MomY,hit.MomZ);
				     auto it_phi = PhiPerTrack.find(hit.TrackID);
				     if(it_phi != PhiPerTrack.end())
				       {
					 h_ParticlePhi->Fill((it_phi->second-TempMom.Phi())*TMath::RadToDeg(),hit.Pname.c_str(),1.);
				       }
				   }
			       }
			   }
	  
			 // for(size_t id = 0;id<event->DaughterNames.size();++id)
			 // 	{
			 // 	  int trackID = event->DaughterTrackID[id];
			 // 	  std::cout<<"decayed : "<<event->DaughterNames[id]<<" #"<<trackID<<"\n";
			 
			 // 	  for(auto& Hits :AllHits)
			 // 	    {
			 // 	      for(auto hit : *Hits)
			 // 		{
			 // 		  //hit.Print();
			 // 		  if(hit.TrackID==trackID)
			 // 		    {
			 // 		      std::cout<<"Branch: "<<Hits->GetBranchName()<<" "<<hit.HitPosX<<" "<<hit.HitPosY<<" "<<hit.HitPosZ<<" "<<hit.LayerID<<" "<<hit.Pname<<"\n";
			 // 		    }
			 // 		}
			 // 	    }
			 // 	}
			 

      
		       }

		     //TObjArray* obj_ret = new TObjArray(3);
		     //obj_ret->Add(h_HitPatternX);
		     //obj_ret->Add(h_HitPatternY);
		     //obj_ret->Add(h_ParticlePhi);
		     
		     //return h_HitPatternX;//obj_ret;
		     return std::make_tuple(h_HitPatternX,h_HitPatternY,h_ParticlePhi);
		   };			    
					    

  // Create the pool of processes
  //ROOT::TProcessExecutor workers(4);

  // Process the TChain
  //auto ObjArray = workers.ProcTree(*Chain, f_Process, "G4Tree");
  //ObjArray->Print();
  //tp.Process(f_Process);

  TH2F* h_HitPatternXM = nullptr;
  TH2F* h_HitPatternYM = nullptr;
  TH2F* h_ParticlePhiM = nullptr;
  std::tie(h_HitPatternXM,h_HitPatternYM,h_ParticlePhiM) = f_Process(reader);


  //TH2F* h_HitPatternXM = static_cast<TH2F*>(ObjArray);
  //TH2F* h_HitPatternYM = static_cast<TH2F*>(ObjArray->At(1));
  //TH2F* h_ParticlePhiM = static_cast<TH2F*>(ObjArray->At(2));
  
  auto plotHist = [](auto h, double min=-1, double max=-1, TString nameX="")
		  {
		    TString nameC(h->GetName());
		    nameC+="Canvas";
		    TCanvas* c = new TCanvas(nameC,nameC,600,600);
		    c->cd();
		    h->Draw();

		    // TString nameY("Counts / ");
		    // double BinW = h->GetXaxis()->GetBinWidth(1);
		    // BinW *= 1.e3;
		    // int BinW2 = BinW*10;
		    // BinW = BinW2/10.;
		    // nameY+=TString::Format("%.1f",BinW);
		    // nameY+=" MeV";
		    // TString FullTitle(";");
		    // FullTitle+=nameX;
		    // FullTitle+=";";
		    // FullTitle+=nameY;
		    // FullTitle+=";";
		    // h->SetTitle(FullTitle);
		    
		    double min1 = min < 0 ? h->GetXaxis()->GetXmin() : min;
		    double max1 = max < 0 ? h->GetXaxis()->GetXmax() : max;
		    h->GetXaxis()->SetRangeUser(min,max);
		    h->Draw();
		    c->Draw("colz");
		  };




  
  //auto h_HitPatternXM = h_HitPatternX.Merge(); 
  //auto h_HitPatternYM = h_HitPatternY.Merge(); 
  //auto h_ParticlePhiM = h_ParticlePhi.Merge(); 

  if(nameOut == "")
    {
      plotHist(h_HitPatternXM);
      plotHist(h_HitPatternYM);
      plotHist(h_ParticlePhiM);
    }
  else
    {
      std::cout << "Saving: "<<nameOut<<"\n";
      TFile* f_out = new TFile(nameOut.c_str(),"RECREATE");
      f_out->cd();
      h_HitPatternXM->Write();
      h_HitPatternYM->Write();
      h_ParticlePhiM->Write();
      f_out->Close();
    }
  std::cout<<"Done !\n";
}

void runDataCoincidence(const std::string& nameList, const std::set<std::string>& ParticleList = std::set<std::string>(), const std::string& nameOut = "" )
{

  //int nthreads = 4;
  //ROOT::EnableImplicitMT(nthreads);
  
  TChain* Chain = new TChain("G4Tree");
  TFile* f_first = nullptr;
  std::cout << "Files :" << nameList << std::endl;
  if(nameList.find(".root") != std::string::npos)
    {
      std::cout << "Load from single file " << nameList << std::endl;
      int temp_nb = Chain->AddFile(nameList.c_str());
      f_first = new TFile(nameList.c_str());
      std::cout << " Loaded " << temp_nb << " files " << std::endl;

    }
  else
    {
      std::cout << "Adding Chain from List files" << std::endl;
      std::ifstream List(nameList.c_str());
      std::string infiles;
      int nb_file = 0;
      while(std::getline(List, infiles))
        {	  
	  std::cout << infiles << std::endl;
	  int temp_nb = Chain->AddFile(infiles.c_str());
          if(nb_file==0)
	    f_first = new TFile(infiles.c_str());
	  ++nb_file;
	}
      std::cout << " Loaded " << nb_file << " files " << std::endl;
    }

  //TTree* tt = (TTree*) f->Get("G4Tree");
  std::vector<std::string>* nameDetInFile = (std::vector<std::string>*)(f_first->Get("nameDet"));
  std::vector<std::string> nameDet ;
  for(auto name : *nameDetInFile)
    nameDet.emplace_back(name);
  f_first->Close();
  //f_first->Delete();
  std::cout<<" load nameDet done !\n";
  for(auto name : nameDet)
    std::cout<<name<<"\n";

  size_t id_Si3 = 0;
  size_t id_FMF2 = 0;
  size_t id_TrckFwd = 0; 
  size_t id_PSCE = 0;
  size_t id_PSFE = 0;
  std::set<size_t> id_MDCs;
  for(size_t id = 0; id < nameDet.size(); ++id)
    {
      if(nameDet[id]=="HypHI_InSi_log3")
	id_Si3 = id;
      if(nameDet[id]=="FMF2_log")
	id_FMF2 = id;
      if(nameDet[id] == "HypHI_TrackFwd_log")
	id_TrckFwd = id;
      if(nameDet[id] == "PSCE")
	id_PSCE = id;
      if(nameDet[id] == "PSBE")
	id_PSFE = id;
      if(auto found = nameDet[id].find("MG") ;  found != std::string::npos)
	id_MDCs.insert(id);
    }

    
  //ROOT::TTreeProcessorMT tp(Chain);
  TTreeReader reader(Chain);
  
  //ROOT::TThreadedObject<TH2F> h_HitPatternXAll("h_HitPatternX","h_HitPatternX",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_HitPatternYAll("h_HitPatternY","h_HitPatternY",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_ParticlePhiAll("h_Phi","h_Phi",360*5,-180,180,30,0,30);
 
  auto f_Process = [&](TTreeReader& reader)
		   {
		     
		     TTreeReaderValue<TG4Sol_Event> Revent(reader,"TG4Sol_Event");
		     std::vector<TTreeReaderArray<TG4Sol_Hit>*> AllHits;
		     for(auto name : nameDet)
		       AllHits.emplace_back(new TTreeReaderArray<TG4Sol_Hit>(reader,name.c_str()));

		     TH1F* h_particleStatus = new TH1F("h_ParticleStatus","h_ParticleStatus",30,0,30);

		     TH2F* h_HitPatternX = new TH2F("h_HitPatternX","h_HitPatternX",2000,-20,20,20,0,20);
		     TH2F* h_HitPatternY = new TH2F("h_HitPatternY","h_HitPatternY",2000,-20,20,20,0,20);
		     
		     TH2F* h_HitPatternX_PSCE = new TH2F("h_HitPatternX_PSCE","h_HitPatternX_PSCE",2000,-20,20,20,0,20);
		     TH2F* h_HitPatternX_PSFE = new TH2F("h_HitPatternX_PSFE","h_HitPatternX_PSFE",2000,-20,20,20,0,20);
		     TH2F* h_HitPatternY_PSCE = new TH2F("h_HitPatternY_PSCE","h_HitPatternY_PSCE",2000,-20,20,20,0,20);
		     TH2F* h_HitPatternY_PSFE = new TH2F("h_HitPatternY_PSFE","h_HitPatternY_PSFE",2000,-20,20,20,0,20);

		     TH2F* h_ParticlePhi = new TH2F("h_Phi","h_Phi",360*5,-180,180,30,0,30);

		     TH2F* h_CoincidenceMDC_TOF =  new TH2F("h_MDC_TOF","h_MDC_TOF",100,0,100,20,0,20);

		     TH2F* h_R_Z_MDC = new TH2F("h_R_Z_MDC","h_R_Z_MDC",500,0.,2.,400.,0.,0.4);
		     TH2F* h_X_Z_MDC = new TH2F("h_X_Z_MDC","h_X_Z_MDC",500,0.,2.,400.,-0.4,0.4);
		     TH2F* h_Y_Z_MDC = new TH2F("h_Y_Z_MDC","h_Y_Z_MDC",500,0.,2.,400.,-0.4,0.4);
		     // auto h_HitPatternX = h_HitPatternX.Get();
		     // auto h_HitPatternY = h_HitPatternY.Get();
		     // auto h_ParticlePhi = h_ParticlePhiAll.Get();
		     
		     const auto Entries = Chain->GetEntries();
		     std::cout << " Entries :" << Entries << std::endl;
		     int timing = 0;
		     int first_event = 0;
		     int view_event = 1;
		     while(reader.Next())
		       {
			 int nb = reader.GetCurrentEntry();
			 auto event = Revent.Get();
			 if(static_cast<int>(static_cast<double>(nb) / static_cast<double>(Entries) * 10) == timing)
			   {
			     std::cout <<" Processing :" << timing * 10 << "%  Event #"<<nb<<" \n";
			     ++timing;
			   }
			 
			 std::set<int> validTrack;
			 std::unordered_map<int, std::string> nameTrack;
			 
			 for(size_t id = 0;id<event->BeamNames.size();++id)
			   {
			     int trackID = event->BeamTrackID[id];
			     std::string nameP = event->BeamNames[id];
			     if(view_event==nb)
			       std::cout<<"beam : "<<event->BeamNames[id]<<" #"<<trackID<<"\n";
			     if(ParticleList.size()>0)
			       {
				 auto it_par = ParticleList.find(event->BeamNames[id]);
				 if(it_par != ParticleList.end())
				   {
				     validTrack.insert(trackID);
				     nameTrack.insert(std::make_pair(trackID,nameP));
				     h_particleStatus->Fill(nameP.c_str(),1.);
				   }
			       }
			     else
			       {
				 validTrack.insert(trackID);
				 nameTrack.insert(std::make_pair(trackID,nameP));
				 h_particleStatus->Fill(nameP.c_str(),1.);
			       }
			   }
			 for(size_t id = 0;id<event->DaughterNames.size();++id)
			   {
			     int trackID = event->DaughterTrackID[id];
			     std::string nameP = event->DaughterNames[id];
			     nameP += "Decay";
			     if(view_event==nb)
			       std::cout<<"beam : "<<event->DaughterNames[id]<<" #"<<trackID<<"\n";
			     if(ParticleList.size()>0)
			       {
				 auto it_par = ParticleList.find(event->DaughterNames[id]);
				 if(it_par != ParticleList.end())
				   {
				     validTrack.insert(trackID);
				     nameTrack.insert(std::make_pair(trackID,nameP));
				   }
			       }
			     else
			       {
				 validTrack.insert(trackID);
				 nameTrack.insert(std::make_pair(trackID,nameP));
			       }
			   }

			 std::unordered_map<int, std::tuple<double,double,double> > trackOnPSCE, trackOnPSFE;
			 for(auto hit : *AllHits[id_PSCE])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     if(it_find != validTrack.end())
			       {
				 trackOnPSCE.insert(std::make_pair(hit.TrackID,std::make_tuple(hit.HitPosX,hit.HitPosY,hit.HitPosZ)));
			       }
			   }

			 for(auto hit : *AllHits[id_PSFE])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     if(it_find != validTrack.end())
			       {
				 trackOnPSFE.insert(std::make_pair(hit.TrackID,std::make_tuple(hit.HitPosX,hit.HitPosY,hit.HitPosZ)));
			       }
			   }
			 
			 std::unordered_map<int, std::vector<std::tuple<int,double,double,double> > > TracksInMDC;
			 
			 for(size_t it_br = 0 ;it_br < AllHits.size(); ++it_br)
			   {
			     if(first_event==0)
			       {
				 TString nameBranch = AllHits[it_br]->GetBranchName();
				 if(it_br == id_FMF2 || it_br == id_TrckFwd)
				   {
				     for(auto index : ROOT::TSeqI(3)) 
				       {
					 TString nameBranch2 = nameBranch;
					 nameBranch2 += "_";
					 nameBranch2 += index;
					 std::cout<<"name "<<it_br<<" : "<<nameBranch2<<"\n";
					 h_HitPatternX->Fill(h_HitPatternX->GetXaxis()->GetXmin()-10,nameBranch2,1.);
					 h_HitPatternY->Fill(h_HitPatternY->GetXaxis()->GetXmin()-10,nameBranch2,1.);
				       }
				   }
				 else
				   {
				     std::cout<<"name "<<it_br<<" : "<<nameBranch<<"\n";
				     h_HitPatternX->Fill(h_HitPatternX->GetXaxis()->GetXmin()-10,nameBranch,1.);
				     h_HitPatternY->Fill(h_HitPatternY->GetXaxis()->GetXmin()-10,nameBranch,1.);
				   }
			       }
			     
			     auto it_MDC = id_MDCs.find(it_br);
			     
			     for(auto hit : *AllHits[it_br])
			       {
				 //hit.Print();
				 auto it_find = validTrack.find(hit.TrackID);
				 
				 if(it_find != validTrack.end())
				   {
				     //std::cout<<"Branch: "<<Hits->GetBranchName()<<" "<<hit.HitPosX<<" "<<hit.HitPosY<<" "<<hit.HitPosZ<<" "<<" "<<hit.LayerID<<" "<<hit.Pname<<"\n";
				     TString nameBranch = AllHits[it_br]->GetBranchName();
				     if(it_br == id_FMF2 || it_br == id_TrckFwd)
				       {
					 nameBranch += "_";
					 nameBranch += hit.LayerID;
				       }
				     h_HitPatternX->Fill(hit.HitPosX,nameBranch,1.);
				     h_HitPatternY->Fill(hit.HitPosY,nameBranch,1.);

				     if(it_MDC != id_MDCs.end())
				       {
					 auto it_trackInMDC = TracksInMDC.find(hit.TrackID);
					 if(it_trackInMDC == TracksInMDC.end())
					   TracksInMDC.insert(std::make_pair(hit.TrackID, std::vector<std::tuple<int,double,double,double>>(1,std::make_tuple(*it_MDC,hit.HitPosX,hit.HitPosY,hit.HitPosZ))));
					 else
					   it_trackInMDC->second.emplace_back(std::make_tuple(*it_MDC,hit.HitPosX,hit.HitPosY,hit.HitPosZ));
				       }
									   
				     if(auto it_findPSCE = trackOnPSCE.find(hit.TrackID) ; it_findPSCE != trackOnPSCE.end())
				       {
					 h_HitPatternX_PSCE->Fill(hit.HitPosX,nameBranch,1.);
					 h_HitPatternY_PSCE->Fill(hit.HitPosY,nameBranch,1.);
				       }
				     if(auto it_findPSFE = trackOnPSFE.find(hit.TrackID) ; it_findPSFE != trackOnPSFE.end())
				       {
					 h_HitPatternX_PSFE->Fill(hit.HitPosX,nameBranch,1.);
					 h_HitPatternY_PSFE->Fill(hit.HitPosY,nameBranch,1.);
				       }
				   }
			       }
			   }

			 for(auto [idTrack, InMDC] : TracksInMDC )
			   {
			     std::string tempName (nameTrack[idTrack].c_str());
			     h_CoincidenceMDC_TOF->Fill( InMDC.size(),tempName.c_str(), 1.);
			     auto f_R = [](double X,double Y) {
			       return TMath::Sqrt(X*X +Y*Y);
			     };
			     if(view_event==nb)
			       {
				 for(auto [idLayer, LayerX, LayerY, LayerZ] : InMDC)
				   {
				     std::cout<<" MDC idTrack"<< idTrack <<" : "<<idLayer<<" "<<LayerX/100.<<" "<<LayerY/100.<<" "<<LayerZ/100.<<"\n";
				     h_R_Z_MDC->Fill(LayerZ/100.,f_R(LayerX/100.,LayerY/100.));
				     h_X_Z_MDC->Fill(LayerZ/100.,LayerX/100.);
				     h_Y_Z_MDC->Fill(LayerZ/100.,LayerY/100.);
				   }
			       }
			     if(auto it_TOF = trackOnPSFE.find(idTrack) ; it_TOF!=trackOnPSFE.end())
			       {
				 std::string tempName2 (tempName);
				 tempName2 += "_onPSFE";
				 h_CoincidenceMDC_TOF->Fill( InMDC.size(),tempName2.c_str(), 1.);
			       }
			     if(auto it_TOF = trackOnPSCE.find(idTrack) ; it_TOF!=trackOnPSCE.end())
			       {
				 std::string tempName2 (tempName);
				 tempName2 += "_onPSCE";
				 h_CoincidenceMDC_TOF->Fill( InMDC.size(),tempName2.c_str(), 1.);
				 if(view_event==nb)
				   {
				     std::cout<<"PSCE :"<<it_TOF->first<<" "<<std::get<0>(it_TOF->second)/100.<<" "<<std::get<1>(it_TOF->second)/100.<<" "<<std::get<2>(it_TOF->second)/100.<<"\n";
				     h_R_Z_MDC->Fill(std::get<2>(it_TOF->second)/100., f_R(std::get<0>(it_TOF->second)/100.,std::get<1>(it_TOF->second)/100.) );
				     h_X_Z_MDC->Fill(std::get<2>(it_TOF->second)/100., std::get<0>(it_TOF->second)/100.);
				     h_Y_Z_MDC->Fill(std::get<2>(it_TOF->second)/100., std::get<1>(it_TOF->second)/100.);
				   }
			       }
			   }
			 std::unordered_map<int,double> PhiPerTrack;
			 
			 for(auto hit : *AllHits[id_Si3])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     
			     if(it_find != validTrack.end())
			       {
				 TVector3 TempMom(hit.MomX,hit.MomY,hit.MomZ);
				 PhiPerTrack.insert(std::make_pair(hit.TrackID,TempMom.Phi()));
			       }
			   }
			 for(auto hit : *AllHits[id_FMF2])
			   {
			     auto it_find = validTrack.find(hit.TrackID);
			     
			     if(it_find != validTrack.end())
			       {
				 if(hit.LayerID == 0)
				   {
				     TVector3 TempMom(hit.MomX,hit.MomY,hit.MomZ);
				     auto it_phi = PhiPerTrack.find(hit.TrackID);
				     if(it_phi != PhiPerTrack.end())
				       {
					 h_ParticlePhi->Fill((it_phi->second-TempMom.Phi())*TMath::RadToDeg(),hit.Pname.c_str(),1.);
				       }
				   }
			       }
			   }

			 if(first_event==0)
			   ++first_event;

			 // for(size_t id = 0;id<event->DaughterNames.size();++id)
			 // 	{
			 // 	  int trackID = event->DaughterTrackID[id];
			 // 	  std::cout<<"decayed : "<<event->DaughterNames[id]<<" #"<<trackID<<"\n";
			 
			 // 	  for(auto& Hits :AllHits)
			 // 	    {
			 // 	      for(auto hit : *Hits)
			 // 		{
			 // 		  //hit.Print();
			 // 		  if(hit.TrackID==trackID)
			 // 		    {
			 // 		      std::cout<<"Branch: "<<Hits->GetBranchName()<<" "<<hit.HitPosX<<" "<<hit.HitPosY<<" "<<hit.HitPosZ<<" "<<hit.LayerID<<" "<<hit.Pname<<"\n";
			 // 		    }
			 // 		}
			 // 	    }
			 // 	}
			 

      
		       }

		     //TObjArray* obj_ret = new TObjArray(3);
		     //obj_ret->Add(h_HitPatternX);
		     //obj_ret->Add(h_HitPatternY);
		     //obj_ret->Add(h_ParticlePhi);
		     
		     //return h_HitPatternX;//obj_ret;
		     return std::make_tuple(h_particleStatus,h_HitPatternX,h_HitPatternY,h_ParticlePhi,h_HitPatternX_PSCE,h_HitPatternY_PSCE,h_HitPatternX_PSFE,h_HitPatternY_PSFE,h_CoincidenceMDC_TOF,h_R_Z_MDC,h_X_Z_MDC,h_Y_Z_MDC);
		   };
					    

  // Create the pool of processes
  //ROOT::TProcessExecutor workers(4);

  // Process the TChain
  //auto ObjArray = workers.ProcTree(*Chain, f_Process, "G4Tree");
  //ObjArray->Print();
  //tp.Process(f_Process);
  TH1F* h_particleStatusM = nullptr;
  TH2F* h_HitPatternXM = nullptr;
  TH2F* h_HitPatternYM = nullptr;
  TH2F* h_ParticlePhiM = nullptr;
  TH2F* h_HitPatternX_PSCEM = nullptr;
  TH2F* h_HitPatternY_PSCEM = nullptr;
  TH2F* h_HitPatternX_PSFEM = nullptr;
  TH2F* h_HitPatternY_PSFEM = nullptr;
  TH2F* h_CoincidenceMDC_TOFM = nullptr;
  TH2F* h_R_Z_TOFM = nullptr;
  TH2F* h_X_Z_TOFM = nullptr;
  TH2F* h_Y_Z_TOFM = nullptr;

  std::tie(h_particleStatusM,h_HitPatternXM,h_HitPatternYM,h_ParticlePhiM,h_HitPatternX_PSCEM,h_HitPatternY_PSCEM,h_HitPatternX_PSFEM,h_HitPatternY_PSFEM,h_CoincidenceMDC_TOFM,h_R_Z_TOFM,h_X_Z_TOFM,h_Y_Z_TOFM) = f_Process(reader);


  //TH2F* h_HitPatternXM = static_cast<TH2F*>(ObjArray);
  //TH2F* h_HitPatternYM = static_cast<TH2F*>(ObjArray->At(1));
  //TH2F* h_ParticlePhiM = static_cast<TH2F*>(ObjArray->At(2));
  
  auto plotHist = [](auto h, double min=-1, double max=-1, TString nameX="")
		  {
		    TString nameC(h->GetName());
		    nameC+="Canvas";
		    TCanvas* c = new TCanvas(nameC,nameC,600,600);
		    c->cd();
		    h->Draw();

		    // TString nameY("Counts / ");
		    // double BinW = h->GetXaxis()->GetBinWidth(1);
		    // BinW *= 1.e3;
		    // int BinW2 = BinW*10;
		    // BinW = BinW2/10.;
		    // nameY+=TString::Format("%.1f",BinW);
		    // nameY+=" MeV";
		    // TString FullTitle(";");
		    // FullTitle+=nameX;
		    // FullTitle+=";";
		    // FullTitle+=nameY;
		    // FullTitle+=";";
		    // h->SetTitle(FullTitle);
		    
		    double min1 = min < 0 ? h->GetXaxis()->GetXmin() : min;
		    double max1 = max < 0 ? h->GetXaxis()->GetXmax() : max;
		    h->GetXaxis()->SetRangeUser(min,max);
		    h->Draw();
		    c->Draw("colz");
		  };




  
  //auto h_HitPatternXM = h_HitPatternX.Merge(); 
  //auto h_HitPatternYM = h_HitPatternY.Merge(); 
  //auto h_ParticlePhiM = h_ParticlePhi.Merge(); 

  if(nameOut == "")
    {
      plotHist(h_particleStatusM);
      plotHist(h_HitPatternXM);
      plotHist(h_HitPatternYM);
      plotHist(h_ParticlePhiM);
      plotHist(h_HitPatternX_PSCEM);
      plotHist(h_HitPatternY_PSCEM);
      plotHist(h_HitPatternX_PSFEM);
      plotHist(h_HitPatternY_PSFEM);
      plotHist(h_CoincidenceMDC_TOFM);
      plotHist(h_R_Z_TOFM);
      plotHist(h_X_Z_TOFM);
      plotHist(h_Y_Z_TOFM);
    }
  else
    {
      std::cout << "Saving: "<<nameOut<<"\n";
      TFile* f_out = new TFile(nameOut.c_str(),"RECREATE");
      f_out->cd();
      h_particleStatusM->Write();

      h_HitPatternXM->Write();
      h_HitPatternYM->Write();
      h_ParticlePhiM->Write();

      h_HitPatternX_PSCEM->Write();
      h_HitPatternY_PSCEM->Write();
      h_HitPatternX_PSFEM->Write();
      h_HitPatternY_PSFEM->Write();
      h_CoincidenceMDC_TOFM->Write();

      h_R_Z_TOFM->Write();
      f_out->Close();
    }
  std::cout<<"Done !\n";
}

  

void runDose(const std::string& nameList)
{
  //int nthreads = 4;
  //ROOT::EnableImplicitMT(

  TChain* Chain = new TChain("G4Tree");
  TFile* f_first = nullptr;
  std::cout << "Files :" << nameList << std::endl;
  if(nameList.find(".root") != std::string::npos)
    {
      std::cout << "Load from single file " << nameList << std::endl;
      int temp_nb = Chain->AddFile(nameList.c_str());
      f_first = new TFile(nameList.c_str());
      std::cout << " Loaded " << temp_nb << " files " << std::endl;

    }
  else
    {
      std::cout << "Adding Chain from List files" << std::endl;
      std::ifstream List(nameList.c_str());
      std::string infiles;
      int nb_file = 0;
      while(std::getline(List, infiles))
        {
	  std::cout << infiles << std::endl;
	  int temp_nb = Chain->AddFile(infiles.c_str());
          if(nb_file==0)
	    f_first = new TFile(infiles.c_str());
	  ++nb_file;
	}
      std::cout << " Loaded " << nb_file << " files " << std::endl;
    }

  //TTree* tt = (TTree*) f->Get("G4Tree");
  std::vector<std::string>* nameDetInFile = (std::vector<std::string>*)(f_first->Get("nameDet"));
  std::vector<std::string> nameDet ;
  for(auto name : *nameDetInFile)
    nameDet.emplace_back(name);
  f_first->Close();
  //f_first->Delete();
  std::cout<<" load nameDet done !\n";
  for(auto name : nameDet)
    std::cout<<name<<"\n";


  std::set<std::string> setNameDet;
  for(size_t id = 0; id < nameDet.size(); ++id)
    {
      if(nameDet[id]=="Si1_Strip_log_x")
	setNameDet.insert(nameDet[id]);
      if(nameDet[id]=="Si1_Strip_log_y")
	setNameDet.insert(nameDet[id]);
      if(nameDet[id]=="Si2_Strip_log_x")
	setNameDet.insert(nameDet[id]);
      if(nameDet[id]=="Si2_Strip_log_y")
	setNameDet.insert(nameDet[id]);
    }

  std::unordered_map<std::string,double> detector_mass = {{"Si1",4.*4.*0.03*2.333},{"Si2",6.*6.*0.03*2.333},{"MiniFiber",10*10*0.5*1.032}};
  std::vector<std::string> orderDetName;

  //ROOT::TTreeProcessorMT tp(Chain);
  TTreeReader reader(Chain);

  //ROOT::TThreadedObject<TH2F> h_HitPatternXAll("h_HitPatternX","h_HitPatternX",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_HitPatternYAll("h_HitPatternY","h_HitPatternY",2000,-20,20,20,0,20);
  //ROOT::TThreadedObject<TH2F> h_ParticlePhiAll("h_Phi","h_Phi",360*5,-180,180,30,0,30);

  std::vector<TTreeReaderArray<TG4Sol_Hit>*> AllHits;
  for(auto name : setNameDet)
    {
      AllHits.emplace_back(new TTreeReaderArray<TG4Sol_Hit>(reader,name.c_str()));
      for(auto it_detmass : detector_mass)
	{
	  if(auto it_s = name.find(it_detmass.first); it_s!=std::string::npos)
	    orderDetName.emplace_back(it_detmass.first);
	}
    }
  const auto Entries = Chain->GetEntries();
  std::cout << " Entries :" << Entries << std::endl;
  int timing = 0;
  int first_event = 0;

  std::unordered_map<std::string,std::array<double,2> > TotalEnergy = {{"Si1",{0.,0.}},{"Si2",{0.,0.}}};

  while(reader.Next())
    {
      int nb = reader.GetCurrentEntry();
      if(static_cast<int>(static_cast<double>(nb) / static_cast<double>(Entries) * 10) == timing)
	{
	  std::cout <<" Processing :" << timing * 10 << "%  Event #"<<nb<<" \n";
	  ++timing;
	}

      for(size_t it_br = 0 ;it_br < AllHits.size(); ++it_br)
	{

	  for(auto hit : *AllHits[it_br])
	    {
	      double tempE = hit.Energy;
	      TotalEnergy[orderDetName[it_br]][0] += tempE;
	      TotalEnergy[orderDetName[it_br]][1] += tempE*tempE;
	    }
	}
    }

  for(auto it_E : TotalEnergy)
    {
      std::cout<<" det :"<<it_E.first<<" "<<it_E.second[0]<<" "<<it_E.second[1]<<"\n";

      double mass = detector_mass[it_E.first]*1e-3; // kg

      double tE = it_E.second[0];
      double tE2 = it_E.second[1];

      double rmsE = tE2 /(double)Entries - tE*tE/(double)Entries/(double)Entries;
      rmsE = TMath::Sqrt(rmsE);
      rmsE *= (double)Entries;

      std::cout<<" --> total E [MeV] = "<< tE <<  " +- "<<rmsE<<" ";

      tE *= TMath::Qe()*1e6; // J
      rmsE *= TMath::Qe()*1e6; //J

      std::cout<<" [J] = "<<tE<<" +- "<<rmsE<<"\n";

      double dose = tE / mass ; // Gray = J/kg
      double rms_dose = rmsE/mass;

      std::cout<<" --> dose [Gray] = "<<dose<<" +- "<<rms_dose<<"\n";
      std::cout<<" --> dose [Gray/reaction] = "<<dose/(double)Entries<<" +- "<<rms_dose/(double)Entries<<"\n";
    }

  std::cout<<"Done !\n";
}


void runDataOld(std::string nameF)
{
  TFile* f = new TFile(nameF.c_str());
  TTree* tt = (TTree*) f->Get("G4Tree");
  std::vector<std::string>* nameDet = (std::vector<std::string>*)(f->Get("nameDet"));

  std::cout<<" load nameDet done !\n";
  for(auto name : *nameDet)
    std::cout<<name<<"\n";

  TG4Sol_Event* event = new TG4Sol_Event;
  tt->SetBranchAddress("TG4Sol_Event",&event);

  std::vector<TClonesArray*> AllHits;
  for(auto name : *nameDet)
    {
      AllHits.emplace_back( new TClonesArray("TG4Sol_Hit",20));
      tt->SetBranchAddress(name.c_str(),&AllHits.back());
    }
  
  for(int i = 0;i<tt->GetEntries();++i)
    {
      tt->GetEntry(i);
      std::cout<<"Event #"<<i<<"\n";
      for(auto nameBeam : event->BeamNames)
	std::cout<<"beam : "<<nameBeam<<"\n";
      
      for(size_t k=0; k<AllHits.size();++k)
	{
	  TClonesArray* Hits = AllHits[k];
	  std::cout<<"Branch: "<<nameDet->at(k)<<"\n";

      for(int j=0;j<Hits->GetEntries();++j)
	{
	  TG4Sol_Hit* hit = dynamic_cast<TG4Sol_Hit*>(Hits->At(j));
	  hit->Print();
	}

      
	}
      
    }
}


void runPhaseSpace(std::string nameF)
{
  TFile* f = new TFile(nameF.c_str());
  TTree* tt = (TTree*) f->Get("G4Tree");
 
  tt->Draw("DaughterMomentums_Z:TMath::ATan2(sqrt(DaughterMomentums_X*DaughterMomentums_X+DaughterMomentums_Y*DaughterMomentums_Y),DaughterMomentums_Z)*TMath::RadToDeg()>>h1(220,0,22,150,0,1.5)","DaughterMasses<1","goff");

  TH2F* h1 = (TH2F*)gDirectory->Get("h1");
  TH2F* h11 = (TH2F*) h1->Clone();
  h11->SetName("piTh");
  h11->SetTitle("piTh;#pi- polar angle (degree);#pi- Pz (GeV/c);");
  h11->SetDirectory(0);

  tt->Draw("DaughterMomentums_Z:TMath::ATan2(sqrt(DaughterMomentums_X*DaughterMomentums_X+DaughterMomentums_Y*DaughterMomentums_Y),DaughterMomentums_Z)*TMath::RadToDeg()>>h2(220,0,11,150,0,15)","DaughterMasses>1","goff");

  TH2F* h2 = (TH2F*)gDirectory->Get("h2");
  TH2F* h21 = (TH2F*) h2->Clone();
  h21->SetName("He3Th");
  h21->SetTitle("He3Th;^{3}He polar angle (degree);^{3}He Pz (GeV/c);");
  h21->SetDirectory(0);

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  h11->Draw("colz");

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h21->Draw("colz");
  
}


int main(int argc,char** argv)
{
  if (argc != 3)
    {
      std::cout<<"E> runData In.root Out.root\n";
      return -1;
    }

  runData(argv[1],{},argv[2]);
}
