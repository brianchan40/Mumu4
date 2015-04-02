#include "TChain.h"
#include "TString.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include <iostream>
#include "TMath.h"
#include "TLegend.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

void rotation(TVector3 mom_muon, TVector3 vertex_muon, TVector3 vertex_electron, TVector3 *three_vectors){
	//factor is to change MeV/c to cm/s through MeV/c * 10^6eV/1MeV * [(5.344286*10^(-28) (kg*m/s))/ (eV/c)] * 1/(1.8835*10^(-28) kg) * 100cm/1m
	Double_t factor = 2.837423 * (TMath::Power(10, 8));
	Double_t theta = TMath::Pi()/2 - TMath::ATan2(mom_muon.Y(), mom_muon.X());
	//cout << "theta: " << theta << endl;

	TVector3 velocity_muon;
	velocity_muon.SetXYZ(factor * mom_muon.X() * TMath::Cos(theta) - factor * mom_muon.Y() * TMath::Sin(theta), factor * mom_muon.X() * TMath::Sin(theta) + factor * mom_muon.Y() * TMath::Cos(theta), factor * mom_muon.Z());
	/*cout << "Old Velocity of Muon: " << factor * mom_muon.X() << ", " << factor * mom_muon.Y() << endl;
	cout << "calculation process: " << endl;
	cout << "factor = " << factor << ", mom_muon.X() = " << mom_muon.X() << ", TMath::Cos(theta) = " << TMath::Cos(theta) << endl;
	cout << "factor * mom_muon.X() = " << factor * mom_muon.X() << ", * TMath::Cos(theta) = " << factor * mom_muon.X() * TMath::Cos(theta) << endl;
	cout << "factor = " << factor << ", mom_muon.Y() = " << mom_muon.Y() << ", TMath::Sin(theta) = " << TMath::Sin(theta) << endl;
	cout << "factor * mom_muon.Y() = " << factor * mom_muon.Y() << ", * TMath::Sin(theta) = " << factor * mom_muon.Y() * TMath::Sin(theta) << endl;
	//cout << "Mag of Old vs. New: " << mom_muon.Mag()*factor << " vs. " << velocity_muon.Mag() << endl;*/
	TVector3 new_vertex_muon;
	new_vertex_muon.SetXYZ(vertex_muon.X() * TMath::Cos(theta) - vertex_muon.Y() * TMath::Sin(theta), vertex_muon.X() * TMath::Sin(theta) + vertex_muon.Y() * TMath::Cos(theta), vertex_muon.Z());
	TVector3 new_vertex_electron;
	new_vertex_electron.SetXYZ(vertex_electron.X() * TMath::Cos(theta) - vertex_electron.Y() * TMath::Sin(theta), vertex_electron.X() * TMath::Sin(theta) + vertex_electron.Y() * TMath::Cos(theta), vertex_electron.Z());
	three_vectors[0] = velocity_muon;
	three_vectors[1] = new_vertex_muon;
	three_vectors[2] = new_vertex_electron;
}

Double_t radiuscalc(TVector3 velocity_muon){
	Double_t radius = TMath::Abs((1.8835*(TMath::Power(10,-28)) * TMath::Sqrt(TMath::Power(velocity_muon.X(), 2) + TMath::Power(velocity_muon.Y(), 2)))/(1.602177*(TMath::Power(10, -19)) * 1.5));
	//cout << "Radius: " << radius << endl;
	return radius;
}

void pos_translation(Double_t radius, TVector3 new_vertex_muon, TVector3 new_vertex_electron, TVector3 *two_vectors){
	TVector3 translation;
	translation.SetXYZ(new_vertex_muon.X() + radius, new_vertex_muon.Y(), new_vertex_muon.Z());
	TVector3 final_vertex_muon = new_vertex_muon - translation;
	TVector3 final_vertex_electron = new_vertex_electron - translation;

	two_vectors[0] = final_vertex_muon;
	two_vectors[1] = final_vertex_electron;
}

void neg_translation(Double_t radius, TVector3 new_vertex_muon, TVector3 new_vertex_electron, TVector3 *two_vectors){
	TVector3 translation;
	translation.SetXYZ(new_vertex_muon.X() - radius, new_vertex_muon.Y(), new_vertex_muon.Z());
	TVector3 final_vertex_muon = new_vertex_muon - translation;
	TVector3 final_vertex_electron = new_vertex_electron - translation;

	two_vectors[0] = final_vertex_muon;
	two_vectors[1] = final_vertex_electron;
}

Double_t pos_closest_distance(TVector3 velocity_muon, TVector3 final_vertex_electron, Double_t radius, Double_t ang_vel){
	Double_t t = final_vertex_electron.Z()/velocity_muon.Z();
	cout << "time: " << t << endl;
	Double_t delta_x = (-1) * radius * TMath::Cos(ang_vel * t) - final_vertex_electron.X();
	Double_t delta_y = radius * TMath::Sin(ang_vel * t) - final_vertex_electron.Y();
	Double_t distance = sqrt(TMath::Power(delta_x, 2) + TMath::Power(delta_y,2));
	return distance;
}

Double_t neg_closest_distance(TVector3 velocity_muon, TVector3 final_vertex_electron, Double_t radius, Double_t ang_vel){
	Double_t t = final_vertex_electron.Z()/velocity_muon.Z();
	cout << "time (neg): " << t << endl;
	Double_t delta_x = radius * TMath::Cos(ang_vel * t) - final_vertex_electron.X();
	Double_t delta_y = radius * TMath::Sin(ang_vel * t) - final_vertex_electron.Y();
	Double_t distance = sqrt(TMath::Power(delta_x, 2) + TMath::Power(delta_y,2));
	return distance;
}

Double_t distancecalc_pos(TVector3 mom_muon, TVector3 vertex_muon, TVector3 vertex_electron){
	cout << "positive" << endl;
	cout << "mom_muon: " << mom_muon.X() << ", " << mom_muon.Y() << ", " << mom_muon.Z() << endl;
	cout << "vertex_muon: " << vertex_muon.X() << ", " << vertex_muon.Y() << ", " << vertex_muon.Z() << endl;
	cout << "vertex_electron: " << vertex_electron.X() << ", " << vertex_electron.Y() << ", " << vertex_electron.Z() << endl;
	TVector3 three_vectors[3];
	rotation(mom_muon, vertex_muon, vertex_electron, three_vectors);
	TVector3 velocity_muon = three_vectors[0];
	TVector3 new_vertex_muon = three_vectors[1];
	TVector3 new_vertex_electron = three_vectors[2];
	cout << "velocity_muon: " << velocity_muon.X() << ", " << velocity_muon.Y() << ", " << velocity_muon.Z() << endl;
	cout << "new_vertex_muon: " << new_vertex_muon.X() << ", " << new_vertex_muon.Y() << ", " << new_vertex_muon.Z() << endl;
	cout << "new_vertex_electron: " << new_vertex_electron.X() << ", " << new_vertex_electron.Y() << ", " << new_vertex_electron.Z() << endl;

	Double_t radius = radiuscalc(velocity_muon);
	//cout << "radius: " << radius << endl;
	
	TVector3 two_vectors[2];
	pos_translation(radius, new_vertex_muon, new_vertex_electron, two_vectors);
	TVector3 final_vertex_muon = two_vectors[0];
	cout << "final_vertex_muon: " << final_vertex_muon.X() << ", " << final_vertex_muon.Y() << ", " << final_vertex_muon.Z() << endl;
	TVector3 final_vertex_electron = two_vectors[1];
	cout << "final_vertex_electron: " << final_vertex_electron.X() << ", " << final_vertex_electron.Y() << ", " << final_vertex_electron.Z() << endl;

	Double_t ang_vel = 1.275935595 * (TMath::Power(10, 9));
	
	Double_t distance = pos_closest_distance(velocity_muon, final_vertex_electron, radius, ang_vel);
	cout << "distance: " << distance << endl;

	return distance;
}

Double_t distancecalc_neg(TVector3 mom_muon, TVector3 vertex_muon, TVector3 vertex_electron){
	cout << "negative" << endl;
	cout << "mom_muon: " << mom_muon.X() << ", " << mom_muon.Y() << ", " << mom_muon.Z() << endl;
	cout << "vertex_muon: " << vertex_muon.X() << ", " << vertex_muon.Y() << ", " << vertex_muon.Z() << endl;
	cout << "vertex_electron: " << vertex_electron.X() << ", " << vertex_electron.Y() << ", " << vertex_electron.Z() << endl;
	TVector3 three_vectors[3];
	rotation(mom_muon, vertex_muon, vertex_electron, three_vectors);
	TVector3 velocity_muon = three_vectors[0];
	TVector3 new_vertex_muon = three_vectors[1];
	TVector3 new_vertex_electron = three_vectors[2];
	cout << "velocity_muon: " << velocity_muon.X() << ", " << velocity_muon.Y() << ", " << velocity_muon.Z() << endl;
	cout << "new_vertex_muon: " << new_vertex_muon.X() << ", " << new_vertex_muon.Y() << ", " << new_vertex_muon.Z() << endl;
	cout << "new_vertex_electron: " << new_vertex_electron.X() << ", " << new_vertex_electron.Y() << ", " << new_vertex_electron.Z() << endl;

	Double_t radius = radiuscalc(velocity_muon);
	
	TVector3 two_vectors[2];
	neg_translation(radius, new_vertex_muon, new_vertex_electron, two_vectors);
	TVector3 final_vertex_muon = two_vectors[0];
	cout << "final_vertex_muon: " << final_vertex_muon.X() << ", " << final_vertex_muon.Y() << ", " << final_vertex_muon.Z() << endl;
	TVector3 final_vertex_electron = two_vectors[1];
	cout << "final_vertex_electron: " << final_vertex_electron.X() << ", " << final_vertex_electron.Y() << ", " << final_vertex_electron.Z() << endl;

	Double_t ang_vel = 1.275935595 * (TMath::Power(10, 9));
	
	Double_t distance = neg_closest_distance(velocity_muon, final_vertex_electron, radius, ang_vel);
	cout << "distance: " << distance << endl;

	return distance;
}

void Helix(){
  	using namespace std;

	//create chain ch1, which only takes the variable ntp1011
	TChain *ch1 = new TChain("ntp1011");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch1->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch1->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch1->Add("tauusernm-mumucompton-r2on-2.root");
  	ch1->Add("tauusernm-mumucompton-r2on-4.root");
  	ch1->Add("tauusernm-mumucompton-r2on-5.root");
  	ch1->Add("tauusernm-mumucompton-r2on-9.root");
  	ch1->Add("tauusernm-mumucompton-r2on-12.root");
  	ch1->Add("tauusernm-mumucompton-r2on-14.root");
  	ch1->Add("tauusernm-mumucompton-r2on-17.root");
  	ch1->Add("tauusernm-mumucompton-r2on-19.root");
  	ch1->Add("tauusernm-mumucompton-r2on-21.root");
  	ch1->Add("tauusernm-mumucompton-r2on-24.root");
  	ch1->Add("tauusernm-mumucompton-r2on-27.root");
  	ch1->Add("tauusernm-mumucompton-r3on-1.root");
  	ch1->Add("tauusernm-mumucompton-r3on-7.root");
  	ch1->Add("tauusernm-mumucompton-r3on-8.root");
  	ch1->Add("tauusernm-mumucompton-r4on-1.root");
  	ch1->Add("tauusernm-mumucompton-r4on-7.root");
  	ch1->Add("tauusernm-mumucompton-r4on-8.root");
  	ch1->Add("tauusernm-mumucompton-r4on-11.root");
  	ch1->Add("tauusernm-mumucompton-r4on-13.root");
  	ch1->Add("tauusernm-mumucompton-r4on-15.root");
  	ch1->Add("tauusernm-mumucompton-r5on-2.root");
  	ch1->Add("tauusernm-mumucompton-r5on-8.root");
  	ch1->Add("tauusernm-mumucompton-r5on-16.root");
  	ch1->Add("tauusernm-mumucompton-r6on-2.root");
  	ch1->Add("tauusernm-mumucompton-r6on-3.root");
  	ch1->Add("tauusernm-mumucompton-r6on-5.root");
  	ch1->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch1->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch1->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}

	//create chain ch2, which only takes the variable ntp1012
	TChain *ch2 = new TChain("ntp1012");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch2->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch2->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch2->Add("tauusernm-mumucompton-r2on-2.root");
  	ch2->Add("tauusernm-mumucompton-r2on-4.root");
  	ch2->Add("tauusernm-mumucompton-r2on-5.root");
  	ch2->Add("tauusernm-mumucompton-r2on-9.root");
  	ch2->Add("tauusernm-mumucompton-r2on-12.root");
  	ch2->Add("tauusernm-mumucompton-r2on-14.root");
  	ch2->Add("tauusernm-mumucompton-r2on-17.root");
  	ch2->Add("tauusernm-mumucompton-r2on-19.root");
  	ch2->Add("tauusernm-mumucompton-r2on-21.root");
  	ch2->Add("tauusernm-mumucompton-r2on-24.root");
  	ch2->Add("tauusernm-mumucompton-r2on-27.root");
  	ch2->Add("tauusernm-mumucompton-r3on-1.root");
  	ch2->Add("tauusernm-mumucompton-r3on-7.root");
  	ch2->Add("tauusernm-mumucompton-r3on-8.root");
  	ch2->Add("tauusernm-mumucompton-r4on-1.root");
  	ch2->Add("tauusernm-mumucompton-r4on-7.root");
  	ch2->Add("tauusernm-mumucompton-r4on-8.root");
  	ch2->Add("tauusernm-mumucompton-r4on-11.root");
  	ch2->Add("tauusernm-mumucompton-r4on-13.root");
  	ch2->Add("tauusernm-mumucompton-r4on-15.root");
  	ch2->Add("tauusernm-mumucompton-r5on-2.root");
  	ch2->Add("tauusernm-mumucompton-r5on-8.root");
  	ch2->Add("tauusernm-mumucompton-r5on-16.root");
  	ch2->Add("tauusernm-mumucompton-r6on-2.root");
  	ch2->Add("tauusernm-mumucompton-r6on-3.root");
  	ch2->Add("tauusernm-mumucompton-r6on-5.root");
  	ch2->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch2->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch2->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}

	//create chain ch3, which only takes the variable ntp1013
	TChain *ch3 = new TChain("ntp1013");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch3->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch3->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch3->Add("tauusernm-mumucompton-r2on-2.root");
  	ch3->Add("tauusernm-mumucompton-r2on-4.root");
  	ch3->Add("tauusernm-mumucompton-r2on-5.root");
  	ch3->Add("tauusernm-mumucompton-r2on-9.root");
  	ch3->Add("tauusernm-mumucompton-r2on-12.root");
  	ch3->Add("tauusernm-mumucompton-r2on-14.root");
  	ch3->Add("tauusernm-mumucompton-r2on-17.root");
  	ch3->Add("tauusernm-mumucompton-r2on-19.root");
  	ch3->Add("tauusernm-mumucompton-r2on-21.root");
  	ch3->Add("tauusernm-mumucompton-r2on-24.root");
  	ch3->Add("tauusernm-mumucompton-r2on-27.root");
  	ch3->Add("tauusernm-mumucompton-r3on-1.root");
  	ch3->Add("tauusernm-mumucompton-r3on-7.root");
  	ch3->Add("tauusernm-mumucompton-r3on-8.root");
  	ch3->Add("tauusernm-mumucompton-r4on-1.root");
  	ch3->Add("tauusernm-mumucompton-r4on-7.root");
  	ch3->Add("tauusernm-mumucompton-r4on-8.root");
  	ch3->Add("tauusernm-mumucompton-r4on-11.root");
  	ch3->Add("tauusernm-mumucompton-r4on-13.root");
  	ch3->Add("tauusernm-mumucompton-r4on-15.root");
  	ch3->Add("tauusernm-mumucompton-r5on-2.root");
  	ch3->Add("tauusernm-mumucompton-r5on-8.root");
  	ch3->Add("tauusernm-mumucompton-r5on-16.root");
  	ch3->Add("tauusernm-mumucompton-r6on-2.root");
  	ch3->Add("tauusernm-mumucompton-r6on-3.root");
  	ch3->Add("tauusernm-mumucompton-r6on-5.root");
  	ch3->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch3->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch3->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}

	//create chain ch5, which only takes the variable ntp1015
	TChain *ch5 = new TChain("ntp1015");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch5->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch5->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch5->Add("tauusernm-mumucompton-r2on-2.root");
  	ch5->Add("tauusernm-mumucompton-r2on-4.root");
  	ch5->Add("tauusernm-mumucompton-r2on-5.root");
  	ch5->Add("tauusernm-mumucompton-r2on-9.root");
  	ch5->Add("tauusernm-mumucompton-r2on-12.root");
  	ch5->Add("tauusernm-mumucompton-r2on-14.root");
  	ch5->Add("tauusernm-mumucompton-r2on-17.root");
  	ch5->Add("tauusernm-mumucompton-r2on-19.root");
  	ch5->Add("tauusernm-mumucompton-r2on-21.root");
  	ch5->Add("tauusernm-mumucompton-r2on-24.root");
  	ch5->Add("tauusernm-mumucompton-r2on-27.root");
  	ch5->Add("tauusernm-mumucompton-r3on-1.root");
  	ch5->Add("tauusernm-mumucompton-r3on-7.root");
  	ch5->Add("tauusernm-mumucompton-r3on-8.root");
  	ch5->Add("tauusernm-mumucompton-r4on-1.root");
  	ch5->Add("tauusernm-mumucompton-r4on-7.root");
  	ch5->Add("tauusernm-mumucompton-r4on-8.root");
  	ch5->Add("tauusernm-mumucompton-r4on-11.root");
  	ch5->Add("tauusernm-mumucompton-r4on-13.root");
  	ch5->Add("tauusernm-mumucompton-r4on-15.root");
  	ch5->Add("tauusernm-mumucompton-r5on-2.root");
  	ch5->Add("tauusernm-mumucompton-r5on-8.root");
  	ch5->Add("tauusernm-mumucompton-r5on-16.root");
  	ch5->Add("tauusernm-mumucompton-r6on-2.root");
  	ch5->Add("tauusernm-mumucompton-r6on-3.root");
  	ch5->Add("tauusernm-mumucompton-r6on-5.root");
  	ch5->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch5->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch5->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}

  	//create chain ch7, which only takes the variable ntp1017
	TChain *ch7 = new TChain("ntp1017");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch7->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch7->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch7->Add("tauusernm-mumucompton-r2on-2.root");
  	ch7->Add("tauusernm-mumucompton-r2on-4.root");
  	ch7->Add("tauusernm-mumucompton-r2on-5.root");
  	ch7->Add("tauusernm-mumucompton-r2on-9.root");
  	ch7->Add("tauusernm-mumucompton-r2on-12.root");
  	ch7->Add("tauusernm-mumucompton-r2on-14.root");
  	ch7->Add("tauusernm-mumucompton-r2on-17.root");
  	ch7->Add("tauusernm-mumucompton-r2on-19.root");
  	ch7->Add("tauusernm-mumucompton-r2on-21.root");
  	ch7->Add("tauusernm-mumucompton-r2on-24.root");
  	ch7->Add("tauusernm-mumucompton-r2on-27.root");
  	ch7->Add("tauusernm-mumucompton-r3on-1.root");
  	ch7->Add("tauusernm-mumucompton-r3on-7.root");
  	ch7->Add("tauusernm-mumucompton-r3on-8.root");
  	ch7->Add("tauusernm-mumucompton-r4on-1.root");
  	ch7->Add("tauusernm-mumucompton-r4on-7.root");
  	ch7->Add("tauusernm-mumucompton-r4on-8.root");
  	ch7->Add("tauusernm-mumucompton-r4on-11.root");
  	ch7->Add("tauusernm-mumucompton-r4on-13.root");
  	ch7->Add("tauusernm-mumucompton-r4on-15.root");
  	ch7->Add("tauusernm-mumucompton-r5on-2.root");
  	ch7->Add("tauusernm-mumucompton-r5on-8.root");
  	ch7->Add("tauusernm-mumucompton-r5on-16.root");
  	ch7->Add("tauusernm-mumucompton-r6on-2.root");
  	ch7->Add("tauusernm-mumucompton-r6on-3.root");
  	ch7->Add("tauusernm-mumucompton-r6on-5.root");
  	ch7->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch7->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch7->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}

  	//create chain ch33, which only takes the variable ntp1033
	TChain *ch33 = new TChain("ntp1033");
	for (int i = 1 ; i <= 4 ; ++i) {
    	ch33->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	for (int i = 8 ; i <= 9 ; ++i) {
    	ch33->AddFile(TString::Format("tauusernm-mumucompton-r1on-%d.root", i));
  	}
  	ch33->Add("tauusernm-mumucompton-r2on-2.root");
  	ch33->Add("tauusernm-mumucompton-r2on-4.root");
  	ch33->Add("tauusernm-mumucompton-r2on-5.root");
  	ch33->Add("tauusernm-mumucompton-r2on-9.root");
  	ch33->Add("tauusernm-mumucompton-r2on-12.root");
  	ch33->Add("tauusernm-mumucompton-r2on-14.root");
  	ch33->Add("tauusernm-mumucompton-r2on-17.root");
  	ch33->Add("tauusernm-mumucompton-r2on-19.root");
  	ch33->Add("tauusernm-mumucompton-r2on-21.root");
  	ch33->Add("tauusernm-mumucompton-r2on-24.root");
  	ch33->Add("tauusernm-mumucompton-r2on-27.root");
  	ch33->Add("tauusernm-mumucompton-r3on-1.root");
  	ch33->Add("tauusernm-mumucompton-r3on-7.root");
  	ch33->Add("tauusernm-mumucompton-r3on-8.root");
  	ch33->Add("tauusernm-mumucompton-r4on-1.root");
  	ch33->Add("tauusernm-mumucompton-r4on-7.root");
  	ch33->Add("tauusernm-mumucompton-r4on-8.root");
  	ch33->Add("tauusernm-mumucompton-r4on-11.root");
  	ch33->Add("tauusernm-mumucompton-r4on-13.root");
  	ch33->Add("tauusernm-mumucompton-r4on-15.root");
  	ch33->Add("tauusernm-mumucompton-r5on-2.root");
  	ch33->Add("tauusernm-mumucompton-r5on-8.root");
  	ch33->Add("tauusernm-mumucompton-r5on-16.root");
  	ch33->Add("tauusernm-mumucompton-r6on-2.root");
  	ch33->Add("tauusernm-mumucompton-r6on-3.root");
  	ch33->Add("tauusernm-mumucompton-r6on-5.root");
  	ch33->Add("tauusernm-mumucompton-r6on-14.root");
  	for (int i = 7 ; i <= 9 ; ++i) {
    	ch33->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	for (int i = 16 ; i <= 20 ; ++i) {
    	ch33->AddFile(TString::Format("tauusernm-mumucompton-r6on-%d.root", i));
  	}
  	//Declaring Variables
  	//1011
	int nevents1 = ch1->GetEntries();
	Int_t nTruthCands1;
	Int_t True_pid1[1000];
	Float_t True_px_lab1[1000], True_py_lab1[1000], True_pz_lab1[1000], True_xVtx1[1000], True_yVtx1[1000], True_zVtx1[1000];
	Int_t True_moth_idxInTruthList1[1000], True_isthep1[1000];
	//1012
	int nevents2 = ch2->GetEntries();
	Int_t nTruthCands2;
	Int_t True_pid2[1000];
	Float_t True_px_lab2[1000], True_py_lab2[1000], True_pz_lab2[1000], True_xVtx2[1000], True_yVtx2[1000], True_zVtx2[1000];
	Int_t True_moth_idxInTruthList2[1000], True_isthep2[1000];
	//1013
	int nevents3 = ch3->GetEntries();
	Int_t nTruthCands3;
	Int_t True_pid3[1000];
	Float_t True_px_lab3[1000], True_py_lab3[1000], True_pz_lab3[1000], True_xVtx3[1000], True_yVtx3[1000], True_zVtx3[1000];
	Int_t True_moth_idxInTruthList3[1000], True_isthep3[1000];
	//1015
	int nevents5 = ch5->GetEntries();
	Int_t nTruthCands5;
	Int_t True_pid5[1000];
	Float_t True_px_lab5[1000], True_py_lab5[1000], True_pz_lab5[1000], True_xVtx5[1000], True_yVtx5[1000], True_zVtx5[1000];
	Int_t True_moth_idxInTruthList5[1000], True_isthep5[1000];
	//1017
	int nevents7 = ch7->GetEntries();
	Int_t nTruthCands7;
	Int_t True_pid7[1000];
	Float_t True_px_lab7[1000], True_py_lab7[1000], True_pz_lab7[1000], True_xVtx7[1000], True_yVtx7[1000], True_zVtx7[1000];
	Int_t True_moth_idxInTruthList7[1000], True_isthep7[1000];
	//1033
	int nevents33 = ch33->GetEntries();
	Int_t nTruthCands33;
	Int_t True_pid33[1000];
	Float_t True_px_lab33[1000], True_py_lab33[1000], True_pz_lab33[1000], True_xVtx33[1000], True_yVtx33[1000], True_zVtx33[1000];
	Int_t True_moth_idxInTruthList33[1000], True_isthep33[1000];

	//Setting variables to branch addresses 
	//1011
	ch1->SetBranchAddress ("nTruthCands", &nTruthCands1);
	ch1->SetBranchAddress ("True_pid", True_pid1);
	ch1->SetBranchAddress ("True_px_lab", True_px_lab1);
	ch1->SetBranchAddress ("True_py_lab", True_py_lab1);
	ch1->SetBranchAddress ("True_pz_lab", True_pz_lab1);
	ch1->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList1);
	ch1->SetBranchAddress ("True_xVtx", True_xVtx1);
	ch1->SetBranchAddress ("True_yVtx", True_yVtx1);
	ch1->SetBranchAddress ("True_zVtx", True_zVtx1);
	ch1->SetBranchAddress ("True_isthep", True_isthep1);

	//1012
	ch2->SetBranchAddress ("nTruthCands", &nTruthCands2);
	ch2->SetBranchAddress ("True_pid", True_pid2);
	ch2->SetBranchAddress ("True_px_lab", True_px_lab2);
	ch2->SetBranchAddress ("True_py_lab", True_py_lab2);
	ch2->SetBranchAddress ("True_pz_lab", True_pz_lab2);
	ch2->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList2);
	ch2->SetBranchAddress ("True_xVtx", True_xVtx2);
	ch2->SetBranchAddress ("True_yVtx", True_yVtx2);
	ch2->SetBranchAddress ("True_zVtx", True_zVtx2);
	ch2->SetBranchAddress ("True_isthep", True_isthep2);

	//1013
	ch3->SetBranchAddress ("nTruthCands", &nTruthCands3);
	ch3->SetBranchAddress ("True_pid", True_pid3);
	ch3->SetBranchAddress ("True_px_lab", True_px_lab3);
	ch3->SetBranchAddress ("True_py_lab", True_py_lab3);
	ch3->SetBranchAddress ("True_pz_lab", True_pz_lab3);
	ch3->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList3);
	ch3->SetBranchAddress ("True_xVtx", True_xVtx3);
	ch3->SetBranchAddress ("True_yVtx", True_yVtx3);
	ch3->SetBranchAddress ("True_zVtx", True_zVtx3);
	ch3->SetBranchAddress ("True_isthep", True_isthep3);

	//1015
	ch5->SetBranchAddress ("nTruthCands", &nTruthCands5);
	ch5->SetBranchAddress ("True_pid", True_pid5);
	ch5->SetBranchAddress ("True_px_lab", True_px_lab5);
	ch5->SetBranchAddress ("True_py_lab", True_py_lab5);
	ch5->SetBranchAddress ("True_pz_lab", True_pz_lab5);
	ch5->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList5);
	ch5->SetBranchAddress ("True_xVtx", True_xVtx5);
	ch5->SetBranchAddress ("True_yVtx", True_yVtx5);
	ch5->SetBranchAddress ("True_zVtx", True_zVtx5);
	ch5->SetBranchAddress ("True_isthep", True_isthep5);

	//1017
	ch7->SetBranchAddress ("nTruthCands", &nTruthCands7);
	ch7->SetBranchAddress ("True_pid", True_pid7);
	ch7->SetBranchAddress ("True_px_lab", True_px_lab7);
	ch7->SetBranchAddress ("True_py_lab", True_py_lab7);
	ch7->SetBranchAddress ("True_pz_lab", True_pz_lab7);
	ch7->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList7);
	ch7->SetBranchAddress ("True_xVtx", True_xVtx7);
	ch7->SetBranchAddress ("True_yVtx", True_yVtx7);
	ch7->SetBranchAddress ("True_zVtx", True_zVtx7);
	ch7->SetBranchAddress ("True_isthep", True_isthep7);

	//1033
	ch33->SetBranchAddress ("nTruthCands", &nTruthCands33);
	ch33->SetBranchAddress ("True_pid", True_pid33);
	ch33->SetBranchAddress ("True_px_lab", True_px_lab33);
	ch33->SetBranchAddress ("True_py_lab", True_py_lab33);
	ch33->SetBranchAddress ("True_pz_lab", True_pz_lab33);
	ch33->SetBranchAddress ("True_moth_idxInTruthList", True_moth_idxInTruthList33);
	ch33->SetBranchAddress ("True_xVtx", True_xVtx33);
	ch33->SetBranchAddress ("True_yVtx", True_yVtx33);
	ch33->SetBranchAddress ("True_zVtx", True_zVtx33);
	ch33->SetBranchAddress ("True_isthep", True_isthep33);

	//Canvas
	TCanvas *c1 = new TCanvas("c1", "Origin Distance vs. Minimum Distance", 1500, 800);
	TCanvas *c2 = new TCanvas("c2", "Minimum Distance", 1500, 800);
	TCanvas *c3 = new TCanvas("c3", "Mother-Daughter", 1500, 800);
	TCanvas *c4 = new TCanvas("c4", "Radius", 1500, 800);
	TCanvas *c5 = new TCanvas("c5", "Mother-Daughter vs. Minimum Distance", 1500, 800);
	TCanvas *c6 = new TCanvas("c6", "1st Check: positive vs. negative", 1500, 800);
  TCanvas *c7 = new TCanvas("c7", "<15cm", 1500, 800);

	c2->Divide(2,1);
	c6->Divide(2,1);

	//Histograms
	TH2D *origin_vs_min = new TH2D("origin_vs_min", "Origin Distance vs. Minimum Distance", 2000, 0, 200, 2000, 0, 200);
	TH1F *min_ion_plot = new TH1F("min_ion_plot", "Minimum Distance (Ionization)", 2000, 0, 200);
	TH1F *min_others_plot = new TH1F("min_others_plot", "Minimum Distance (Gamma)", 2000, 0, 200);
	TH1F *mother_daughter = new TH1F("mother_daughter", "Distance between Mother-Daughter", 1000, 0, 500);
	TH2D *origin_vs_mother_daughter = new TH2D("origin_vs_mother_daughter", "Origin Distance vs. Distance between Mother-Daughter", 2000, 0, 200, 2000, 0, 200);
	TH1F *radii = new TH1F("radii", "Radius Values", 6000, 0, 1200);
	TH1F *pos_mother_daughter = new TH1F("pos_mother_daughter", "Distance between Mother-Daughter (Positive Muons)", 1000, 0, 500);
	TH1F *neg_mother_daughter = new TH1F("neg_mother_daughter", "Distance between Mother-Daughter (Negative Muons)", 1000, 0, 500);
  TH2D *limited_origin_vs_mother_daughter = new TH2D("limited_origin_vs_mother_daughter", "Origin Distance vs. Distance between Mother-Daughter (<15cm)", 2000, 0, 200, 2000, 0, 200);

	//1011
	for(int i = 0; i < nevents1; i++){
		ch1->GetEntry(i);
		if (i%1000==0) cerr << i << " / " << nevents1 << endl;

		//pick out the muons
		int id_muon1 = -1;
		int id_muon2 = -1;
		int id_muon3 = -1;
		TVector3 mom_muon1;
    	TVector3 vertex_muon1;
    	TVector3 mom_muon2;
    	TVector3 vertex_muon2;
    	TVector3 mom_muon3;
    	TVector3 vertex_muon3;

		for(int k = 0; k < nTruthCands1; k++){
			if((TMath::Abs(True_pid1[k]) == 13) && (True_xVtx1[k] != -100) && (True_yVtx1[k] != -100) && (True_zVtx1[k] != -100)){
				if(id_muon1 == -1){
					id_muon1 = k;
					mom_muon1.SetXYZ(True_px_lab1[k] * 1000, True_py_lab1[k] * 1000, True_pz_lab1[k] * 1000);
					vertex_muon1.SetXYZ(True_xVtx1[k], True_yVtx1[k], True_zVtx1[k]);
				}
				else if((id_muon1 != -1) && (id_muon2 == -1)){
					id_muon2 = k;
					mom_muon2.SetXYZ(True_px_lab1[k] * 1000, True_py_lab1[k] * 1000, True_pz_lab1[k] * 1000);
          			vertex_muon2.SetXYZ(True_xVtx1[k], True_yVtx1[k], True_zVtx1[k]);
				}
				else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
					id_muon3 = k;
					mom_muon3.SetXYZ(True_px_lab1[k] * 1000, True_py_lab1[k] * 1000, True_pz_lab1[k] * 1000);
          			vertex_muon3.SetXYZ(True_xVtx1[k], True_yVtx1[k], True_zVtx1[k]);
				}
        		else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
          			cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
        		}
			}
		}

		for(int j = 0; j < nTruthCands1; j++){
			
		    //vertex
		    TVector3 origin(True_xVtx1[j], True_yVtx1[j], True_zVtx1[j]);

		    bool distance_valid = true;
		    if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
		    	distance_valid = false;
		    }

		    //compare distances
		    Double_t distance1 = 100000000;
		    Double_t distance2 = 100000000;
		    Double_t distance3 = 100000000;
        	
			if((distance_valid == true) && (True_pid1[j] == 11)){
		   		if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid1[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid1[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid1[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid1[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid1[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid1[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

          		Double_t distances[] = {distance1, distance2, distance3};

          		Double_t distance_min = TMath::MinElement(3, distances);

          		if(distance_min != 100000000){
					if(True_isthep1[j] == 202){
              			//cout << "MINIMUM DISTANCE: " << distance_min << endl;
              			min_ion_plot->Fill(distance_min);
              			if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
            		}
            		else{
              			min_others_plot->Fill(distance_min);
            		}
				}
			}

			if((distance_valid == true) && (True_pid1[j] == 11) && (True_isthep1[j] == 202)){
				int moth_id = True_moth_idxInTruthList1[j];
				if(True_pid1[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab1[moth_id] * 1000, True_py_lab1[moth_id] * 1000, True_pz_lab1[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx1[moth_id], True_yVtx1[moth_id], True_zVtx1[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid1[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab1[moth_id] * 1000, True_py_lab1[moth_id] * 1000, True_pz_lab1[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx1[moth_id], True_yVtx1[moth_id], True_zVtx1[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
        }
    }

  	//1012
  	for(int i = 0; i < nevents2; i++){
    	ch2->GetEntry(i);
    	if (i%1000==0) cerr << i << " / " << nevents2 << endl;

	    //pick out the muons
	    int id_muon1 = -1;
	    int id_muon2 = -1;
	    int id_muon3 = -1;
	    TVector3 mom_muon1;
	    TVector3 vertex_muon1;
	    TVector3 mom_muon2;
	    TVector3 vertex_muon2;
	    TVector3 mom_muon3;
	    TVector3 vertex_muon3;

    	for(int k = 0; k < nTruthCands2; k++){
      		if((TMath::Abs(True_pid2[k]) == 13) && (True_xVtx2[k] != -100) && (True_yVtx2[k] != -100) && (True_zVtx2[k] != -100)){
        		if(id_muon1 == -1){
          			id_muon1 = k;
			        mom_muon1.SetXYZ(True_px_lab2[k] * 1000, True_py_lab2[k] * 1000, True_pz_lab2[k] * 1000);
			        vertex_muon1.SetXYZ(True_xVtx2[k], True_yVtx2[k], True_zVtx2[k]);
			    }
        		else if((id_muon1 != -1) && (id_muon2 == -1)){
          			id_muon2 = k;
          			mom_muon2.SetXYZ(True_px_lab2[k] * 1000, True_py_lab2[k] * 1000, True_pz_lab2[k] * 1000);
          			vertex_muon2.SetXYZ(True_xVtx2[k], True_yVtx2[k], True_zVtx2[k]);
        		}
        		else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
          			id_muon3 = k;
          			mom_muon3.SetXYZ(True_px_lab2[k] * 1000, True_py_lab2[k] * 1000, True_pz_lab2[k] * 1000);
          			vertex_muon3.SetXYZ(True_xVtx2[k], True_yVtx2[k], True_zVtx2[k]);
        		}
        		else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
          			cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
        		}
        	}
    	}

    	for(int j = 0; j < nTruthCands2; j++){
      
        	//vertex
        	TVector3 origin(True_xVtx2[j], True_yVtx2[j], True_zVtx2[j]);

        	//check flags in vertex
        	bool distance_valid = true;
        	if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
          		distance_valid = false;
        	}

        	//compare distances
	        Double_t distance1 = 100000000;
	        Double_t distance2 = 100000000;
	        Double_t distance3 = 100000000;
        
        	if((distance_valid == true) && (True_pid2[j] == 11)){
          		if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid2[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid2[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid2[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid2[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid2[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid2[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

          		Double_t distances[] = {distance1, distance2, distance3};

          		Double_t distance_min = TMath::MinElement(3, distances);

          		if(distance_min != 100000000){
            		if(True_isthep2[j] == 202){
            			min_ion_plot->Fill(distance_min);
              			if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              		}
		            else{
		            	min_others_plot->Fill(distance_min);
		            }
				}
			}

			if((distance_valid == true) && (True_pid2[j] == 11) && (True_isthep2[j] == 202)){
				int moth_id = True_moth_idxInTruthList2[j];
				if(True_pid2[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab2[moth_id] * 1000, True_py_lab2[moth_id] * 1000, True_pz_lab2[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx2[moth_id], True_yVtx2[moth_id], True_zVtx2[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid2[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab2[moth_id] * 1000, True_py_lab2[moth_id] * 1000, True_pz_lab2[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx2[moth_id], True_yVtx2[moth_id], True_zVtx2[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
        }
    }

  	//1013
  	for(int i = 0; i < nevents3; i++){
    	ch3->GetEntry(i);
    	if (i%1000==0) cerr << i << " / " << nevents3 << endl;

	    //pick out the muons
	    int id_muon1 = -1;
	    int id_muon2 = -1;
	    int id_muon3 = -1;
	    TVector3 mom_muon1;
	    TVector3 vertex_muon1;
	    TVector3 mom_muon2;
	    TVector3 vertex_muon2;
	    TVector3 mom_muon3;
	    TVector3 vertex_muon3;

	    for(int k = 0; k < nTruthCands3; k++){
	      if((TMath::Abs(True_pid3[k]) == 13) && (True_xVtx3[k] != -100) && (True_yVtx3[k] != -100) && (True_zVtx3[k] != -100)){
	        if(id_muon1 == -1){
	          id_muon1 = k;
	          mom_muon1.SetXYZ(True_px_lab3[k] * 1000, True_py_lab3[k] * 1000, True_pz_lab3[k] * 1000);
	          vertex_muon1.SetXYZ(True_xVtx3[k], True_yVtx3[k], True_zVtx3[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 == -1)){
	          id_muon2 = k;
	          mom_muon2.SetXYZ(True_px_lab3[k] * 1000, True_py_lab3[k] * 1000, True_pz_lab3[k] * 1000);
	          vertex_muon2.SetXYZ(True_xVtx3[k], True_yVtx3[k], True_zVtx3[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
	          id_muon3 = k;
	          mom_muon3.SetXYZ(True_px_lab3[k] * 1000, True_py_lab3[k] * 1000, True_pz_lab3[k] * 1000);
	          vertex_muon3.SetXYZ(True_xVtx3[k], True_yVtx3[k], True_zVtx3[k]);
	        }
	        else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
	          cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
	        }
	      }
	    }

    	for(int j = 0; j < nTruthCands3; j++){
      
        	//vertex
        	TVector3 origin(True_xVtx3[j], True_yVtx3[j], True_zVtx3[j]);

	        //check flags in vertex
	        bool distance_valid = true;
	        if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
	          distance_valid = false;
	        }

	        //compare distances
	        Double_t distance1 = 100000000;
	        Double_t distance2 = 100000000;
	        Double_t distance3 = 100000000;
	        
	        if((distance_valid == true) && (True_pid3[j] == 11)){
	          	if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid3[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid3[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid3[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid3[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid3[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid3[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

	          Double_t distances[] = {distance1, distance2, distance3};

	          Double_t distance_min = TMath::MinElement(3, distances);

	          if(distance_min != 100000000){
	            if(True_isthep3[j] == 202){
	            	min_ion_plot->Fill(distance_min);
	              		if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
	            }
	            else{
	              min_others_plot->Fill(distance_min);
	            }
			  }
			}

			if((distance_valid == true) && (True_pid3[j] == 11) && (True_isthep3[j] == 202)){
				int moth_id = True_moth_idxInTruthList3[j];
				if(True_pid3[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab3[moth_id] * 1000, True_py_lab3[moth_id] * 1000, True_pz_lab3[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx3[moth_id], True_yVtx3[moth_id], True_zVtx3[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid3[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab3[moth_id] * 1000, True_py_lab3[moth_id] * 1000, True_pz_lab3[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx3[moth_id], True_yVtx3[moth_id], True_zVtx3[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
        }
    }

  	//1015
  	for(int i = 0; i < nevents5; i++){
    	ch5->GetEntry(i);
    	if (i%1000==0) cerr << i << " / " << nevents5 << endl;

    	//pick out the muons
	    int id_muon1 = -1;
	    int id_muon2 = -1;
	    int id_muon3 = -1;
	    TVector3 mom_muon1;
	    TVector3 vertex_muon1;
	    TVector3 mom_muon2;
	    TVector3 vertex_muon2;
	    TVector3 mom_muon3;
	    TVector3 vertex_muon3;

	    for(int k = 0; k < nTruthCands5; k++){
	      if((TMath::Abs(True_pid5[k]) == 13) && (True_xVtx5[k] != -100) && (True_yVtx5[k] != -100) && (True_zVtx5[k] != -100)){
	        if(id_muon1 == -1){
	          id_muon1 = k;
	          mom_muon1.SetXYZ(True_px_lab5[k] * 1000, True_py_lab5[k] * 1000, True_pz_lab5[k] * 1000);
	          vertex_muon1.SetXYZ(True_xVtx5[k], True_yVtx5[k], True_zVtx5[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 == -1)){
	          id_muon2 = k;
	          mom_muon2.SetXYZ(True_px_lab5[k] * 1000, True_py_lab5[k] * 1000, True_pz_lab5[k] * 1000);
	          vertex_muon2.SetXYZ(True_xVtx5[k], True_yVtx5[k], True_zVtx5[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
	          id_muon3 = k;
	          mom_muon3.SetXYZ(True_px_lab5[k] * 1000, True_py_lab5[k] * 1000, True_pz_lab5[k] * 1000);
	          vertex_muon3.SetXYZ(True_xVtx5[k], True_yVtx5[k], True_zVtx5[k]);
	        }
	        else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
	          cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
	        }
	      }
	    }

    	for(int j = 0; j < nTruthCands5; j++){
      
        	//vertex
        	TVector3 origin(True_xVtx5[j], True_yVtx5[j], True_zVtx5[j]);

        	//check flags in vertex
        	bool distance_valid = true;
        	if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
          		distance_valid = false;
        	}

	        //compare distances
	        Double_t distance1 = 100000000;
	        Double_t distance2 = 100000000;
	        Double_t distance3 = 100000000;
	        
        	if((distance_valid == true) && (True_pid5[j] == 11)){
          		if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid5[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid5[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid5[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid5[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid5[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid5[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

	          Double_t distances[] = {distance1, distance2, distance3};

	          Double_t distance_min = TMath::MinElement(3, distances);

	          if(distance_min != 100000000){
	            if(True_isthep5[j] == 202){
	            	min_ion_plot->Fill(distance_min);
	            		if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
	            }
	            else{
	              min_others_plot->Fill(distance_min);
	            }

	          }
			}

			if((distance_valid == true) && (True_pid5[j] == 11) && (True_isthep5[j] == 202)){
				int moth_id = True_moth_idxInTruthList5[j];
				if(True_pid5[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab5[moth_id] * 1000, True_py_lab5[moth_id] * 1000, True_pz_lab5[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx5[moth_id], True_yVtx5[moth_id], True_zVtx5[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid5[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab5[moth_id] * 1000, True_py_lab5[moth_id] * 1000, True_pz_lab5[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx5[moth_id], True_yVtx5[moth_id], True_zVtx5[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
        }
    }

  	//1017
  	for(int i = 0; i < nevents7; i++){
    	ch7->GetEntry(i);
    	if (i%1000==0) cerr << i << " / " << nevents7 << endl;

	    //pick out the muons
	    int id_muon1 = -1;
	    int id_muon2 = -1;
	    int id_muon3 = -1;
	    TVector3 mom_muon1;
	    TVector3 vertex_muon1;
	    TVector3 mom_muon2;
	    TVector3 vertex_muon2;
	    TVector3 mom_muon3;
	    TVector3 vertex_muon3;

	    for(int k = 0; k < nTruthCands7; k++){
	      if((TMath::Abs(True_pid7[k]) == 13) && (True_xVtx7[k] != -100) && (True_yVtx7[k] != -100) && (True_zVtx7[k] != -100)){
	        if(id_muon1 == -1){
	          id_muon1 = k;
	          mom_muon1.SetXYZ(True_px_lab7[k] * 1000, True_py_lab7[k] * 1000, True_pz_lab7[k] * 1000);
	          vertex_muon1.SetXYZ(True_xVtx7[k], True_yVtx7[k], True_zVtx7[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 == -1)){
	          id_muon2 = k;
	          mom_muon2.SetXYZ(True_px_lab7[k] * 1000, True_py_lab7[k] * 1000, True_pz_lab7[k] * 1000);
	          vertex_muon2.SetXYZ(True_xVtx7[k], True_yVtx7[k], True_zVtx7[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
	          id_muon3 = k;
	          mom_muon3.SetXYZ(True_px_lab7[k] * 1000, True_py_lab7[k] * 1000, True_pz_lab7[k] * 1000);
	          vertex_muon3.SetXYZ(True_xVtx7[k], True_yVtx7[k], True_zVtx7[k]);
	        }
	        else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
	          cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
	        }
	      }
	    }

    	for(int j = 0; j < nTruthCands7; j++){
      
	        //vertex
	        TVector3 origin(True_xVtx7[j], True_yVtx7[j], True_zVtx7[j]);

	        //check flags in vertex
	        bool distance_valid = true;
	        if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
	          distance_valid = false;
	        }

	        //compare distances
	        Double_t distance1 = 100000000;
	        Double_t distance2 = 100000000;
	        Double_t distance3 = 100000000;
	        Double_t factor1;
	        Double_t factor2;
	        Double_t factor3;
	        
	        if((distance_valid == true) && (True_pid7[j] == 11)){
          		if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid7[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid7[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid7[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid7[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid7[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid7[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

		        Double_t distances[] = {distance1, distance2, distance3};

		        Double_t distance_min = TMath::MinElement(3, distances);

		        if(distance_min != 100000000){
		          if(True_isthep7[j] == 202){
		          	min_ion_plot->Fill(distance_min);
		            	if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
		          }
		          else{
		              min_others_plot->Fill(distance_min);
		          }
		        }
			}

			if((distance_valid == true) && (True_pid7[j] == 11) && (True_isthep7[j] == 202)){
				int moth_id = True_moth_idxInTruthList7[j];
				if(True_pid7[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab7[moth_id] * 1000, True_py_lab7[moth_id] * 1000, True_pz_lab7[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx7[moth_id], True_yVtx7[moth_id], True_zVtx7[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid7[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab7[moth_id] * 1000, True_py_lab7[moth_id] * 1000, True_pz_lab7[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx7[moth_id], True_yVtx7[moth_id], True_zVtx7[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
    	}
  	}

  	//1033
  	for(int i = 0; i < nevents33; i++){
    	ch33->GetEntry(i);
    	if (i%1000==0) cerr << i << " / " << nevents33 << endl;

	    //pick out the muons
	    int id_muon1 = -1;
	    int id_muon2 = -1;
	    int id_muon3 = -1;
	    TVector3 mom_muon1;
	    TVector3 vertex_muon1;
	    TVector3 mom_muon2;
	    TVector3 vertex_muon2;
	    TVector3 mom_muon3;
	    TVector3 vertex_muon3;

	    for(int k = 0; k < nTruthCands33; k++){
	      if((TMath::Abs(True_pid33[k]) == 13) && (True_xVtx33[k] != -100) && (True_yVtx33[k] != -100) && (True_zVtx33[k] != -100)){
	        if(id_muon1 == -1){
	          id_muon1 = k;
	          mom_muon1.SetXYZ(True_px_lab33[k] * 1000, True_py_lab33[k] * 1000, True_pz_lab33[k] * 1000);
	          vertex_muon1.SetXYZ(True_xVtx33[k], True_yVtx33[k], True_zVtx33[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 == -1)){
	          id_muon2 = k;
	          mom_muon2.SetXYZ(True_px_lab33[k] * 1000, True_py_lab33[k] * 1000, True_pz_lab33[k] * 1000);
	          vertex_muon2.SetXYZ(True_xVtx33[k], True_yVtx33[k], True_zVtx33[k]);
	        }
	        else if((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 == -1)){
	          id_muon3 = k;
	          mom_muon3.SetXYZ(True_px_lab33[k] * 1000, True_py_lab33[k] * 1000, True_pz_lab33[k] * 1000);
	          vertex_muon3.SetXYZ(True_xVtx33[k], True_yVtx33[k], True_zVtx33[k]);
	        }
	        else if ((id_muon1 != -1) && (id_muon2 != -1) && (id_muon3 != -1)){
	          cout << "ERROR!!! MORE THAN 3 MUONS IN EVENT " << i << endl;
	        }
	      }
	    }

    	for(int j = 0; j < nTruthCands33; j++){
      
	        //vertex
	        TVector3 origin(True_xVtx33[j], True_yVtx33[j], True_zVtx33[j]);

	        //check flags in vertex
	        bool distance_valid = true;
	        if((origin.X() == -100) || (origin.Y() == -100) || (origin.Z() == -100)){
	          distance_valid = false;
	        }

	        //compare distances
	        Double_t distance1 = 100000000;
	        Double_t distance2 = 100000000;
	        Double_t distance3 = 100000000;
	        
        	if((distance_valid == true) && (True_pid33[j] == 11)){
          		if(id_muon1 != -1){
		   			//cout << "True_pid1[id_muon1] = " << True_pid1[id_muon1] << endl;
            		if(True_pid33[id_muon1] == -13){
            			distance1 = distancecalc_pos(mom_muon1, vertex_muon1, origin);
            		}
            		else if(True_pid33[id_muon1] == 13){
            			distance1 = distancecalc_neg(mom_muon1, vertex_muon1, origin);
            		}
            	}

          		if(id_muon2 != -1){
		   			//cout << "True_pid1[id_muon2] = " << True_pid1[id_muon2] << endl;
            		if(True_pid33[id_muon2] == -13){
            			distance2 = distancecalc_pos(mom_muon2, vertex_muon2, origin);
            		}
            		else if(True_pid33[id_muon2] == 13){
            			distance2 = distancecalc_neg(mom_muon2, vertex_muon2, origin);
            		}
            	}

          		if(id_muon3 != -1){
		   			//cout << "True_pid1[id_muon3] = " << True_pid1[id_muon3] << endl;
            		if(True_pid33[id_muon3] == -13){
            			distance3 = distancecalc_pos(mom_muon3, vertex_muon3, origin);
            		}
            		else if(True_pid33[id_muon3] == 13){
            			distance3 = distancecalc_neg(mom_muon3, vertex_muon3, origin);
            		}
            	}

          		Double_t distances[] = {distance1, distance2, distance3};

          		Double_t distance_min = TMath::MinElement(3, distances);

          		if(distance_min != 100000000){
            			if(True_isthep33[j] == 202){
            				min_ion_plot->Fill(distance_min);
              				if(distance_min == distance1){
              				TVector3 difference = origin - vertex_muon1;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance2){
              				TVector3 difference = origin - vertex_muon2;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
              			else if(distance_min == distance3){
              				TVector3 difference = origin - vertex_muon3;
              				origin_vs_min->Fill(distance_min, difference.Mag());
              			}
            		}
            		else{
              			min_others_plot->Fill(distance_min);
            		}
				}

        	}

        	if((distance_valid == true) && (True_pid33[j] == 11) && (True_isthep33[j] == 202)){
				int moth_id = True_moth_idxInTruthList33[j];
				if(True_pid33[moth_id] == -13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab33[moth_id] * 1000, True_py_lab33[moth_id] * 1000, True_pz_lab33[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx33[moth_id], True_yVtx33[moth_id], True_zVtx33[moth_id]);
					Double_t min_distance = distancecalc_pos(momentum_muon, vertex_mother_muon, origin);
					
					mother_daughter->Fill(min_distance);
					pos_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
				else if(True_pid33[moth_id] == 13){
					TVector3 momentum_muon;
					TVector3 vertex_mother_muon;
					momentum_muon.SetXYZ(True_px_lab33[moth_id] * 1000, True_py_lab33[moth_id] * 1000, True_pz_lab33[moth_id] * 1000);
					vertex_mother_muon.SetXYZ(True_xVtx33[moth_id], True_yVtx33[moth_id], True_zVtx33[moth_id]);
					Double_t min_distance = distancecalc_neg(momentum_muon, vertex_mother_muon, origin);

					mother_daughter->Fill(min_distance);
					neg_mother_daughter->Fill(min_distance);

					TVector3 difference = origin - vertex_mother_muon;
					origin_vs_mother_daughter->Fill(min_distance, difference.Mag());

          if(TMath::Sqrt(TMath::Power(origin.Y(), 2) + TMath::Power(origin.X(), 2)) <= 15){
            limited_origin_vs_mother_daughter->Fill(min_distance, difference.Mag());
          }

					Double_t factor = 2.837423 * (TMath::Power(10, 8));
					Double_t theta = TMath::Pi()/2 - TMath::ATan2(momentum_muon.Y(), momentum_muon.X());
					TVector3 velocity_muon;
					velocity_muon.SetXYZ(factor * momentum_muon.X() * TMath::Cos(theta) - factor * momentum_muon.Y() * TMath::Sin(theta), factor * momentum_muon.X() * TMath::Sin(theta) + factor * momentum_muon.Y() * TMath::Cos(theta), factor * momentum_muon.Z());
					Double_t radius = radiuscalc(velocity_muon);
					radii->Fill(radius);
				}
			}
    	}
    }

    //Draw
    c1->cd();
    origin_vs_min->GetXaxis()->SetTitle("minimum distance (cm)");
    origin_vs_min->GetYaxis()->SetTitle("origin distance (cm)");
    origin_vs_min->Draw();

    c2->cd(1);
    min_ion_plot->GetXaxis()->SetTitle("minimum distance (cm)");
    min_ion_plot->GetYaxis()->SetTitle("Number");
    min_ion_plot->Draw();

    c2->cd(2);
    min_others_plot->GetXaxis()->SetTitle("minimum distance (cm)");
    min_others_plot->GetYaxis()->SetTitle("Number");
    min_others_plot->Draw();

    c3->cd();
    mother_daughter->GetXaxis()->SetTitle("minimum distance (cm)");
    mother_daughter->GetYaxis()->SetTitle("Number");
    mother_daughter->Draw();
    Double_t small = mother_daughter->Integral(0, 20);
    Double_t large = mother_daughter->Integral(21, 1000);
    cout << "Number from 0 to 10 cm: " << small << endl;
    cout << "Number from 11 to 500 cm: " << large << endl;

    c4->cd();
    radii->GetXaxis()->SetTitle("radius (cm)");
    radii->GetYaxis()->SetTitle("Number");
    radii->Draw();

    c5->cd();
   	origin_vs_mother_daughter->GetXaxis()->SetTitle("minimum distance (cm)");
   	origin_vs_mother_daughter->GetYaxis()->SetTitle("origin distance (cm)");
   	origin_vs_mother_daughter->Draw();

   	c6->cd(1);
   	pos_mother_daughter->GetXaxis()->SetTitle("minimum distance (cm)");
   	pos_mother_daughter->GetYaxis()->SetTitle("Number");
   	pos_mother_daughter->Draw();

   	c6->cd(2);
   	neg_mother_daughter->GetXaxis()->SetTitle("minimum distance (cm)");
   	neg_mother_daughter->GetYaxis()->SetTitle("Number");
   	neg_mother_daughter->Draw();

    c7->cd();
    limited_origin_vs_mother_daughter->GetXaxis()->SetTitle("minimum distance (cm)");
    limited_origin_vs_mother_daughter->GetYaxis()->SetTitle("origin distance (cm)");
    limited_origin_vs_mother_daughter->Draw();
}