#include<fstream>
#include <string>
#include "math.h"
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
using namespace std;




//=====================================================================================================================================================//
//s509 Position Variables 
//

#define TRIG(i) (1 << (i - 1))

/*
   On-spill + FOOT. 
   TRIG_LMU_OUT( 1) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu;
   TRIG_LMU_OUT( 2) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_tofd;
   TRIG_LMU_OUT( 3) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_and;
   TRIG_LMU_OUT( 4) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_and and not in_califa_veto;
   TRIG_LMU_OUT( 5) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_or;
   TRIG_LMU_OUT( 6) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_or and not in_califa_veto;
   TRIG_LMU_OUT( 6) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_neuland;

   On-spill - FOOT. 
   TRIG_LMU_OUT( 7) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu;
   TRIG_LMU_OUT( 8) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_tofd;
   TRIG_LMU_OUT( 9) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_and;
   TRIG_LMU_OUT(10) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_and and not in_califa_veto;
   TRIG_LMU_OUT(11) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_or;
   TRIG_LMU_OUT(12) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_or and not in_califa_veto;
   TRIG_LMU_OUT(12) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_neuland;

   Off-spill. 
   TRIG_LMU_OUT(13) = not BEAM_GATE_AUX and in_califa_or;
   TRIG_LMU_OUT(14) = not BEAM_GATE_AUX and in_neuland;
   TRIG_LMU_OUT(15) = not BEAM_GATE_AUX and in_tofd;

   TRIG_LMU_OUT(16) = not BEAM_GATE_AUX and in_rpc; 
   */
double  mm = pow(10,3);
double  mwz_0 = 0.83*mm;
double  mwz_1 = 1.250*mm;
double  mwz_t = 2.773*mm;
double  fib_ang = 14;//deg
double  fib32_z = 4.862*mm;// from bending
double  fib30_z = 5.148*mm;// from bending//yy
double  fib33_z = 6.462*mm;// from bending
double  fib31_z = 6.5939*mm;// from bending
double  tofd_z = 11.394*mm;// from bending 
double  tofd_bar = 27;//mm

//=====================================================================================================================================================//
//FRS varibales
double FRSZ=0;
double FRSAQ=0;
double FRSB=0;
double FRSG=0;
double FRSBR=0;
double FRST=0;
double FRSX=0;
double FRSM=0;
double FRSZ_G = 7;

//=====================================================================================================================================================//
//Hist Init 
//Int_t Energy_sup = 3500;//s509
Int_t Energy_sup = 2100;//S522

//=====================================================================================================================================================//
//R3BRoot Variables
int setPoint;
Int_t TPAT;

//=====================================================================================================================================================//
//To make More Pretty Plots 
void Hist(){
}



double recon(Double_t dep_e, Double_t theta)
{
	Double_t theta_deg = theta * TMath::RadToDeg();
	if (theta_deg<55.5) { // iphos
		if (dep_e<135 || dep_e > 265)
			return -1000;
		return 12369.8 - 194.693 * dep_e +1.20531 * pow(dep_e,2) - 0.00335935 * pow(dep_e,3) +0.00000353193 * pow(dep_e,4);
	} else if (theta_deg < 70.4) { // barrel 1
		if (dep_e<105 || dep_e > 215)
			return -1000;
		return 5381.41 - 69.8561*dep_e + 0.288122*pow(dep_e,2) - 0.000135574*pow(dep_e,3) - 0.000000950755*pow(dep_e,4);
	} else { // barrel 2
		if (dep_e<100 || dep_e > 205)
			return -1000;
		return 3321.43 - 19.5399*dep_e - 0.19495*pow(dep_e,2) + 0.00196785*pow(dep_e,3) - 0.00000440662*pow(dep_e,4);
	}
}

void MinDist(double ax_right, double bx_right, double ay_right, double by_right, double ax_left, double bx_left, double ay_left, double by_left, double &xv, double &yv, double &zv, double &distmin)
{

	double a1, a2, b1, b2, ap1, ap2, bp1, bp2;

	a1 = bx_left;
	a2 = by_left;
	b1 = ax_left;
	b2 = ay_left;

	ap1 = bx_right;
	ap2 = by_right;
	bp1 = ax_right;
	bp2 = ay_right;

	double alpha, beta, A, B, C;

	alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
	beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);

	A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
	B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
	C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));


	double sol1, solf1;
	double x1,y1,z1,xp,yp,zp;


	sol1 = -(A*alpha + C)/(A*beta + B);
	solf1 = alpha + beta* sol1;


	x1 = a1 + b1*sol1;
	y1 = a2 + b2*sol1;
	z1 = sol1;
	xp = ap1 + bp1*solf1;
	yp = ap2 + bp2*solf1;
	zp = solf1;

	xv = (x1+xp)/2.;
	yv = (y1+yp)/2.;
	zv = (z1+zp)/2.;


	distmin = TMath::Sqrt((x1-xp)*(x1-xp) + (y1-yp)*(y1-yp) + (z1-zp)*(z1-zp));


}
//=====================================================================================================================================================//
//Main()
void vertex(){

	TStopwatch timer;
	timer.Start();
	Hist();

	TFile*in =  new TFile("../unpacked/131_mdf_0001.root","READ");    	

	if(in){std::cout<<"File Found" << std::endl;} 
	if(!in){std::cout<<"File Not Found" << std::endl;} 
	TTree*tree = (TTree*)in->Get("evt");
	if(tree){std::cout<< "Tree Found"<<std::endl;}
	if(!tree){std::cout<< "Tree Not Found"<<std::endl;}
	int nevents = tree->GetEntries();
	std::cout << "Number Of Raw Events: "<< nevents << std::endl; 

	R3BEventHeader*DataCA = new R3BEventHeader();
	TBranch*branchData = tree->GetBranch("EventHeader.");
	if(branchData){std::cout<< "EventHeader  Loaded"<<std::endl;}
	if(!branchData){std::cout<< "EventHeader Not Loaded"<<std::endl;}
	branchData->SetAddress(&DataCA);

	TClonesArray*hitDataFRS = new TClonesArray("R3BFrsData",10);
	TBranch*branchFRSHitData = tree->GetBranch("FrsData");
	if(branchFRSHitData){std::cout<<"FrsData Loaded"<<std::endl;}
	if(!branchFRSHitData){std::cout<<"FrsData Not Loaded"<<std::endl;}
	branchFRSHitData->SetAddress(&hitDataFRS);

	TClonesArray*hitDataMWPC0 = new TClonesArray("R3BMwpcHitData",2);
	TBranch*branchMWPC0HitData  = tree->GetBranch("Mwpc0HitData");
	if(branchMWPC0HitData){std::cout<< "Mwpc0HitData Loaded"<<std::endl;}
	if(!branchMWPC0HitData){std::cout<< "Mwpc0HitData Not Loaded"<<std::endl;}
	branchMWPC0HitData->SetAddress(&hitDataMWPC0);

	TClonesArray*hitDataMWPC1 = new TClonesArray("R3BMwpcHitData",2);
	TBranch*branchMWPC1HitData  = tree->GetBranch("Mwpc1HitData");
	if(branchMWPC1HitData){std::cout<< "MwpC1HitData Loaded"<<std::endl;}
	if(!branchMWPC1HitData){std::cout<< "MwpC1HitData Not  Loaded"<<std::endl;}
	branchMWPC1HitData->SetAddress(&hitDataMWPC1);

	TClonesArray*LosTCalData = new TClonesArray("R3BLosTCalData");
	TBranch*branchLosTCalData  = tree->GetBranch("LosTCal");
	if(branchLosTCalData){std::cout<< "LosTCalData Loaded"<<std::endl;}
	if(!branchLosTCalData){std::cout<< "LosTCalData Not  Loaded"<<std::endl;}
	branchLosTCalData->SetAddress(&LosTCalData);

	TClonesArray *CalifaClusterData = new TClonesArray("R3BCalifaClusterData");
	tree->GetBranch("CalifaClusterData")->SetAutoDelete(kFALSE);
	if(CalifaClusterData){std::cout<<"CalifaHitData  Loaded"<<std::endl;}
	if(!CalifaClusterData){std::cout<<"CalifaHitData Not Loaded"<<std::endl;}
	tree->SetBranchAddress("CalifaClusterData", &CalifaClusterData);

	TClonesArray* TofdHitData = new TClonesArray("R3BTofdHitData",10);
	TBranch* branchTofdHit = tree->GetBranch("TofdHit");
	if(branchTofdHit){std::cout<<"TofdHitData  Loaded"<<std::endl;}
	if(!branchTofdHit){std::cout<<"TofdHitData Not Loaded"<<std::endl;}
	branchTofdHit->SetAddress(&TofdHitData);

	TClonesArray *FootHitData = new TClonesArray("R3BFootHitData");
	tree->GetBranch("FootHitData")->SetAutoDelete(kFALSE);
	if(FootHitData){std::cout<< "FootHitData Loaded"<<std::endl;}
	if(!FootHitData){std::cout<< "FootHitData Not Loaded"<<std::endl;}
	tree->SetBranchAddress("FootHitData", &FootHitData);



	//cut normali
	TCutG *cutg3 = new TCutG("cut103f",7);
	cutg3->SetVarX("FOOT10_X");
	cutg3->SetVarY("FOOT3_X");
	cutg3->SetTitle("Graph");
	cutg3->SetFillStyle(1000);
	cutg3->SetPoint(0,70.8277,57.3048);
	cutg3->SetPoint(1,44.7736,35.7081);
	cutg3->SetPoint(2,44.84,31.7815);
	cutg3->SetPoint(3,47.8309,31.536);
	cutg3->SetPoint(4,71.3594,51.5375);
	cutg3->SetPoint(5,71.0271,54.1144);
	cutg3->SetPoint(6,70.8277,57.3048);

	TCutG *cutg5 = new TCutG("CUTG5",6);
	cutg5->SetVarX("FOOT5_Y");
	cutg5->SetVarY("FOOT12_Y");
	cutg5->SetTitle("Graph");
	cutg5->SetFillStyle(1000);
	cutg5->SetPoint(0,55.0653,41.3125);
	cutg5->SetPoint(1,-50.0653,-27.2083);
	cutg5->SetPoint(2,-49.7971,-39.1875);
	cutg5->SetPoint(3,55.0653,29.3333);
	cutg5->SetPoint(4,55.0653,37);
	cutg5->SetPoint(5,55.0653,41.3125);

	TCutG *cutg6 = new TCutG("CUTG6",6);
	cutg6->SetVarX("FOOT9_Y");
	cutg6->SetVarY("FOOT6_Y");
	cutg6->SetTitle("Graph");
	cutg6->SetFillStyle(1000);
	cutg6->SetPoint(0,52.3834,34.6042);
	cutg6->SetPoint(1,-49.5289,-19.0625);
	cutg6->SetPoint(2,-48.7243,-33.9167);
	cutg6->SetPoint(3,53.188,22.1458);
	cutg6->SetPoint(4,53.188,25.9792);
	cutg6->SetPoint(5,52.3834,34.6042);


	TCutG *cutg8 = new TCutG("CUTG8",8);
	cutg8->SetVarX("FOOT8_X");
	cutg8->SetVarY("FOOT11_X");
	cutg8->SetTitle("Graph");
	cutg8->SetFillStyle(1000);
	cutg8->SetPoint(0,-49.6979,-32.6667);
	cutg8->SetPoint(1,-71.6711,-52.4);
	cutg8->SetPoint(2,-71.4141,-58.4);
	cutg8->SetPoint(3,-68.1374,-58.6667);
	cutg8->SetPoint(4,-45.2647,-37.6);
	cutg8->SetPoint(5,-45.3289,-31.8667);
	cutg8->SetPoint(6,-47.9632,-32.4);
	cutg8->SetPoint(7,-49.6979,-32.6667);

	double sx =0.;//4.26095;
	double sy = 0.;

	double angley=0.;//1.08506; //0.08506;
	double x5=63.395+sx;
	double x12=38.323+sx;
	double x9=-63.935+sx;
	double x6=-39.354+sx;
	double z5=527.848;
	double z12=490.752;
	double z9=527.136;
	double z6=490.413;

	double px5=20*TMath::Cos(-74.64723447*TMath::DegToRad());
	double px12=20*TMath::Cos(-74.67117924*TMath::DegToRad());
	double px6=20*TMath::Cos(75.09256426426*TMath::DegToRad());
	double px9=20*TMath::Cos(75.09398359*TMath::DegToRad());
	double pz5=20*TMath::Sin(-74.64723447*TMath::DegToRad());
	double pz12=20*TMath::Sin(-74.67117924*TMath::DegToRad());
	double pz6=20*TMath::Sin(75.09256426426*TMath::DegToRad());
	double pz9=20*TMath::Sin(75.09398359*TMath::DegToRad());

	//cout<<"px12"<<px12<<endl;
	//cout<<"pz12"<<pz12<<endl;

	double s1x5=63.395+sx-px5;
	double s1x12=38.323+sx-px12;
	double s1x9=-63.935+sx+px9;
	double s1x6=-39.354+sx+px6;
	double s1z5=527.848-pz5;
	double s1z12=490.752-pz12;
	double s1z9=527.136+pz9;
	double s1z6=490.413+pz6;

	double c5;
	double c12;
	double c6;
	double c9;

	TVector3 pos5(x5,0,z5);
	TVector3 pos12(x12,0,z12);
	TVector3 pos6(x6,0,z6);
	TVector3 pos9(x9,0,z9);

	TVector3 pos5_1(s1x5,0,s1z5);
	TVector3 pos12_1(s1x12,0,s1z12);
	TVector3 pos6_1(s1x6,0,s1z6);
	TVector3 pos9_1(s1x9,0,s1z9);

	pos5.RotateY(TMath::DegToRad()*angley);
	pos12.RotateY(TMath::DegToRad()*angley);
	pos9.RotateY(TMath::DegToRad()*angley);
	pos6.RotateY(TMath::DegToRad()*angley);

	pos5_1.RotateY(TMath::DegToRad()*angley);
	pos12_1.RotateY(TMath::DegToRad()*angley);
	pos9_1.RotateY(TMath::DegToRad()*angley);
	pos6_1.RotateY(TMath::DegToRad()*angley);

	c5= ((pos5.X()-pos5_1.X())/(pos5.Z()-pos5_1.Z()));
	c12= ((pos12.X()-pos12_1.X())/(pos12.Z()-pos12_1.Z()));
	c6= ((pos6.X()-pos6_1.X())/(pos6.Z()-pos6_1.Z()));
	c9= ((pos9.X()-pos9_1.X())/(pos9.Z()-pos9_1.Z()));

	double td=439.338211;
	double angley_r=0.89515684;
	double fADet_old[20]={0.,0.,0.,0.,0.,-0.2563176063,0.2506446273,0.,0.,0.2506316425, 0.,0.,-0.2563132822,0.,0.,0.,0.,0.,0.,0 };
	double fADet[20]={0.,0.,0.,0.,0.,c5,c6,0.,0.,c9, 0.,0.,c12,0.,0.,0.,0.,0.,0.,0 };
	double fBDet[20]={0.,0.,0.,0.,0.,pos5.X()-c5*(pos5.Z()),pos6.X()-c6*(pos6.Z()),0.,0.,pos9.X()-c9*(pos9.Z()),0.,0.,pos12.X()-c12*(pos12.Z()),0.,0.,0.,0.,0.,0.,0.};

	//cout<<"x5"<< pos5.X()<<endl;
	//cout<<"z5"<< pos5.Z()<<endl;

	//cout<<"x6"<< pos6.X()<<endl;
	//cout<<"z6"<< pos6.Z()<<endl;
	//cout<<"x9"<< pos9.X()<<endl;
	//cout<<"z9"<< pos9.Z()<<endl;
	//cout<<"x12"<< pos12.X()<<endl;
	//cout<<"z12"<< pos12.Z()<<endl;
	//cout<<"x5"<< pos5.X<<endl;

	//cout<<"c5"<< c5<<endl;
	//cout<<"c12"<< c12<<endl;
	//cout<<"c6"<< c6<<endl;
	//cout<<"c9"<< c9<<endl;

	//cout<<"b5: "<<pos5.X()-c5*(pos5.Z())<<endl;
	//cout<<"b12: "<<pos12.X()-c12*(pos12.Z())<<endl;
	//cout<<"b6: "<<pos6.X()-c6*(pos6.Z())<<endl;
	//cout<<"b9: "<<pos9.X()-c9*(pos9.Z())<<endl;


	int goodmwevents=0; int goodfootevents=0;
	double S_aX_mwpc=0; double S_aX_beam=0;
	double A_aX_mwpc=0; double A_aX_beam=0;
	double SNum_p=0; double SDen1=0; double SDen2=0;
	double P=0; double P_t=0; double SNum_p_t=0;
	double SDen2_t=0; int goodfootevents_t=0;
	double S_aX_beam_t=0; double A_aX_beam_t=0;
	double bX_mwpco=-999.;
	Double_t aX_beam = -999; Double_t bX_beam = -999;
	Double_t aY_beam = -999; Double_t bY_beam = -999;
	Double_t aX_right = -999; Double_t bX_right = -999;
	Double_t aY_right = -999; Double_t bY_right = -999;
	Double_t aX_left = -999; Double_t bX_left = -999;
	Double_t aY_left = -999; Double_t bY_left = -999;
	double Foot6_Z=999.;
	double Foot13_Z=999.;
	double Foot10_Z=999.;
	double Foot7_Z=999.;
	Double_t mw0x = -500.; Double_t mw0y = -500.; Double_t mw1x = -500.; Double_t mw1y = -500.;
	Double_t mw0z = -1083.661789-24; Double_t mw1z = -2250.661789-24;
	Double_t aX_mwpc = -999; Double_t bX_mwpc = -999;
	Double_t aY_mwpc = -999; Double_t bY_mwpc = -999;
	Double_t Q = -999 , Q_ = -999;
	bool p2pt=false;
	bool t1=false;
	Double_t M_x = -999;
	Double_t M_y = -999;
	Double_t x_0 = -999; Double_t y_1 = -999; Double_t y_15 = -999; Double_t x_14 = -999;
	Double_t dist_min1=-999999;
	Double_t dist_min15=-99999;
	//Double_t mw0z = -1523; Double_t mw1z = -2690;
	int count_v=0;
	int detId;
	int stripId;
	double energy;
	int i_init = 0;
	int Mult1t=0;
	int Mult2t=0;
	int Mult3t=0;
	int Mult4t=0;
	int Det_Id=0;
	int Bar_Id=0;
	double EH=0;	
	int Califa_mult;     
	vector <double> ts;	
	vector <double> Et0;
	vector <double> Zt0;
	vector <double> xt0;
	vector <double> yt0;
	vector <double> etat0;
	vector <double> Multt0;
	vector <double> Nst0;
	vector <double> Et1;
	vector <double> xt1;
	vector <double> yt1;
	vector <double> etat1;
	vector <double> Multt1;
	vector <double> Nst1;
	vector <double> Et14;
	vector <double> xt14;
	vector <double> yt14;
	vector <double> etat14;
	vector <double> Multt14;
	vector <double> Nst14;
	vector <double> Et15;
	vector <double> xt15;
	vector <double> yt15;
	vector <double> etat15;
	vector <double> Multt15;
	vector <double> Nst15;
	vector <double> Et10;
	vector <double> zt10;
	vector <double> xt10;
	vector <double> etat10;
	vector <double> Multt10;
	vector <double> Nst10;
	vector <double> Et5;
	vector <double> yt5;
	vector <double> etat5;
	vector <double> Multt5;
	vector <double> Nst5;
	vector <double> Et3;
	vector <double> xt3;
	vector <double> zt3;
	vector <double> etat3;
	vector <double> Multt3;
	vector <double> Nst3;
	vector <double> Et12;
	vector <double> yt12;
	vector <double> etat12;
	vector <double> Multt12;
	vector <double> Nst12;
	vector <double> Et8;
	vector <double> zt8;
	vector <double> xt8;
	vector <double> etat8;
	vector <double> Multt8;
	vector <double> Nst8;
	vector <double> Et11;
	vector <double> xt11;
	vector <double> zt11;
	vector <double> etat11;
	vector <double> Multt11;
	vector <double> Nst11;
	vector <double> Et9;
	vector <double> yt9;
	vector <double> etat9;
	vector <double> Multt9;
	vector <double> Nst9;
	vector <double> Et6;
	vector <double> yt6;
	vector <double> etat6;
	vector <double> Multt6;
	vector <double> Nst6;
	vector <double> lost;
	vector <int> Multe0;
	vector <int> Multe1;
	vector <int> Multe14;
	vector <int> Multe15;
	vector <double> Et0r;
	vector <double> Et1r;
	vector <double> zt0;
	vector <double> zt1;
	vector <double> zt14;
	vector <double> zt15;
	vector <int> S_id;
	vector <int> D_id;
	vector <double> hit_energy;
	vector <double> hit_theta;
	vector <double> hit_phi;
	vector <int> nbofcrystalhits;
	int multhits;
	vector <double> frs_aq;
	vector <double> frs_z;
	vector <double> frs_b;
	vector <double> ns;
	vector <double> nf;
	vector <int> c_id;
	vector <int> c_t;
	vector <double> Zt;
	vector <double> mw0xt;
	vector <double> mw0yt;
	vector <double> mw1xt;
	vector <double> mw1yt;
	vector <int> Bar_Idt;
	vector <int> Det_Idt;
	vector <int> Mult1;
	vector <int> Mult2;
	vector <int> Mult3;
	vector <int> Mult4;
	vector <double> theta1;
	vector <double> theta2;
	vector <double> phi1;
	vector <double> phi2;
	vector <double> opa_c;
	vector <double> opa_f;
	vector <double> opa_f_t;
	vector <double> e1c;
	vector <double> e2c;
	vector <double> dphi_f;
	vector <double> phi_fl;
	vector <double> phi_fr;
	vector <double> phi_fps;
	vector <double> foot8_x;
	vector <double> foot11_x;
	vector <double> Et0n;
	vector <double> Et1n;
	vector <double> Et14n;
	vector <double> Et15n;
	vector <double> pro01;
	vector <double> pro45;
	vector <double> Et0new;
	vector <double> Et1new;
	vector <double> Et14new;
	vector <double> Et15new;
	vector <double> AX_MW;
	vector <double> BX_MW;
	vector <double> AY_MW;
	vector <double> BY_MW;
	vector <double> AX_beam;
	vector <double> AY_beam;
	vector <double> BX_beam;
	vector <double> BY_beam;
	vector <double> z_v;
	vector <double> x_v;
	vector <double> y_v;
	vector <double> z_vm_l;
	vector <double> x_vm_l;
	vector <double> y_vm_l;
	vector <double> z_vm_lb;
	vector <double> x_vm_lb;
	vector <double> y_vm_lb;
	vector <double> d;
	vector <double> dm_l;
	vector <double> t_pl;
	vector <double> t_pr;
	vector <double> t_ps;
	vector <double> AX_right;
	vector <double> AY_right;
	vector <double> BX_right;
	vector <double> BY_right;
	vector <double> AX_left;
	vector <double> AY_left;
	vector <double> BX_left;
	vector <double> BY_left;
	vector <double> zt5;
	vector <double> zt12;
	vector <double> zt6;
	vector <double> zt9;
	vector <double> zv_t;
	vector <double> yv_t;
	vector <double> xv_t;
	vector <double> z_vm_r;
	vector <double> x_vm_r;
	vector <double> y_vm_r;
	vector <double> dm_r;
	vector <double> theta_rv;
	vector <double> theta_lv;
	vector <double> phi_rv;
	vector <double> phi_lv;
	vector <double> z_vm_t;
	vector <double> x_vm_t;
	vector <double> y_vm_t;
	vector <double> dm_t;
	vector <double> z_vm_tt;
	vector <double> x_vm_tt;
	vector <double> y_vm_tt;
	vector <double> dm_tt;


	//Starting new variables 
	vector <double> ml;
	vector <double> mr;
	vector <int> eventn;
	vector <double> theta_p_l;
	vector <double> theta_p_r;
	vector <double> phi_p_l;
	vector <double> phi_p_r;
	vector <double> theta_p_l_m;
	vector <double> theta_p_r_m;
	vector <double> phi_p_l_m;
	vector <double> phi_p_r_m;
	vector <double> theta_p_l_m_t;
	vector <double> theta_p_r_m_t;
	vector <double> phi_p_l_m_t;
	vector <double> phi_p_r_m_t;
	vector <double> theta_rn;
	vector <double> theta_ln;
	vector <double> phi_ln;
	vector <double> phi_rn;
	vector <double> opa_rl;
	vector <double> opa_lr;
	vector <double> opa_t;
	vector <double>x_extt;
	vector <double> dpl;
	vector <double> dpr;
	vector <double> dtl;
	vector <double> dtr;
	vector <double> Mm;
	vector <double> Pmt;
	vector <double> Pmtf;
	vector <double> Pmx;
	vector <double> Pmy;
	vector <double> Pmz;
	vector <double> Em;
	vector <double> Erc;
	vector <double> Elc;
	vector <double> Mmf;
	vector <double> Pmxf;
	vector <double> Pmyf;
	vector <double> Pmzf;
	vector <double> Emf;
	TLorentzVector p3f, p4f,pmissf;
	TLorentzVector p3, p4, ptg, pmiss, pbeam;
	double Emiss, mmiss, Emissf, mmissf;
	bool p2pt_n=false;
	int mult0n=0;
	int mult1n=0;
	int mult14n=0;
	int mult15n=0;
	int mult0n_n=0;
	int mult1n_n=0;
	int mult14n_n=0;
	int mult15n_n=0;
	bool t2=false;
	bool t7 = false;
	bool t8 = false;
	bool t13 = false;
	bool t14 = false;
	bool t15 = false;
	bool t16 = false;
	bool t3 = false;
	bool t4=false;
	Int_t TPATn;

	//New File definition
	TFile *q = new TFile("vertex_0001_ts.root","RECREATE");
	TTree *vertex = new TTree ("vertex","vertex"); 

	vertex->Branch("FRS_Aq",&frs_aq);
	vertex->Branch("FRS_Z",&frs_z);
	vertex->Branch("FRS_B",&frs_b);
	vertex->Branch("FOOT3_E",&Et3);
	vertex->Branch("Nf",&nf);
	vertex->Branch("Ns",&ns);
	vertex->Branch("Crystal_ID",&c_id);
	vertex->Branch("p2p",&p2pt_n);
	vertex->Branch("TPAT",&TPATn);
	vertex->Branch("T1",&t1);
	vertex->Branch("T2",&t2);
	vertex->Branch("T3",&t3);
	vertex->Branch("T4",&t4);
	vertex->Branch("Strip_Id",&S_id);
	vertex->Branch("Det_Id",&D_id);
	vertex->Branch("TOFD_Z",&Zt);
	vertex->Branch("TOFD_Mult1",&Mult1);
	vertex->Branch("TOFD_Mult2",&Mult2);
	vertex->Branch("TOFD_Mult3",&Mult3);
	vertex->Branch("TOFD_Mult4",&Mult4);
	vertex->Branch("TOFD_DetId",&Det_Idt);
	vertex->Branch("TOFD_BarId",&Bar_Idt);
	vertex->Branch("Z_vertex_m_t",&z_vm_t);
	vertex->Branch("Y_vertex_m_t",&y_vm_t);
	vertex->Branch("X_vertex_m_t",&x_vm_t);
	vertex->Branch("Dist_m_t",&dm_t);
	vertex->Branch("Theta_Califa_r",&theta_rv);
	vertex->Branch("Theta_Califa_l",&theta_lv);
	vertex->Branch("Phi_c_r",&phi_rv);
	vertex->Branch("Phi_c_l",&phi_lv);
	vertex->Branch("Energy_Califa_r",&Erc);
	vertex->Branch("Energy_Califa_l",&Elc);
	vertex->Branch("Opening_angle_Califa",&opa_c);
	vertex->Branch("Opening_angle_Foot_t",&opa_f_t);
	vertex->Branch("E_loss_p",&hit_energy);
	vertex->Branch("NbOfCrystalHits",&nbofcrystalhits);
	vertex->Branch("NbHits",&multhits);
	vertex->Branch("Theta_proton_left",&t_pl);
	vertex->Branch("Theta_proton_right",&t_pr);
	vertex->Branch("Theta_protons",&t_ps);
	vertex->Branch("Mult_left",&ml);
	vertex->Branch("Mult_right",&mr);
	vertex->Branch("Event",&eventn);
	vertex->Branch("Theta_p_l",&theta_p_l_m_t);
	vertex->Branch("Theta_p_r",&theta_p_r_m_t);
	vertex->Branch("Phi_p_l",&phi_p_l_m_t);
	vertex->Branch("Phi_p_r",&phi_p_r_m_t);
	vertex->Branch("Opening_angle_t",&opa_t);
	vertex->Branch("Dpl",&dpl);
	vertex->Branch("Dpr",&dpr);
	vertex->Branch("Dtr",&dtr);
	vertex->Branch("Dtl",&dtl);
	vertex->Branch("Pt",&Pmt);
	vertex->Branch("Ptf",&Pmtf);
	vertex->Branch("Mmiss",&Mm);
	vertex->Branch("Emiss",&Em);
	vertex->Branch("Px",&Pmx);
	vertex->Branch("Py",&Pmy);
	vertex->Branch("Pz",&Pmz);
	vertex->Branch("Mmissf",&Mmf);
	vertex->Branch("Emissf",&Emf);
	vertex->Branch("Pxf",&Pmxf);
	vertex->Branch("Pyf",&Pmyf);
	vertex->Branch("Pzf",&Pmzf);
	vertex->Branch("TimeStamp",&ts);
	//nevents=100000;

	for(int i = 0; i<nevents ; i++){

		//cout<<"i-> "<<i<<endl;
		TPATn=TPAT;
		p2pt_n=p2pt;
		eventn.clear();
		eventn.push_back(i);

		//bool p2pt=false;
		bool p2pt_n=false;
		int mult0n=0;
		int mult1n=0;
		int mult14n=0;
		int mult15n=0;
		int mult0n_n=0;
		int mult1n_n=0;
		int mult14n_n=0;
		int mult15n_n=0;

		std::clog << "Analysis Info    :   " << (i-i_init)*100/(nevents-i_init) << " % of events treated, " << count_v << " Good vertex found " << "\r";
		//std::clog << "Analysis Info    :   " << count_v << " Good vertex found" << "\r";
		//cout<<"TPAT--> "<<TPAT<<endl;
		t4=false;

		t1 = (TRIG(1)) & TPAT;
		t2 = (TRIG(2)) & TPAT;
		t3 = (TRIG(3)) & TPAT;
		t4 = (TRIG(4)) & TPAT;
		t7 = (TRIG(7)) & TPAT;
		t8 = (TRIG(8)) & TPAT;
		t13 = (TRIG(13)) & TPAT;
		t14 = (TRIG(14)) & TPAT;
		t15 = (TRIG(15)) & TPAT;
		t16 = (TRIG(16)) & TPAT;

		//clear variables inside code
		Mult1.clear();
		Mult2.clear();
		Mult3.clear();
		Mult4.clear();
		Zt.clear();
		c_t.clear();
		Bar_Idt.clear();
		Det_Idt.clear();
		mw0xt.clear();
		mw0yt.clear();
		mw1xt.clear();
		mw1yt.clear();
		Et0.clear();
		xt0.clear();
		yt0.clear();
		zt0.clear();
		zt1.clear();
		zt14.clear();
		zt15.clear();
		etat0.clear();
		Multt0.clear();
		Nst0.clear();
		Et1.clear();
		xt1.clear();
		yt1.clear();
		etat1.clear();
		Multt1.clear();
		Nst1.clear();
		Et14.clear();
		xt14.clear();
		yt14.clear();
		etat14.clear();
		Multt14.clear();
		Nst14.clear();
		Et15.clear();
		xt15.clear();
		yt15.clear();
		etat15.clear();
		Multt15.clear();
		Nst15.clear();
		Multe0.clear();
		Multe1.clear();
		Multe14.clear();
		Multe15.clear();
		p2pt=false;
		Et0r.clear();
		Et1r.clear();
		Zt0.clear();
		S_id.clear();
		D_id.clear();
		frs_b.clear();
		frs_aq.clear();
		frs_z.clear();
		Et8.clear();
		xt8.clear();
		zt8.clear();
		zt11.clear();
		zt3.clear();
		zt10.clear();
		etat8.clear();
		Multt8.clear();
		Et11.clear();
		xt11.clear();
		etat11.clear();
		Multt11.clear();
		Et9.clear();
		yt9.clear();
		etat9.clear();
		Multt9.clear();
		Et6.clear();
		yt6.clear();
		etat6.clear();
		Multt6.clear();
		Et10.clear();
		xt10.clear();
		etat10.clear();
		Multt10.clear();
		Et3.clear();
		xt3.clear();
		etat3.clear();
		Multt3.clear();
		Et12.clear();
		yt12.clear();
		etat12.clear();
		Multt12.clear();
		Et5.clear();
		yt5.clear();
		etat5.clear();
		Multt5.clear();
		lost.clear();
		ns.clear();
		nf.clear();
		c_id.clear();
		hit_energy.clear();
		hit_theta.clear();
		hit_phi.clear();
		nbofcrystalhits.clear();
		multhits=-1;
		Nst5.clear();
		Nst12.clear();
		Nst3.clear();	
		Nst10.clear();
		Nst11.clear();
		Nst8.clear();
		Nst6.clear();
		Nst9.clear();

		// Clear of new variables
		phi_rv.clear();
		phi_lv.clear();
		theta_rv.clear();
		theta_lv.clear();
		theta_p_l_m_t.clear();
		theta_p_r_m_t.clear();
		phi_p_l_m_t.clear();
		phi_p_r_m_t.clear();
		theta_p_l.clear();
		theta_p_r.clear();
		phi_p_l.clear();
		phi_p_r.clear();
		opa_t.clear();
		AX_MW.clear();
		AY_MW.clear();
		BX_MW.clear();
		BY_MW.clear();
		AX_beam.clear();
		AY_beam.clear();
		BX_beam.clear();
		BY_beam.clear();
		phi1.clear();
		phi2.clear();
		theta1.clear();
		theta2.clear();
		opa_c.clear();
		z_vm_t.clear();
		x_vm_t.clear();
		y_vm_t.clear();	
		dm_t.clear();
		t_pl.clear();
		t_pr.clear();
		t_ps.clear();
		AX_right.clear();
		AY_right.clear();
		BX_right.clear();
		BY_right.clear();
		AX_left.clear();
		AY_left.clear();
		BX_left.clear();
		BY_left.clear();
		zt5.clear();
		zt12.clear();
		zt6.clear();
		zt9.clear();
		ml.clear();
		mr.clear();
		opa_f.clear();
		opa_f_t.clear();
		dtl.clear();
		dtr.clear();
		dpl.clear();
		dpr.clear();
		Mm.clear();
		Pmx.clear();
		Pmy.clear();
		Pmz.clear();
		Em.clear();
		Erc.clear();
		Elc.clear();
		Mmf.clear();
		Pmxf.clear();
		Pmyf.clear();
		Pmzf.clear();
		Emf.clear();
		e1c.clear();
		e2c.clear();
		Pmt.clear();
		Pmtf.clear();

		hitDataFRS->Clear();
		hitDataMWPC0->Clear();
		hitDataMWPC1->Clear();
		TofdHitData->Clear();
		LosTCalData->Clear();
		FootHitData->Clear();	
		CalifaClusterData->Clear();

		double x_mwt=-999.;
		double y_mwt=-999.;
		double x_ft=-999.;
		double Phi1=999.;
		double Phi2=999.;
		double Theta1=999;
		double Theta2=999.;
		double Energy1=999.;
		double Energy2=999.;
		double openingAngle=999.;
		tree->GetEntry(i);

		//For Loop Over Array Size want all to be same size for hit level -> "one event in an event"
		Int_t CalifaHitsPerEvent = CalifaClusterData->GetEntriesFast();
		Int_t frsHitsPerEvent = hitDataFRS->GetEntriesFast();
		Int_t mwpc0HitsPerEvent = hitDataMWPC0->GetEntriesFast();
		Int_t mwpC1HitsPerEvent= hitDataMWPC1->GetEntriesFast();
		Int_t LoshitPerEvent = LosTCalData->GetEntriesFast();
		Int_t tofdHitsPerEvent = TofdHitData->GetEntriesFast();
		Int_t footHitsPerEvent = FootHitData->GetEntriesFast();

		// bool p2p = (TRIG(4) | TRIG(10)) & DataCA->GetTpat();
		p2pt = (TRIG(4)) & DataCA->GetTpat();
		t1 = (TRIG(1)) & DataCA->GetTpat();
		TPAT = DataCA->GetTpat();
		EH = DataCA->GetTimeStamp();
		ts.push_back(EH);
		R3BFrsData**frsData = new R3BFrsData*[frsHitsPerEvent];
		R3BCalifaClusterData** califaData = new R3BCalifaClusterData*[CalifaHitsPerEvent];
		R3BMwpcHitData** mwpc0Data = new R3BMwpcHitData*[mwpc0HitsPerEvent];
		R3BMwpcHitData** mwpC1Data = new R3BMwpcHitData*[mwpC1HitsPerEvent];
		R3BTofdHitData** tofdData = new R3BTofdHitData*[tofdHitsPerEvent];
		R3BLosTCalData** LosData = new R3BLosTCalData*[LoshitPerEvent];
		R3BFootHitData** foothData = new R3BFootHitData*[footHitsPerEvent];
		Double_t Energy[1000];
		Double_t Phi[1000];
		Double_t Theta[1000];
		Double_t crystalmult[1000];
		Double_t NS[1000];
		Double_t NF[1000];
		Int_t cid;
		Int_t C_T[1000];
		multhits = CalifaHitsPerEvent;
		int mul_cluster_califa = multhits;
		for (Int_t j = 0; j<CalifaHitsPerEvent; j++){
			califaData[j] = static_cast<R3BCalifaClusterData*>(CalifaClusterData->At(j));		
			Energy[j] = califaData[j]->GetEnergy();
			Phi[j] = TMath::RadToDeg()* califaData[j]->GetPhi();
			Theta[j] = TMath::RadToDeg()* califaData[j]->GetTheta();
			NS[j] = califaData[j]->GetNs();
			NF[j] = califaData[j]->GetNf();
			ns.push_back(NS[j]);
			nf.push_back(NF[j]);
			crystalmult[j] = califaData[j]->GetNbOfCrystalHits();
			hit_energy.push_back(Energy[j]);
			hit_theta.push_back(Theta[j]);		
			hit_phi.push_back(Phi[j]);
			nbofcrystalhits.push_back(crystalmult[j]);
			C_T[j] = califaData[j]->GetClusterType();	
			cid = ((R3BCalifaClusterData*)califaData[j])->GetCrystalList().at(0);
			c_id.push_back(cid);
			c_t.push_back(C_T[j]);	
		}

		if(tofdHitsPerEvent>0){
			for(Int_t j = 0 ; j<tofdHitsPerEvent;j++){			
				tofdData[j] = static_cast<R3BTofdHitData*>(TofdHitData->At(j));
				Bar_Id = tofdData[j]->GetBarId();
				Det_Id = tofdData[j]->GetDetId();
				Q = tofdData[j]->GetEloss();
				Det_Idt.push_back(Det_Id);
				Bar_Idt.push_back(Bar_Id);
				Zt.push_back(Q);
				if(Det_Id==1){Mult1t++;}
				if(Det_Id==2){Mult2t++;}
				if(Det_Id==3){Mult3t++;}
				if(Det_Id==4){Mult4t++;}

			}
		}

		int Plane = 0;
		if(frsHitsPerEvent>0){
			for(Int_t j = 0 ;j<frsHitsPerEvent;j++){
				frsData[j] = static_cast<R3BFrsData*>(hitDataFRS->At(j));
				FRSAQ  = frsData[j]->GetAq();
				FRSZ  = frsData[j]->GetZ();
				FRSB  = frsData[j]->GetBeta();
				FRSBR  = frsData[j]->GetBrho();
				FRST = frsData[j]->GetTof();
				FRSX = frsData[j]->GetXS2();
				FRSM = FRSAQ*FRSZ;
				frs_aq.push_back(FRSAQ);
				frs_z.push_back(FRSZ);
				frs_b.push_back(FRSB);
			}
		}
		if(footHitsPerEvent>0){
			for (Int_t i = 0; i < footHitsPerEvent; i++){
				foothData[i]=static_cast<R3BFootHitData*>(FootHitData->At(i));
				detId = foothData[i]->GetDetId() - 1;
				//stripId = foothData[i]->GetStripId() - 1;
				energy = foothData[i]->GetEnergy();
				D_id.push_back(detId);
				if(detId==0){  
					xt0.push_back(foothData[i]->GetPosLab().X()); 
					zt0.push_back(foothData[i]->GetPosLab().Z());
					Et0.push_back(foothData[i]->GetEnergy());
				}
				if(detId==14){  
					xt14.push_back(foothData[i]->GetPosLab().X());
					zt14.push_back(foothData[i]->GetPosLab().Z());
					Et14.push_back(foothData[i]->GetEnergy());
				}

				if(detId==8){  
					xt8.push_back(foothData[i]->GetPosLab().X());
					zt8.push_back(foothData[i]->GetPosLab().Z());
					Et8.push_back(foothData[i]->GetEnergy());
				}
				if(detId==11){  
					xt11.push_back(foothData[i]->GetPosLab().X()); 
					zt11.push_back(foothData[i]->GetPosLab().Z());
					Et11.push_back(foothData[i]->GetEnergy());
				}
				if(detId==10){ 
					xt10.push_back(foothData[i]->GetPosLab().X()); 
					zt10.push_back(foothData[i]->GetPosLab().Z());
					Et10.push_back(foothData[i]->GetEnergy());
				}
				if(detId==3){  
					xt3.push_back(foothData[i]->GetPosLab().X()); 
					zt3.push_back(foothData[i]->GetPosLab().Z());
					Et3.push_back(foothData[i]->GetEnergy());
				}
				if(detId==1){  
					yt1.push_back(foothData[i]->GetPosLab().Y()); 
					Et1.push_back(foothData[i]->GetEnergy());
				}
				if(detId==15){  
					yt15.push_back(foothData[i]->GetPosLab().Y()); 
					Et15.push_back(foothData[i]->GetEnergy());
				}
				if(detId==5){  
					yt5.push_back(foothData[i]->GetPosLab().Y()); 
					Et5.push_back(foothData[i]->GetEnergy());
				}
				if(detId==12){  
					yt12.push_back(foothData[i]->GetPosLab().Y()); 
					Et12.push_back(foothData[i]->GetEnergy());
				}
				if(detId==6){  
					yt6.push_back(foothData[i]->GetPosLab().Y()); 
					Et6.push_back(foothData[i]->GetEnergy());
				}
				if(detId==9){  
					yt9.push_back(foothData[i]->GetPosLab().Y()); 
					Et9.push_back(foothData[i]->GetEnergy());
				}
			}
		}



		int count_los=0;
		if(LoshitPerEvent>0){
			for(int l=0; l<LoshitPerEvent; l++){

				LosData[l]=static_cast<R3BLosTCalData*>(LosTCalData->At(l));
				count_los++;
			}

		}
		lost.push_back(count_los);

		for (Int_t j = 0; j<mwpc0HitsPerEvent; j++){
			mwpc0Data[j]=static_cast<R3BMwpcHitData*>(hitDataMWPC0->At(j));
			mw0x = mwpc0Data[j]->GetX() -0.6;
			mw0y = mwpc0Data[j]->GetY();
			// MW0_spot->Fill(mw0x,mw0y);
			mw0xt.push_back(mw0x);
			mw0yt.push_back(mw0y);
		}
		for (Int_t j = 0; j<mwpC1HitsPerEvent; j++){

			mwpC1Data[j]=static_cast<R3BMwpcHitData*>(hitDataMWPC1->At(j));
			mw1x = mwpC1Data[j]->GetX();
			mw1y = mwpC1Data[j]->GetY() + 17.0;
			//MW_spot->Fill(mw1x,mw1y);
			mw1xt.push_back(mw1x);
			mw1yt.push_back(mw1y);
		}

		if(TPAT<4000 /*&& t4*/){

			//MWPC TRACKING
			double mw0xn=-999;
			double mw0yn=-999;
			double mw1xn=-999;
			double mw1yn=-999;
			double xmw0=-999;
			double xmw1=-999;
			double frsaq=-999;
			double frsz=-999;
			if(frs_z.size()>0 && frs_aq.size()>0){
				frsz=frs_z.at(0);
				frsaq=frs_aq.at(0);

			}			
			double Zmin = 1. , Zmax = 8. , AQmin = 1. , AQmax = 4. ;//all

			//double Zmin = 5.4 , Zmax = 6.6 , AQmin = 2.675 , AQmax = 2.694 ;//16C s522
			//double Zmin = 4.5 , Zmax = 5.5 , AQmin = 4.6 , AQmax = 4.65 ;//B s509
			//double Zmin = 5.6 , Zmax = 6.2 , AQmin = 2.740 , AQmax = 2.755 ;//18C s509
			//double Zmin = 6.3 , Zmax = 7.4 , AQmin = 2.715 , AQmax = 2.730 ;//19N s509
			if(frsz<Zmax && frsz>Zmin && frsaq<AQmax && frsaq>AQmin){
				//CALIFA CrystalID To Position xyx
				//CALIFA ID_TO_POSITION Variables
				Double_t CALIFA_IDtoX[2432];
				Double_t CALIFA_IDtoY[2432];
				Double_t CALIFA_IDtoZ[2432];

				//Load CALIFA Map
				ifstream infile("/feynman/home/dphn/lena/al268379/work/s522/foot/vertex/CALIFA_Mapping.dat");
				Int_t a;
				Double_t b, c, d;

				while (infile >> a >> b >> c >> d){
					//cout << a << "   " << b << "   " << c << "   " << d << endl; 
					CALIFA_IDtoX[a-1] = b;
					CALIFA_IDtoY[a-1] = c;
					CALIFA_IDtoZ[a-1] = d;
					//cout << a << "   " << CALIFA_IDtoX[a-1] << "   " << CALIFA_IDtoY[a-1] << "   " << CALIFA_IDtoZ[a-1] << endl; 
				}

				infile.close();

				//if(mw0x>-500. && mw0y>-500. && mw1xt>-500. && mw1yt>-500){

				TVector3 master_califa[2];
				double maxE1 = 0., maxE2 = 0.;
				int crystalhits1 = 0, crystalhits2 = 0;

				if(mul_cluster_califa>1 && hit_energy.size()>1 && hit_theta.size()>1 && hit_phi.size()>1){
					maxE1 = hit_energy.at(0);
					master_califa[0].SetMagThetaPhi(1.,hit_theta.at(0)*TMath::DegToRad(),hit_phi.at(0)*TMath::DegToRad());
					maxE2 = hit_energy.at(1);
					master_califa[1].SetMagThetaPhi(1.,hit_theta.at(1)*TMath::DegToRad(),hit_phi.at(1)*TMath::DegToRad());
				}

				//////////////////////////////
				if(mul_cluster_califa>1){

					Energy1 = maxE1;
					Energy2 = maxE2;
					e1c.push_back(Energy1);
					e2c.push_back(Energy2);
					Phi1 = master_califa[0].Phi()*TMath::RadToDeg();
					Phi2 = master_califa[1].Phi()*TMath::RadToDeg();
					phi1.push_back(Phi1);
					phi2.push_back(Phi2);

				}

				if((mul_cluster_califa>1 && TMath::Abs(Phi1-Phi2)>150. && TMath::Abs(Phi1-Phi2)<210. && Energy1>40000. && Energy2 > 40000. )){

					if(mul_cluster_califa>1){

						Theta1 = master_califa[0].Theta()*TMath::RadToDeg();
						Theta2 = master_califa[1].Theta()*TMath::RadToDeg();
						theta1.push_back(Theta1);
						theta2.push_back(Theta2);
						openingAngle = TMath::Sin(TMath::DegToRad()*Theta1)*TMath::Sin(TMath::DegToRad()*Theta2)*TMath::Cos(TMath::DegToRad()*Phi2-TMath::DegToRad()*Phi1) + TMath::Cos(TMath::DegToRad()*Theta1)*TMath::Cos(TMath::DegToRad()*Theta2);
						openingAngle = TMath::RadToDeg()*TMath::ACos(openingAngle);
						opa_c.push_back(openingAngle);

					}

					int c_id_l=-999.;
					int c_id_r=-999.;
					double theta_l=-999.;
					double theta_r=-999.;
					double phi_r=-999.;
					double phi_l=-999.;
					double El=-999;
					double Er=-999;

					if(TMath::Abs(Phi1)>90. && TMath::Abs(Phi2)<90. && Theta1>0 && Theta2>0 && c_id.size()>1){
						theta_r=Theta1;
						theta_l=Theta2;
						Er=Energy1;
						El=Energy2;
						Erc.push_back(Er);
						Elc.push_back(El);
						theta_rv.push_back(theta_r);
						theta_lv.push_back(theta_l);
						phi_r=Phi1;
						phi_l=Phi2;
						phi_rv.push_back(phi_r);
						phi_lv.push_back(phi_l);
						c_id_r=c_id.at(0);
						c_id_l=c_id.at(1);
					}

					if(TMath::Abs(Phi1)<90. && TMath::Abs(Phi2)>90. && Theta1>0 && Theta2>0 && c_id.size()>1){
						theta_r=Theta2;
						theta_l=Theta1;
						Er=Energy2;
						El=Energy1;	
						Erc.push_back(Er);
						Elc.push_back(El);	
						theta_rv.push_back(theta_r);
						theta_lv.push_back(theta_l);
						c_id_r=c_id.at(1);
						c_id_l=c_id.at(0);
						phi_r=Phi2;
						phi_l=Phi1;
						phi_rv.push_back(phi_r);
						phi_lv.push_back(phi_l);

					}

					for(int a=0; a<mw0xt.size(); a++){
						for(int b=0; b<mw1xt.size(); b++){	

							mw0xn=mw0xt.at(a);
							mw1xn=mw1xt.at(b);
							aX_mwpc = (-mw0xn + mw1xn)/(mw0z - mw1z);
							bX_mwpc = -mw1xn - aX_mwpc*(mw1z);
							AX_MW.push_back(aX_mwpc);
							BX_MW.push_back(bX_mwpc);

						}
					}

					for(int c=0; c<mw0yt.size(); c++){
						for(int d=0; d<mw1yt.size(); d++){

							mw0yn=mw0yt.at(c);
							mw1yn=mw1yt.at(d);
							aY_mwpc = (mw0yn -  mw1yn)/(mw0z - mw1z);
							bY_mwpc = mw1yn - aY_mwpc*(mw1z);
							AY_MW.push_back(aY_mwpc);
							BY_MW.push_back(bY_mwpc);

						}
					}

					if(aX_mwpc*441+bX_mwpc>-7){
						double Mult_beam=-999.;

						if(xt0.size()>0 && yt1.size()>0 && xt14.size()>0 && yt15.size()>0 && zt14.size()>0 && zt0.size()>0 && zt1.size()>0 && zt15.size()>0){
							for(int n=0; n<xt0.size(); n++){
								for(int m=0; m<xt14.size();m++){
									for(int r=0; r<yt1.size(); r++){
										for(int h=0; h<yt15.size();h++){

											aX_beam = ((xt14.at(m))-(xt0.at(n)))/(zt14.at(m) - zt0.at(n));
											bX_beam = (xt0.at(n)+sx)-aX_beam*zt0.at(n);
											AX_beam.push_back(aX_beam);
											BX_beam.push_back(bX_beam);	
											Mult_beam++;
											aY_beam = ((yt15.at(h))-(yt1.at(r)))/(zt15.at(h) - zt1.at(r));
											bY_beam = (yt1.at(r)+sy)-aY_beam*zt1.at(r);
											AY_beam.push_back(aY_beam);
											BY_beam.push_back(bY_beam);

										}
									}
								}
							}
						}

						int Count_right=0.;
						int Count_right_x=0;
						int Count_right_y=0;
						int Mult_left=0;
						int Mult_right=0;
						int Mult_right_x=0.;
						int Mult_right_y=0;
						int Mult_left_x=0.;
						int Mult_left_y=0;
						TVector3 vprot_right;
						TVector3 vprot_left;
						double xr10=0;
						double xr3=0;
						double xr10_m=0;
						double xr3_m=0;
						double theta_p_rf=999.;
						double phi_p_rf=999.;	
						double theta_p_lf=999.;
						double phi_p_lf=999.;	
						double opa_foot=999.;
						double opa_foot_t=999.;
						double theta_la=-999.;
						double theta_ra=-999.;
						double phi_la=-999.;
						double phi_ra=-999.;
						TVector3 pra,pla;

						if(xt10.size()>0 && xt3.size()>0 && yt5.size()>0 && yt12.size()>0){
							for(int i=0; i<xt10.size(); i++){
								for(int j=0; j<xt3.size();j++){
									if(cutg3->IsInside(xt10.at(i),xt3.at(j))){
										Count_right_x++;
										for(int r=0; r<yt5.size(); r++){
											for(int h=0; h<yt12.size();h++){
												if(cutg5->IsInside(yt5.at(r),yt12.at(h))){
													//Try cndition on proton energy
													//if(Et10.at(i)<15 && Et3.at(j)<15 && Et5.at(r)<25 && Et12.at(h)<25){
													aX_right = (xt10.at(i)-xt3.at(j))/( zt10.at(i)- zt3.at(j));
													bX_right = (xt3.at(j)+sx)-aX_right*(zt3.at(j));
													if(aX_right>0.255 && aX_right<1.57){	
														AX_right.push_back(aX_right);
														BX_right.push_back(bX_right);
														theta_p_r.push_back(theta_ra);
														phi_p_r.push_back(phi_ra);
														Foot6_Z = (bX_right-fBDet[5])/(fADet[5]-aX_right);
														Foot13_Z = (bX_right-fBDet[12])/(fADet[12]-aX_right);
														if((Foot6_Z>474.6 && Foot6_Z<571.6) && (Foot13_Z<533.4 && Foot13_Z>454.9) && aX_right>0){
															zt5.push_back(Foot6_Z);
															zt12.push_back(Foot13_Z);	
															aY_right = (yt5.at(r)-yt12.at(h))/( Foot6_Z - Foot13_Z);
															bY_right = (yt12.at(h)+sy)-aY_right*Foot13_Z;
															//Insert condition aY!
															//if(abs(aY_right)<0.69){
															AY_right.push_back(aY_right);
															BY_right.push_back(bY_right);
															xr10=aX_right*400+bX_right;
															xr3=aX_right*600+bX_right;	
															TLine *l_r=new TLine (400,xr10,600,xr3);

															TVector3 Posr_1a(aX_right*Foot6_Z+bX_right,aY_right*Foot6_Z+bY_right,Foot6_Z);
															TVector3 Posr_2a(aX_right*Foot13_Z+bX_right,aY_right*Foot13_Z+bY_right,Foot13_Z);
															pra=Posr_1a-Posr_2a;

															//pr=Posr_1-Posr;
															//pl=Posl_1-Posl;       
															theta_ra=pra.Theta();
															phi_ra=pra.Phi();
															Count_right_y++;
															//} phi cond
														}//cuttr
													}//cutcinq
													//}cut_energy
												}
											}

										}
									}
								}

							}
						}	
						Mult_right_x = Count_right_x;
						Mult_right_y = Count_right_y;
						//vector multiplicity
						int Count_left_x=0;
						int Count_left_y=0;
						double xl11=0;
						double xl8=0;
						double xl11_m=0;
						double xl8_m=0;
						double xl11_ml=0;
						double xl8_ml=0;
						double xr3_ml=0;
						double xr10_ml=0;
						if(xt8.size()>0 && xt11.size()>0 && yt9.size()>0 && yt6.size()>0){
							for(int i=0; i<xt8.size(); i++){
								for(int j=0; j<xt11.size();j++){
									if(cutg8->IsInside(xt8.at(i),xt11.at(j))){
										Count_left_x++;
										for(int r=0; r<yt9.size(); r++){
											for(int h=0; h<yt6.size();h++){
												if(cutg6->IsInside(yt9.at(r),yt6.at(h))){
													//if(Et8.at(i)<15 && Et11.at(j)<15 && Et9.at(r)<25 && Et6.at(h)<25){
														aX_left = (xt8.at(i)-xt11.at(j))/( zt8.at(i)- zt11.at(j));
														bX_left = (xt11.at(j)+sx)-aX_left*(zt11.at(j));
														if(aX_left<-0.255 && aX_left>-1.57){	
															AX_left.push_back(aX_left);
															BX_left.push_back(bX_left);
															theta_p_l.push_back(theta_la);
															phi_p_l.push_back(phi_la);
															Foot10_Z = (bX_left-fBDet[9])/(fADet[9]-aX_left);
															Foot7_Z = (bX_left-fBDet[6])/(fADet[6]-aX_left);
															if((Foot7_Z>440.34 && Foot7_Z<548.9) && (Foot10_Z<571.3 && Foot10_Z>482.6) && aX_left<0){
																aY_left = (yt9.at(r)-yt6.at(h))/( Foot10_Z - Foot7_Z);
																bY_left = (yt6.at(h)+sy)-aY_left*Foot7_Z;
																//if(abs(aY_left)<0.69){
																	AY_left.push_back(aY_left);
																	BY_left.push_back(bY_left);
																	TVector3 Posl_1a(aX_left*Foot7_Z+bX_left,aY_left*Foot7_Z+bY_left,Foot7_Z);
																	TVector3 Posl_2a(aX_left*Foot10_Z+bX_left,aY_left*Foot10_Z+bY_left,Foot10_Z);
																	pla=Posl_2a-Posl_1a;
																	//pr=Posr_1-Posr;
																	//pl=Posl_1-Posl;       
																	theta_la=pla.Theta();
																	phi_la=pla.Phi();
																	Count_left_y++;
																//} phi cond
															}//cut
														}//cut
													//}cut_energy
												}
											}

										}
									}
								}
							}
						}

						Mult_left_x = Count_left_x;
						Mult_left_y = Count_left_y;
						//vector multiplicity

						double openingAngle_f=999.;
						opa_foot = TMath::Sin(theta_ra)*TMath::Sin(theta_la)*TMath::Cos(phi_la-phi_ra) + TMath::Cos(theta_ra)*TMath::Cos(theta_la);
						openingAngle_f = TMath::RadToDeg()*TMath::ACos(opa_foot);
						opa_f.push_back(openingAngle_f);

						// Selection of the right FOOT trck using CALIFA xyz Crystal information
						double xv, yv, zv, dist;
						double dist_min=99999.;
						double dist_minn=9999.;
						double xv_min=99999.;
						double yv_min=99999.;
						double zv_min=99999.;
						int p_min=-99999.;
						int q_min=-99999.;
						int min_d_f=0;
						int min_d_f_r=0;
						int min_d_f_t=0;
						int min_d_f_t1=0;
						int min_d_f_t2=0;

						double dist_min1=99999.;
						double dist_min2=99999.;
						Mult_left=Mult_left_y;
						Mult_right=Mult_right_y;
						ml.push_back(Mult_left);
						mr.push_back(Mult_right);
						TVector3 plm,plr,prm,prl;
						double Foot13_Zrn=-999.;
						double Foot7_Zln=-999.;
						double Foot10_Zln=-999.;
						double Foot6_Zrn=-999.;

						TVector3 pr,pl;
						double Foot13_Zn=-999.;
						double Foot7_Zn=-999.;
						double Foot10_Zn=-999.;
						double Foot6_Zn=-999.;

						TVector3 prt,plt;
						double Foot13_Znt=-999.;
						double Foot7_Znt=-999.;
						double Foot10_Znt=-999.;
						double Foot6_Znt=-999.;


						int dist_min_i=0;
						if(Mult_left>0 && Mult_right>0){

							for(Int_t p=0 ; p<Mult_right; p++){
								for(Int_t q=0 ; q<Mult_left; q++){
									Double_t temp_dist = 0.;
									Foot13_Znt = (BX_right.at(p)-fBDet[12])/(fADet[12]-AX_right.at(p));
									Foot7_Znt = (BX_left.at(q)-fBDet[6])/(fADet[6]-AX_left.at(q));
									Foot6_Znt = (BX_right.at(p)-fBDet[5])/(fADet[5]-AX_right.at(p));
									Foot10_Znt = (BX_left.at(q)-fBDet[9])/(fADet[9]-AX_left.at(q));
									TVector3 Posr_1t(AX_right.at(p)*Foot6_Znt+BX_right.at(p),AY_right.at(p)*Foot6_Znt+BY_right.at(p),Foot6_Znt);
									TVector3 Posl_1t(AX_left.at(q)*Foot10_Znt+BX_left.at(q),AY_left.at(q)*Foot10_Znt+BY_left.at(q),Foot10_Znt);
									TVector3 Posrt(AX_right.at(p)*Foot13_Znt+BX_right.at(p),AY_right.at(p)*Foot13_Znt+BY_right.at(p),Foot13_Znt);
									TVector3 Poslt(AX_left.at(q)*Foot7_Znt+BX_left.at(q),AY_left.at(q)*Foot7_Znt+BY_left.at(q),Foot7_Znt);
									//prt=Posrt-Vertext;
									//plt=Poslt-Vertext;
									prt=Posr_1t-Posrt;
									plt=Posl_1t-Poslt;       
									double theta_rtt=prt.Theta()*TMath::RadToDeg();
									double theta_ltt=plt.Theta()*TMath::RadToDeg();
									double phi_rtt=prt.Phi()*TMath::RadToDeg();
									double phi_ltt=plt.Phi()*TMath::RadToDeg();
									//if(abs(theta_-theta_ltt)<5 && abs(theta_r-theta_rtt)<5){}	

									TVector3 c_t_rn;
									TVector3 c_t_ln;
									TVector3 c_pos_ln(CALIFA_IDtoX[c_id_l-2433]*10,CALIFA_IDtoY[c_id_l-2433]*10,CALIFA_IDtoZ[c_id_l-2433]*10+490);
									TVector3 c_pos_rn(CALIFA_IDtoX[c_id_r-2433]*10,CALIFA_IDtoY[c_id_r-2433]*10,CALIFA_IDtoZ[c_id_r-2433]*10+490);
									c_t_rn=c_pos_rn-Posrt;
									c_t_ln=c_pos_ln-Poslt;

									double c_theta_ln=c_t_ln.Theta()*TMath::RadToDeg();
									double c_phi_ln=c_t_ln.Phi()*TMath::RadToDeg();
									double c_theta_rn=c_t_rn.Theta()*TMath::RadToDeg();
									double c_phi_rn=c_t_rn.Phi()*TMath::RadToDeg();

									double dphi_ln= -999.;
									double dphi_rn= -999.;
									double dtheta_ln= -999.;
									double dtheta_rn= -999.;


									if(c_theta_ln<90 && c_theta_rn<90 && abs(c_phi_rn)>90 && abs(c_phi_ln)<90){
										//double dphi_l= phi_lt-c_phi_r;
										//double dphi_r= phi_rt-c_phi_l;
										//double dtheta_l=theta_lt-c_theta_l;
										//double dtheta_r=theta_rt-c_theta_r;
										dphi_ln= c_phi_ln;
										dphi_rn= c_phi_rn;
										dtheta_ln=c_theta_ln;
										dtheta_rn=c_theta_rn;
									}

									MinDist(AX_right.at(p),BX_right.at(p),AY_right.at(p),BY_right.at(p),AX_left.at(q),BX_left.at(q),AY_left.at(q),BY_left.at(q),xv,yv,zv,dist);
									temp_dist = dist; //+ removed
									if(theta_l>0 && theta_r>0){

										if(abs(dtheta_ln-theta_ltt)<10 && abs(dtheta_rn-theta_rtt)<10 && dist_min>temp_dist && TMath::Abs(xv)<30. && TMath::Abs(yv)<20./* xv>-20  && xv<7*/){
											dist_min = temp_dist;
											//dist_minn = dist_min;
											p_min = p;
											q_min = q;
											min_d_f_t=1;
										}
									}
								}
							}

							if(dist_min<9999. && p_min>=0 && q_min>=0  && min_d_f_t==1){
								MinDist(AX_right.at(p_min),BX_right.at(p_min),AY_right.at(p_min),BY_right.at(p_min),AX_left.at(q_min),BX_left.at(q_min),AY_left.at(q_min),BY_left.at(q_min),xv_min,yv_min,zv_min,dist_min);

								double theta_prmt=999.;
								double theta_plmt=999.;
								double phi_prmt=999.;
								double phi_plmt=999.;
								if(dist_min<10){	
									if(AX_right.at(p_min)>0 && AX_left.at(q_min)<0){

										//cout<<"---------*BINGOOO*----------"<<endl;
										count_v++;	
										z_vm_t.push_back(zv_min);
										y_vm_t.push_back(yv_min);
										x_vm_t.push_back(xv_min);
										dm_t.push_back(dist_min);
										TVector3 Vertex(xv_min,yv_min,zv_min);

										Foot13_Zn = (BX_right.at(p_min)-fBDet[12])/(fADet[12]-AX_right.at(p_min));
										Foot7_Zn = (BX_left.at(q_min)-fBDet[6])/(fADet[6]-AX_left.at(q_min));
										Foot6_Zn = (BX_right.at(p_min)-fBDet[5])/(fADet[5]-AX_right.at(p_min));
										Foot10_Zn = (BX_left.at(q_min)-fBDet[9])/(fADet[9]-AX_left.at(q_min));
										TVector3 Posr_1(AX_right.at(p_min)*Foot6_Zn+BX_right.at(p_min),AY_right.at(p_min)*Foot6_Zn+BY_right.at(p_min),Foot6_Zn);
										TVector3 Posl_1(AX_left.at(q_min)*Foot10_Zn+BX_left.at(q_min),AY_left.at(q_min)*Foot10_Zn+BY_left.at(q_min),Foot10_Zn);
										TVector3 Posr(AX_right.at(p_min)*Foot13_Zn+BX_right.at(p_min),AY_right.at(p_min)*Foot13_Zn+BY_right.at(p_min),Foot13_Zn);
										TVector3 Posl(AX_left.at(q_min)*Foot7_Zn+BX_left.at(q_min),AY_left.at(q_min)*Foot7_Zn+BY_left.at(q_min),Foot7_Zn);
										pr=Posr-Vertex;
										pl=Posl-Vertex;
										//	pr=Posr_1-Posr;
										//	pl=Posl_1-Posl;       
										double theta_rt=pr.Theta();
										double theta_lt=pl.Theta();
										double phi_rt=pr.Phi();
										double phi_lt=pl.Phi();

										theta_p_r_m_t.push_back(theta_lt*TMath::RadToDeg());
										phi_p_r_m_t.push_back(phi_lt*TMath::RadToDeg());

										theta_p_l_m_t.push_back(theta_rt*TMath::RadToDeg());
										phi_p_l_m_t.push_back(phi_rt*TMath::RadToDeg());
										double openingAngle_ft=999.;
										opa_foot_t = TMath::Sin(theta_lt)*TMath::Sin(theta_rt)*TMath::Cos(phi_rt-phi_lt) + TMath::Cos(theta_lt)*TMath::Cos(theta_rt);
										openingAngle_ft = TMath::RadToDeg()*TMath::ACos(opa_foot_t);
										opa_f_t.push_back(openingAngle_ft);
										opa_t.push_back(openingAngle_ft);
										dist_min_i=1;

										TVector3 cen(0,0,441);
										TVector3 c_t_r;
										TVector3 c_t_l;
										TVector3 c_pos_l(CALIFA_IDtoX[c_id_l-2433]*10,CALIFA_IDtoY[c_id_l-2433]*10,CALIFA_IDtoZ[c_id_l-2433]*10+498);
										TVector3 c_pos_r(CALIFA_IDtoX[c_id_r-2433]*10,CALIFA_IDtoY[c_id_r-2433]*10,CALIFA_IDtoZ[c_id_r-2433]*10+498);
										c_t_r=c_pos_r-Vertex;
										c_t_l=c_pos_l-Vertex;

										double mp = 0.938; // GeV
										double m_16C = 14.9176; // GeV
										double m_12C = 11.178; // GeV
										double c_theta_l=c_t_l.Theta();
										double c_phi_l=c_t_l.Phi();
										double c_theta_r=c_t_r.Theta();
										double c_phi_r=c_t_r.Phi();
										double rec_leftf = recon(El/1000.,theta_rt)/1000.;
										double rec_rightf = recon(Er/1000.,theta_rt)/1000.;
										double p3_magf = TMath::Sqrt(pow(rec_leftf+mp,2)-pow(mp,2));
										double p4_magf = TMath::Sqrt(pow(rec_rightf+mp,2)-pow(mp,2));
										TVector3 P3f,P4f;
										P3f.SetMagThetaPhi(p3_magf,theta_rt,phi_rt);
										P4f.SetMagThetaPhi(p4_magf,theta_lt,phi_lt);
										p3f.SetVectM(P3f,mp);
										p4f.SetVectM(P4f,mp);
										ptg.SetXYZM(0.,0.,0.,mp);
										pmissf = p3f+p4f-ptg;
										pbeam.SetXYZM(0.,0.,TMath::Sqrt(pow(1.244029*16+m_16C,2)-pow(m_16C,2)),m_16C);
										//pbeam.SetXYZM(0.,0.,TMath::Sqrt(pow(1.25*12+m_12C,2)-pow(m_12C,2)),m_12C);
										pmissf.Boost(-pbeam.BoostVector());	
										Emissf=mp-pmissf.E();
										mmissf=pmissf.M2();
										Emf.push_back(Emissf);
										Mmf.push_back(mmissf);
										Pmtf.push_back(pmissf.Mag());
										Pmxf.push_back(pmissf.X());
										Pmyf.push_back(pmissf.Y());
										Pmzf.push_back(pmissf.Z());

										if(c_theta_l*TMath::RadToDeg()<90 && c_theta_r*TMath::RadToDeg()<90 && abs(c_phi_r*TMath::RadToDeg())>90 && abs(c_phi_l*TMath::RadToDeg())<90){

											double dphi_l= c_phi_l;
											double dphi_r= c_phi_r;
											double dtheta_l=c_theta_l;
											double dtheta_r=c_theta_r;

											dpl.push_back(dphi_l*TMath::RadToDeg());
											dpr.push_back(dphi_r*TMath::RadToDeg());
											dtl.push_back(dtheta_l*TMath::RadToDeg());
											dtr.push_back(dtheta_r*TMath::RadToDeg());

											double rec_left = recon(El/1000.,dtheta_l)/1000.;
											double rec_right = recon(Er/1000.,dtheta_r)/1000.;

											double p3_mag = TMath::Sqrt(pow(rec_left+mp,2)-pow(mp,2));
											double p4_mag = TMath::Sqrt(pow(rec_right+mp,2)-pow(mp,2));
											TVector3 P3,P4;
											P3.SetMagThetaPhi(p3_mag,dtheta_l,dphi_l);
											P4.SetMagThetaPhi(p4_mag,dtheta_r,dphi_r);
											p3.SetVectM(P3,mp);
											p4.SetVectM(P4,mp);
											pmiss = p3+p4-ptg;
											pmiss.Boost(-pbeam.BoostVector());	
											Emiss=mp-pmiss.E();
											mmiss=pmiss.M2();
											Em.push_back(Emiss);
											Mm.push_back(mmiss);
											Pmt.push_back(pmiss.Mag());
											Pmx.push_back(pmiss.X());
											Pmy.push_back(pmiss.Y());
											Pmz.push_back(pmiss.Z());

										}

									}
								}
							}


						}

					}
				}//TPAT
			}//mwpcaX
			}//FRS	
			vertex->Fill();
		}
		vertex->Write();
		q->Close();
		Double_t rtime = timer.RealTime();
		Double_t ctime = timer.CpuTime();
		Float_t cpuUsage = ctime / rtime;
		//cout << "CPU used: " << cpuUsage << endl;
		cout << endl;
		cout << "Real time " << rtime/60 << " min, CPU time " << ctime/60 <<"min" << endl << endl;
		cout << "Macro finished successfully." << endl;
		cout << "Number of vertex found is--> " << count_v << endl;
		cout << "Percentage of vertex found is--> " << (count_v/nevents)*100 << endl;
		//cout << "s522 16C without correction => p2p true + TCutG" << endl;

	}




