#ifndef tnp_weight_h
#define tnp_weight_h

#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++
//
// - MuID: (tnp_weight_muid_pbpb)   Preliminary
//   * idx = 0: nominal
//   * idx = -1: syst variation,  +1 sigma
//   * idx = -2: syst variation,  -1 sigma
//   * idx = +1: stat variation,  +1 sigma
//   * idx = +2: stat variation,  -1 sigma
//
// - Trigger: (tnp_weight_trg_pbpb)  Preliminary
//   * idx = 0:  nominal
//   * idx = -1: TnP syst variation,  +1 sigma
//   * idx = -2: TnP syst variation,  -1 sigma
//   * idx = +1: stat variation,  +1 sigma
//   * idx = +2: stat variation,  -1 sigma
//   * idx = -3: run variation,  +1 sigma
//   * idx = -4: run variation,  -1 sigma

// For all:
//   * idx = +200: tnp efficiency for data
//   * idx = +300: tnp efficiency for MC

// THE INDIVIDUAL SFs
// ++++++++++++++++++
double tnp_weight_muid_pbpb(double eta, int idx=0);
double tnp_weight_trig_pbpb(double pt, double eta, int idx=0);



///////////////////////////////////////////////////
//                 M u I D    P b P b                //
///////////////////////////////////////////////////
double tnp_weight_muid_pbpb(double eta, int idx)
{
   double x = eta;
   double syst = 0.6e-2;  //preliminary

   double num=1,den=1;
   

   // MC
   if (x > -2.4&&x <= -2.1) den = 0.994717;
   if (x > -2.1&&x <= -1.6) den = 0.993653;
   if (x > -1.6&&x <= -1.2) den = 0.982687;
   if (x > -1.2&&x <= -0.9) den = 0.962992;
   if (x > -0.9&&x <= -0.6) den = 0.970103;
   if (x > -0.6&&x <= -0.3) den = 0.982073;
   if (x > -0.3&&x <= 0) den = 0.968219;
   if (x > 0 && x <= 0.3) den = 0.961703;
   if (x > 0.3&&x <= 0.6) den = 0.978947;
   if (x > 0.6&&x <= 0.9) den = 0.970324;
   if (x > 0.9&&x <= 1.2) den = 0.956433;
   if (x > 1.2&&x <= 1.6) den = 0.982577;
   if (x > 1.6&&x <= 2.1) den = 0.995146;
   if (x > 2.1&&x <= 2.4) den = 0.994106;


   // data
   if (idx <= 0 || idx > 10) { // nominal
	   if (x > -2.4&&x <= -2.1) num = 0.983734;
	   if (x > -2.1&&x <= -1.6) num = 0.993892;
	   if (x > -1.6&&x <= -1.2) num = 0.979222;
	   if (x > -1.2&&x <= -0.9) num = 0.955529;
	   if (x > -0.9&&x <= -0.6) num = 0.966578;
	   if (x > -0.6&&x <= -0.3) num = 0.983958;
	   if (x > -0.3&&x <= 0) num = 0.959491;
	   if (x > 0 && x <= 0.3) num = 0.959744;
	   if (x > 0.3&&x <= 0.6) num = 0.976402;
	   if (x > 0.6&&x <= 0.9) num = 0.96839;
	   if (x > 0.9&&x <= 1.2) num = 0.961177;
	   if (x > 1.2&&x <= 1.6) num = 0.98078;
	   if (x > 1.6&&x <= 2.1) num = 0.991852;
	   if (x > 2.1&&x <= 2.4) num = 0.993082;
   }
   else if (idx == 1) { // stat up
	   if (x > -2.4&&x <= -2.1) num = 0.98710825;
	   if (x > -2.1&&x <= -1.6) num = 0.99575701;
	   if (x > -1.6&&x <= -1.2) num = 0.98250043;
	   if (x > -1.2&&x <= -0.9) num = 0.96029298;
	   if (x > -0.9&&x <= -0.6) num = 0.97057735;
	   if (x > -0.6&&x <= -0.3) num = 0.98705592;
	   if (x > -0.3&&x <= 0) num = 0.96370511;
	   if (x > 0 && x <= 0.3) num = 0.96396652;
	   if (x > 0.3&&x <= 0.6) num = 0.97980996;
	   if (x > 0.6&&x <= 0.9) num = 0.97234624;
	   if (x > 0.9&&x <= 1.2) num = 0.9659164;
	   if (x > 1.2&&x <= 1.6) num = 0.98388187;
	   if (x > 1.6&&x <= 2.1) num = 0.99371692;
	   if (x > 2.1&&x <= 2.4) num = 0.9961665;
   }
   else if (idx == 2) { // stat down
	   if (x > -2.4&&x <= -2.1) num = 0.98075691;
	   if (x > -2.1&&x <= -1.6) num = 0.99224693;
	   if (x > -1.6&&x <= -1.2) num = 0.97618383;
	   if (x > -1.2&&x <= -0.9) num = 0.95101985;
	   if (x > -0.9&&x <= -0.6) num = 0.96279449;
	   if (x > -0.6&&x <= -0.3) num = 0.98108946;
	   if (x > -0.3&&x <= 0) num = 0.95550071;
	   if (x > 0 && x <= 0.3) num = 0.95575316;
	   if (x > 0.3&&x <= 0.6) num = 0.97322907;
	   if (x > 0.6&&x <= 0.9) num = 0.96468663;
	   if (x > 0.9&&x <= 1.2) num = 0.9567025;
	   if (x > 1.2&&x <= 1.6) num = 0.9779218;
	   if (x > 1.6&&x <= 2.1) num = 0.99018962;
	   if (x > 2.1&&x <= 2.4) num = 0.9905424;
   }

   if (idx == 200) den = 1.;
   if (idx == 300) num = den * den;


   double syst_factor = 1.;
   if (idx == -1) syst_factor = 1 + syst;
   if (idx == -2) syst_factor = 1 - syst;
   return (num / den)*syst_factor;
}


///////////////////////////////////////////////////
//               T R G      P b P b              //
///////////////////////////////////////////////////
double tnp_weight_trig_pbpb(double pt, double eta, int idx)
{
  double num=1.0, den=1.0;

  // MC
  if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
    if (pt > 15 && pt <= 20) den = 0.940563;
    else if (pt > 20 && pt <= 30) den = 0.953007;
    else if (pt > 30 && pt <= 50) den = 0.964625;
    else if (pt > 50 && pt <= 80) den = 0.966917;
    else if (pt > 80 && pt <= 9999) den = 0.958959;
  }
  else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
    if (pt > 15 && pt <= 20) den = 0.9281;
    else if (pt > 20 && pt <= 30) den = 0.940129;
    else if (pt > 30 && pt <= 50) den = 0.953978;
    else if (pt > 50 && pt <= 80) den = 0.958659;
    else if (pt > 80 && pt <= 9999) den = 0.961068;
  }
  else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
    if (pt > 15 && pt <= 20) den = 0.890849;
    else if (pt > 20 && pt <= 30) den = 0.914955;
    else if (pt > 30 && pt <= 50) den = 0.934417;
    else if (pt > 50 && pt <= 80) den = 0.943679;
    else if (pt > 80 && pt <= 9999) den = 0.947818;
  }

  // data
  if (idx <= 0 || idx > 10) { // nominal
    if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
      if (pt > 15 && pt <= 20) num = 0.875759;
      else if (pt > 20 && pt <= 30) num = 0.92399;
      else if (pt > 30 && pt <= 50) num = 0.948799;
      else if (pt > 50 && pt <= 80) num = 0.949345;
      else if (pt > 80 && pt <= 9999) num = 0.940723;
    }
    else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
      if (pt > 15 && pt <= 20) num = 0.893974;
      else if (pt > 20 && pt <= 30) num = 0.906569;
      else if (pt > 30 && pt <= 50) num = 0.926756;
      else if (pt > 50 && pt <= 80) num = 0.940292;
      else if (pt > 80 && pt <= 9999) num = 0.928692;
    }
    else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
      if (pt > 15 && pt <= 20) num = 0.828836;
      else if (pt > 20 && pt <= 30) num = 0.871728;
      else if (pt > 30 && pt <= 50) num = 0.898965;
      else if (pt > 50 && pt <= 80) num = 0.91405;
      else if (pt > 80 && pt <= 9999) num = 0.889446;
    }
  }
  else if (idx == 1) { // stat up
    if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
      if (pt > 15 && pt <= 20) num = 0.890676;
      else if (pt > 20 && pt <= 30) num = 0.928769;
      else if (pt > 30 && pt <= 50) num = 0.952176;
      else if (pt > 50 && pt <= 80) num = 0.954073;
      else if (pt > 80 && pt <= 9999) num = 0.953832;
    }
    else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
      if (pt > 15 && pt <= 20) num = 0.908596;
      else if (pt > 20 && pt <= 30) num = 0.913711;
      else if (pt > 30 && pt <= 50) num = 0.929705;
      else if (pt > 50 && pt <= 80) num = 0.94689;
      else if (pt > 80 && pt <= 9999) num = 0.948762;
    }
    else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
      if (pt > 15 && pt <= 20) num = 0.85611;
      else if (pt > 20 && pt <= 30) num = 0.886324;
      else if (pt > 30 && pt <= 50) num = 0.905823;
      else if (pt > 50 && pt <= 80) num = 0.929945;
      else if (pt > 80 && pt <= 9999) num = 0.940054;
    }
  }
  else if (idx == 2) { // stat down
    if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
      if (pt > 15 && pt <= 20) num = 0.859607;
      else if (pt > 20 && pt <= 30) num = 0.919003;
      else if (pt > 30 && pt <= 50) num = 0.945421;
      else if (pt > 50 && pt <= 80) num = 0.944336;
      else if (pt > 80 && pt <= 9999) num = 0.925502;
    }
    else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
      if (pt > 15 && pt <= 20) num = 0.878284;
      else if (pt > 20 && pt <= 30) num = 0.899137;
      else if (pt > 30 && pt <= 50) num = 0.923753;
      else if (pt > 50 && pt <= 80) num = 0.933253;
      else if (pt > 80 && pt <= 9999) num = 0.904548;
    }
    else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
      if (pt > 15 && pt <= 20) num = 0.798669;
      else if (pt > 20 && pt <= 30) num = 0.856249;
      else if (pt > 30 && pt <= 50) num = 0.891814;
      else if (pt > 50 && pt <= 80) num = 0.896099;
      else if (pt > 80 && pt <= 9999) num = 0.819963;
    }
  }
  if (idx == -1) { // TnP fit syst up
    if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
      if (pt > 15 && pt <= 20) num = 0.901332;
      else if (pt > 20 && pt <= 30) num = 0.925191;
      else if (pt > 30 && pt <= 50) num = 0.949999;
      else if (pt > 50 && pt <= 80) num = 0.950461;
      else if (pt > 80 && pt <= 9999) num = 0.9415;
    }
    else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
      if (pt > 15 && pt <= 20) num = 0.90741;
      else if (pt > 20 && pt <= 30) num = 0.909254;
      else if (pt > 30 && pt <= 50) num = 0.9299;
      else if (pt > 50 && pt <= 80) num = 0.944354;
      else if (pt > 80 && pt <= 9999) num = 0.933213;
    }
    else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
      if (pt > 15 && pt <= 20) num = 0.848026;
      else if (pt > 20 && pt <= 30) num = 0.87791;
      else if (pt > 30 && pt <= 50) num = 0.905147;
      else if (pt > 50 && pt <= 80) num = 0.920232;
      else if (pt > 80 && pt <= 9999) num = 0.895628;
    }
  }
  else if (idx == -2) { // TnP fit syst down
    if (fabs(eta) > 0 && fabs(eta) <= 1.2) { 
      if (pt > 15 && pt <= 20) num = 0.850186;
      else if (pt > 20 && pt <= 30) num = 0.92279;
      else if (pt > 30 && pt <= 50) num = 0.947598;
      else if (pt > 50 && pt <= 80) num = 0.948229;
      else if (pt > 80 && pt <= 9999) num = 0.939947;
    }
    else if (fabs(eta) > 1.2 && fabs(eta) <= 2.1) { 
      if (pt > 15 && pt <= 20) num = 0.880538;
      else if (pt > 20 && pt <= 30) num = 0.903885;
      else if (pt > 30 && pt <= 50) num = 0.923613;
      else if (pt > 50 && pt <= 80) num = 0.93623;
      else if (pt > 80 && pt <= 9999) num = 0.92417;
    }
    else if (fabs(eta) > 2.1 && fabs(eta) <= 2.4) { 
      if (pt > 15 && pt <= 20) num = 0.809645;
      else if (pt > 20 && pt <= 30) num = 0.865546;
      else if (pt > 30 && pt <= 50) num = 0.892783;
      else if (pt > 50 && pt <= 80) num = 0.907868;
      else if (pt > 80 && pt <= 9999) num = 0.883264;
    }
  }
  else if (idx == -3) { num *= (1.0+0.006); } // run-dependence syst up
  else if (idx == -4) { num *= (1.0-0.006); } // run-dependence syst down

  if (idx == 200) den = 1.0;
  if (idx == 300) num = den * den;

  return (num/den);
}

#endif //#ifndef tnp_weight_h
