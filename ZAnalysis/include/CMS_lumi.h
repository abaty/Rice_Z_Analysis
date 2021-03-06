#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
TString extraText   = "Supplementary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
//float lumiTextSize     = 0.6;
float lumiTextSize     = 0.4;
float lumiTextOffset   = 0.2;
//float cmsTextSize      = 0.75;
float cmsTextSize      = 0.4;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "1.7 nb^{-1} (5.02 TeV PbPb)";
TString lumi_sqrtS090 = "1.7 nb^{-1} (5.02 TeV PbPb 0-100%)";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10, float fontMultiplier = 1.0, bool doCMS = true, bool doLumi = true, bool doExtraText = false, bool do090Lumi = false, bool increaseExtraOffset = false );

