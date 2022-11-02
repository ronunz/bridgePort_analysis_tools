#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TPolyMarker.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TExec.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TTree.h"

#include <vector>
#include "time.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Define Global Variables here
//const char *filename="GeGI_20200117-1_Simulation8.root"; // Filename of simulated GeGI data
const char *filename="GeGI_20200117-1_Simulation9.root"; // Updated - MOST RECENT SIMULATION - Filename of simulated GeGI data
//const char *filename2="RawEventFileCs137.root"; // Filename of measured GeGI data
const char *filename2="20200604-1_Eu154155_10min.root";
const char *SRMfilename="SRM.txt";

const char *Strip_X_Name="hEnergyDepositGeDSSDx_%d"; // Name of the x-strips
const char *Strip_Y_Name="hEnergyDepositGeDSSDy_%d"; // Name of the y-strips
const char *GeGI_Strip_X_Name="GeGI Channel %d"; // Name of the x-strips from GeGI data
const char *GeGI_Strip_Y_Name="GeGI Channel %d"; // Name of the y-strips from GeGI data

TFile *fin = new TFile(filename); // Load Simulation data
TFile *fin2 = new TFile(filename2); // Load Measured data (GeGI)

// Source information   - Efficiency Calculations
double Photon_Emission=1778278108; // Should match number of GPS from Macro
                                   // Number of photons emitted from 94.1 uCi Cs-137 over 10 minutes
double Cs_Peak_Location=661.66;
//double Peak_Location = Cs_Peak_Location;
double Peak_Location = 86;
double Peak_Width = 10;

// Global Histogram Variables
TH1D *h1; // Histograms for simulated x-channel
TH1D *h2; // Histograms for simulated y-channel
TH1D *GeGI_Histograms[32]; //Histograms for measured x AND y channels. 0-15 y, 16-31 x.
TH1D *GeGI_Total; // Histogram for sum of all strips
TH1D *GeGI_Total_X; // Histogram for sum all x strips
TH1D *GeGI_Total_Y; // HIstogram for sum all y strips
TH2F *GeGI_Image; // 2D Histogram to create detector images

// Declaration of Functions
TF1 *f1[16]; // Fitting for simulated x strips
TF1 *f2[16]; // Fitting for simulated y strips
TF1 *f_Total; // Fitting for sum of channels histogram
TF1 *GeGI_f1[16]; // Fitting for mesaured x strips
TF1 *GeGI_f2[16]; // Fitting for measured y strips
TF1 *GeGI_f_Total_X; // Fitting for sum of x channel histograms
TF1 *GeGI_f_Total_Y; // Fitting for sum of y channel histograms

// Canvas Declartions
TCanvas *c1;

// Declaration of Fit Paramater for each of the 16 X and Y SIMULATED strips
double SigmaX[16];
double SigmaX_Error[16];
double MeanX[16];
double MeanX_Error[16];
double AmplitudeX[16];
double AmplitudeX_Error[16];
double AreaX[16];
double AreaX_Error[16];
double ResolutionX[16];
double ResolutionX_Error[16];
double EfficiencyX[16];
double EfficiencyX_Error[16];
double SigmaY[16];
double SigmaY_Error[16];
double MeanY[16];
double MeanY_Error[16];
double AmplitudeY[16];
double AmplitudeY_Error[16];
double AreaY[16];
double AreaY_Error[16];
double ResolutionY[16];
double ResolutionY_Error[16];
double EfficiencyY[16];
double EfficiencyY_Error[16];
double Strip_Number[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // Need an array with the listed strip number for the TGraphErrors

// Declaration of Fit Paramters for each of the 16 X and Y MEASURED strips
double GeGI_SigmaX[16];
double GeGI_SigmaX_Error[16];
double GeGI_MeanX[16];
double GeGI_MeanX_Error[16];
double GeGI_AmplitudeX[16];
double GeGI_AmplitudeX_Error[16];
double GeGI_AreaX[16];
double GeGI_AreaX_Error[16];
double GeGI_ResolutionX[16];
double GeGI_ResolutionX_Error[16];
double GeGI_EfficiencyX[16];
double GeGI_EfficiencyX_Error[16];
double GeGI_SigmaY[16];
double GeGI_SigmaY_Error[16];
double GeGI_MeanY[16];
double GeGI_MeanY_Error[16];
double GeGI_AmplitudeY[16];
double GeGI_AmplitudeY_Error[16];
double GeGI_AreaY[16];
double GeGI_AreaY_Error[16];
double GeGI_ResolutionY[16];
double GeGI_ResolutionY_Error[16];
double GeGI_EfficiencyY[16];
double GeGI_EfficiencyY_Error[16];
// Mannual Calculations
double GeGI_AreaX_Manual[16];
double GeGI_EfficiencyX_Manual[16];
double GeGI_AreaY_Manual[16];
double GeGI_EfficiencyY_Manual[16];


// Declaration of Variables for the TOTAL of SIMULATED strips
double Total_Amplitude;
double Total_Amplitude_Error;
double Total_Mean;
double Total_Mean_Error;
double Total_Sigma;
double Total_Sigma_Error;
double Total_Area;
double Total_Area_Error;
double Total_Resolution;
double Total_Resolution_Error;
double Total_Efficiency;
double Total_Efficiency_Error;
double Total_Area_From_Strips;
double Total_Efficiency_From_Strips;


// Declaration of variables for the TOAL of MEASURED strips (x and y)
double GeGI_Total_Amplitude_X;
double GeGI_Total_Amplitude_X_Error;
double GeGI_Total_Mean_X;
double GeGI_Total_Mean_X_Error;
double GeGI_Total_Sigma_X;
double GeGI_Total_Sigma_X_Error;
double GeGI_Total_Area_X;
double GeGI_Total_Area_X_Error;
double GeGI_Total_Resolution_X;
double GeGI_Total_Resolution_X_Error;
double GeGI_Total_Efficiency_X;
double GeGI_Total_Efficiency_X_Error;

double GeGI_Total_Amplitude_Y;
double GeGI_Total_Amplitude_Y_Error;
double GeGI_Total_Mean_Y;
double GeGI_Total_Mean_Y_Error;
double GeGI_Total_Sigma_Y;
double GeGI_Total_Sigma_Y_Error;
double GeGI_Total_Area_Y;
double GeGI_Total_Area_Y_Error;
double GeGI_Total_Resolution_Y;
double GeGI_Total_Resolution_Y_Error;
double GeGI_Total_Efficiency_Y;
double GeGI_Total_Efficiency_Y_Error;

// Declaration of varible for SRM Input
int SRM_Num_of_Isotopes;
double SRM_Source_HL_Activity[5][2]; // Source half-life and acitivy informaiton
double SRM_Emission[5][20][3];  // 5 potential isotopes, 20 possible emission energies, value 1 - energy, value 2 probability, value 3 error  
//double SRM_Emission_Probability[5][20];  //  Emission probabilities
//double SRM_Emission_Probability_Error[5][20]; // Emission probabiity error
struct tm SRM_Date; // Reference date, year, month, day, hour
struct tm Run_Date;
double Time_Difference; // Time between SRM and Run data to adjust for source decay


// Constants
double Channel_Energy_Width = 0.5; // Determined when constructing the ICET_simAnalysisExample histograms (6000,0,3000) 0.5 keV/channel
int Number_Of_Channels=16;
double Pi=3.14152; // Value of Pi

// Efficiency Calculations for simulated data
void Fitting_Efficiency_Calcs_X();
void Fitting_Efficiency_Calcs_Y();

// Simulated Efficiency calculations for strip totals or all interactions within detector volme
void Fitting_Efficiency_Calcs_Total();
void Channel_Sum_Area_Efficiency();

// GeGI Efficiency fitting calculations
void Fitting_GeGI_Efficiency_Calcs_X();
void Fitting_GeGI_Efficiency_Calcs_Y();

// GeGI Efficiency calculations for strip totals or all interactings within detector volume for both x and y sides
void Fitting_GeGI_Efficiency_Calcs_Total_X();
void Fitting_GeGI_Efficiency_Calcs_Total_Y();

// Manual integration of peak areas for GeGI data (due to poor gaussian fittings to data)
void Fitting_GeGI_Efficiency_Manual_Calcs_X();
void Fitting_GeGi_Efficiency_Manual_Calcs_Y();

// Function for constructing GeGI histograms from RootConvert application
void Make_GeGI_Histograms();

// Function for calling the specific functions for plotting
void Plotting_Histograms();

// Ploting the fficiencies for simulated data
void Efficiency_Plot_X();
void Efficiency_Plot_Y();
void Efficiency_Plot_XY();

// Plot the efficiencies for GeGI data
void GeGI_Efficiency_Plot_X();
void GeGI_Efficiency_Plot_Y();
void GeGI_Efficiency_Plot_XY();

// Comparing x and y strip efficiencies on same canvas (x on one and y on another canvas)
void Comparison_Efficiency_Plot_X();
void Comparison_Efficiency_Plot_Y();

// Functions for plotting all 16 channels on single canvas
void Simulated_Multiplot_Strips_X();
void Simulated_Miltiplot_Strips_Y();

void GeGI_Multiplot_Strips_X();
void GeGI_Multiplot_Strips_Y();

// Function for inport of SRM File

void SRM_Input();


//-------------------------------------------------------------------------------------------------
//
//                                        MAIN FUNCTION
//
//------------------------------------------------------------------------------------------------
void GeGI_Analyzer(){
  
  // Open file
  cout << endl << endl << endl;
  cout << "  Simulation File: " << filename << endl;
  cout << "  Measurement File: " << filename2 << endl;
  cout << endl << endl;

  // Open SRM File
  SRM_Input();
  
  //------------------------
  // SIMULATION CALCULATIONS
  //------------------------

  //-----------------
  // Simulation fittings
  //-----------------
  //Fitting_Efficiency_Calcs_X(); // Fits for all x channels
  //Fitting_Efficiency_Calcs_Y(); // Fits for all y channels
  //Fitting_Efficiency_Calcs_Total(); // Fits for sum (hEnergyDepositcal)
  
  //Channel_Sum_Area_Efficiency(); // Calculating and comparing data from Total and Sum of Strips

  //------------------
  // GEGI CALCULATIONS
  //------------------
  
  //--------------------------------------------------
  // Make GeGI Histograms - REQUIRED FOR GeGI ANALYSIS
  //--------------------------------------------------
  Make_GeGI_Histograms();
 

  //---------------
  // GeGI Fittings
  //--------------
  Fitting_GeGI_Efficiency_Calcs_Y();  
  Fitting_GeGI_Efficiency_Calcs_X();
  
  //Fitting_GeGI_Efficiency_Calcs_Total_X();
  //Fitting_GeGI_Efficiency_Calcs_Total_Y();

  // Does not use guassian fittings, but utilizes summed areas within the spectrum following methods 
  // in Gilmore Gamma Spec. Book
  //Fitting_GeGI_Efficiency_Manual_Calcs_X();
  //Fitting_GeGi_Efficiency_Manual_Calcs_Y();

  //---------
  // PLOTTING
  //---------
  Plotting_Histograms();
  

}










//-------------------------------------------------------------------------------------------------
//
//                               SIMULATED EFFICIENCY CALCULATIONS
//
//-------------------------------------------------------------------------------------------------

void Fitting_Efficiency_Calcs_X(){

  for(int Strip_X_Number=0; Strip_X_Number<Number_Of_Channels; Strip_X_Number++){

    // Define the histogram and get data from the already made histogram
    h1 = (TH1D*)fin->Get(Form(Strip_X_Name,Strip_X_Number));      

    // Setup up new function
    f1[Strip_X_Number] = new TF1(Form(Strip_X_Name,Strip_X_Number),"gaus(0)",Peak_Location-Peak_Width,Peak_Location+Peak_Width);
 
    // // Set initial paramters to assist with initial conditions
    f1[Strip_X_Number]->SetParameter(0,1);
    f1[Strip_X_Number]->SetParameter(1,Peak_Location);
    f1[Strip_X_Number]->SetParameter(2,5);

    // WE NEED TO SET RANGE AND ADJUST FOR THE PHOTOPEAKS (SEE COMMENTS BELOW)

    // // Fit the histogram parameters
    h1->Fit(Form(Strip_X_Name,Strip_X_Number),"RMNQ"); // look at fit options
  
    // // Get gaussian amplitude
    AmplitudeX[Strip_X_Number] = f1[Strip_X_Number]->GetParameter(0);
    AmplitudeX_Error[Strip_X_Number] = f1[Strip_X_Number]->GetParError(0);  
  
    // // Get the gaussian offset
    MeanX[Strip_X_Number] = f1[Strip_X_Number]->GetParameter(1);
    MeanX_Error[Strip_X_Number] = f1[Strip_X_Number]->GetParError(1);
  
    // // Get the gaussian standard deviation
    SigmaX[Strip_X_Number] = f1[Strip_X_Number]->GetParameter(2);
    SigmaX_Error[Strip_X_Number] = f1[Strip_X_Number]->GetParError(2);

    // // Calculate the area
    AreaX[Strip_X_Number] = sqrt(2*Pi)*AmplitudeX[Strip_X_Number]*SigmaX[Strip_X_Number]/Channel_Energy_Width;
    AreaX_Error[Strip_X_Number] = sqrt(2*Pi)*sqrt(pow((SigmaX[Strip_X_Number]*AmplitudeX_Error[Strip_X_Number]),2)+
						  pow((AmplitudeX[Strip_X_Number]*SigmaX_Error[Strip_X_Number]),2))/Channel_Energy_Width;

    // //Calculate the resolution
    ResolutionX[Strip_X_Number] = 100.0*2.3548*SigmaX[Strip_X_Number]/MeanX[Strip_X_Number];
    ResolutionX_Error[Strip_X_Number] = 2.3545*sqrt(pow(SigmaX_Error[Strip_X_Number]/MeanX[Strip_X_Number],2)+
						    pow(MeanX_Error[Strip_X_Number]/pow(MeanX[Strip_X_Number],2),2));

    // Calculate the Absolute Efficiency
    EfficiencyX[Strip_X_Number] = 100*AreaX[Strip_X_Number]/Photon_Emission;
    EfficiencyX_Error[Strip_X_Number] = 100*AreaX_Error[Strip_X_Number]/Photon_Emission;
    
    //cout << EfficiencyX[Strip_X_Number] << endl;
    
  }

  return;
}


void Fitting_Efficiency_Calcs_Y(){

  int Strip_Y_Number=0; 

  for(int Strip_Y_Number=0; Strip_Y_Number<Number_Of_Channels; Strip_Y_Number++){

    // Define the histogram and get data from the already made histogram
    h2 = (TH1D*)fin->Get(Form(Strip_Y_Name,Strip_Y_Number));    

    // Setup up new function
    f2[Strip_Y_Number] = new TF1(Form(Strip_Y_Name,Strip_Y_Number),"gaus(0)",Peak_Location-Peak_Width,Peak_Location+Peak_Width);
 
    // Set initial paramters to assist with initial conditions
    f2[Strip_Y_Number]->SetParameter(0,1);
    f2[Strip_Y_Number]->SetParameter(1,Peak_Location);
    f2[Strip_Y_Number]->SetParameter(2,5);

    // Fit the histogram parameters
    h2->Fit(Form(Strip_Y_Name,Strip_Y_Number),"RMNQ");
  
    // Get gaussian amplitude
    AmplitudeY[Strip_Y_Number] = f2[Strip_Y_Number]->GetParameter(0);
    AmplitudeY_Error[Strip_Y_Number] = f2[Strip_Y_Number]->GetParError(0);
    
    // Get the gaussian offset
    MeanY[Strip_Y_Number] = f2[Strip_Y_Number]->GetParameter(1);
    MeanY_Error[Strip_Y_Number] = f2[Strip_Y_Number]->GetParError(1);
  
    // Get the gaussian standard deviation
    SigmaY[Strip_Y_Number] = f2[Strip_Y_Number]->GetParameter(2);
    SigmaY_Error[Strip_Y_Number] = f2[Strip_Y_Number]->GetParError(2);

    // // Calculate the area
    AreaY[Strip_Y_Number] = sqrt(2*Pi)*AmplitudeY[Strip_Y_Number]*SigmaY[Strip_Y_Number]/Channel_Energy_Width;
    AreaY_Error[Strip_Y_Number] = sqrt(2*Pi)*sqrt(pow((SigmaY[Strip_Y_Number]*AmplitudeY_Error[Strip_Y_Number]),2)+
    						  pow((AmplitudeY[Strip_Y_Number]*SigmaY_Error[Strip_Y_Number]),2))/Channel_Energy_Width;
    
    // //Calculate the resolution
    ResolutionY[Strip_Y_Number] = 100.0*2.3548*SigmaY[Strip_Y_Number]/MeanY[Strip_Y_Number];
    ResolutionY_Error[Strip_Y_Number] = 2.3545*sqrt(pow(SigmaY_Error[Strip_Y_Number]/MeanY[Strip_Y_Number],2)+
    						    pow(MeanY_Error[Strip_Y_Number]/pow(MeanY[Strip_Y_Number],2),2));

    // Calculate the Absolute Efficiency
    EfficiencyY[Strip_Y_Number] = 100*AreaY[Strip_Y_Number]/Photon_Emission;
    EfficiencyY_Error[Strip_Y_Number] = 100*AreaY_Error[Strip_Y_Number]/Photon_Emission;
    
  }
  return;
}


void Fitting_Efficiency_Calcs_Total(){

  // Define the histogram and get data from the already made histogram
  h1 = (TH1D*)fin->Get("hEnergyDepositcal");      

  // Setup up new function
  f_Total = new TF1("hEnergyDepositcal","gaus(0)",Peak_Location-Peak_Width,Peak_Location+Peak_Width);
 
  // Set initial paramters to assist with initial conditions
  f_Total->SetParameter(0,1);
  f_Total->SetParameter(1,Peak_Location);
  f_Total->SetParameter(2,5);

  // Fit the histogram parameters
  h1->Fit("hEnergyDepositcal","RMNQ");

  // Get gaussian amplitude
  Total_Amplitude = f_Total->GetParameter(0);
  Total_Amplitude_Error = f_Total->GetParError(0);  
  
  // Get the gaussian offset
  Total_Mean = f_Total->GetParameter(1);
  Total_Mean_Error = f_Total->GetParError(1);
  
  // Get the gaussian standard deviation
  Total_Sigma = f_Total->GetParameter(2);
  Total_Sigma_Error = f_Total->GetParError(2);

  // Calculate the area
  Total_Area = sqrt(2*Pi)*Total_Amplitude*Total_Sigma/Channel_Energy_Width;
  Total_Area_Error = sqrt(2*Pi)*sqrt(pow((Total_Sigma*Total_Amplitude_Error),2)+
				     pow((Total_Amplitude*Total_Sigma_Error),2))/Channel_Energy_Width;

  // Calculate the resolution
  Total_Resolution = 100.0*2.3548*Total_Sigma/Total_Mean;
  Total_Resolution_Error = 2.3545*sqrt(pow(Total_Sigma_Error/Total_Mean,2)+
   				    pow(Total_Mean_Error/pow(Total_Mean,2),2));

  // Calculate the Absolute Efficiency
  Total_Efficiency = 100*Total_Area/Photon_Emission;
  Total_Efficiency_Error = 100*Total_Area_Error/Photon_Emission;
 
  return;
}

void Channel_Sum_Area_Efficiency(){

  // SUM OVER ALL X CHANNELS FOR COMPARISON
  for(int Strip_X_Number=0; Strip_X_Number<Number_Of_Channels; Strip_X_Number++){
    Total_Area_From_Strips += AreaX[Strip_X_Number];
  }

  Total_Efficiency_From_Strips=100*Total_Area_From_Strips/Photon_Emission;
  
  //Efficiency Calculations from summing over strips
  cout << "Efficiency from summing over X Simulated Strips" << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "Efficiency from summed X strips: " << Total_Efficiency_From_Strips << endl << endl;

  //Efficiency calulations from hEnergyDepositcal
  cout << "Efficiency from hEnergydepositcal" << endl;
  cout << "---------------------------------" << endl;
  cout << "Efficiency from hEnergyDepostical: " << Total_Efficiency << endl << endl;


  // Total_Efficiency_From_Strips = 100*Total_Area_From_Strips/Photon_Emission;
  //cout << "Total Area from Total: " << Total_Area << endl;
  //cout << "Total Area from Strips: " << Total_Area_From_Strips << endl;
  
  //cout << "Total Efficiency from Total: " << Total_Efficiency << endl;
  //cout << "Total Efficiency From Strips: " << Total_Efficiency_From_Strips << endl;


}













//------------------------------------------------------------------------------------------------------------
//
//                                   GeGI ANALYSIS
//
//------------------------------------------------------------------------------------------------------------

// This function generates the histograms used for analyzing the GeGI Data. This function uses data from the
// root file generated from rootconvert, which converts raw Binary or ASCII format GeGI data into root trees.
// This function in turn takes each individual event and places it ito the correct histogram
void Make_GeGI_Histograms(){

  // Get Tree
  TTree *GeGI_Data = (TTree*)fin2->Get("data");

  // Define variables used in tree
  int Max_Events=40, gegi_eventnum[Max_Events], gegi_eventcount;
  double gegi_E[Max_Events], gegi_slow[Max_Events], gegi_fastpred[Max_Events], gegi_fastsucc[Max_Events];
  double gegi_posx, gegi_posy;
  int gegi_fastT[Max_Events], gegi_channel[Max_Events];
  long long gegi_globalT[Max_Events];
  
  // Give address of branches to GeGI data
  GeGI_Data->SetBranchAddress("gegi_eventcount",&gegi_eventcount);
  GeGI_Data->SetBranchAddress("gegi_E",&gegi_E);
  GeGI_Data->SetBranchAddress("gegi_globalT",&gegi_globalT);
  GeGI_Data->SetBranchAddress("gegi_eventnum",&gegi_eventnum);
  GeGI_Data->SetBranchAddress("gegi_fastT",&gegi_fastT);
  GeGI_Data->SetBranchAddress("gegi_channel",&gegi_channel);
  GeGI_Data->SetBranchAddress("gegi_slow",&gegi_slow);
  GeGI_Data->SetBranchAddress("gegi_fastpred",&gegi_fastpred);
  GeGI_Data->SetBranchAddress("gegi_fastsucc",&gegi_fastsucc);
  GeGI_Data->SetBranchAddress("gegi_posx",&gegi_posx);
  GeGI_Data->SetBranchAddress("gegi_posy",&gegi_posy);
 
  // Setup this unfilled histograms and set ranges derining x and y channels
  // NOTE: If you find that x and y are flipped compared to simulation, CHANGE THIS HERE
  // NOTE: NOTE THIS SECTION NOW REFLECTS PROPER X AND Y FOR GeGI DATA. SIMULATION NEEDS TO CHANGE
  int Channel_Number_Adjust=0;
  for(int i=0;i<32;i++){
    if(i<16){
      GeGI_Histograms[i] = new TH1D(Form(GeGI_Strip_Y_Name,i),Form(GeGI_Strip_Y_Name,i),6000,0,3000);
    } else{
      Channel_Number_Adjust = i-16;
      GeGI_Histograms[i] = new TH1D(Form(GeGI_Strip_X_Name,i),Form(GeGI_Strip_X_Name,i),6000,0,3000);
    }
  }

  GeGI_Total_X = new TH1D("GeGI_Total_X","GeGI Total X",6000,0,3000); 
  GeGI_Total_Y = new TH1D("GeGI_Total_Y","GeGI Total Y",6000,0,3000); 
  //GeGI_Image = new TH2F("GeGI_Image", "GeGI_Image",160,-40,40,160,-40,40);
  GeGI_Image = new TH2F("GeGI_Image","GeGI_Image",16,-0,16,16,0,16);
  
  // Get total number of entries in the tree
  int nevents = GeGI_Data->GetEntries();

  // Setup stuff for imaging histograms
  int Number_of_ROI = 10;
  double ROI[10][10];

  //Fill histograms with data from tree including putting them into the propper channel
  for (int i=0; i<nevents; i++){
    GeGI_Data->GetEntry(i);
    double Temp_Energy_X=0;
    double Temp_Energy_Y=0;
    double Temp_Pos_X=0;
    double Temp_Pos_Y=0;
    double Strip_Width=5;

    // Fill histograms with every entry, not just the first one
    for(int j=0; j<gegi_eventcount; j++){
      GeGI_Histograms[gegi_channel[j]]->Fill(gegi_slow[j]);
      
      if(gegi_channel[j]>=16 && gegi_channel[j]<32){
      	Temp_Energy_X += gegi_slow[j]; 
	if(gegi_eventcount<=2){
	//Temp_Pos_X = ((gegi_channel[j]-16)+gegi_fastsucc[j]/(gegi_fastsucc[j]+gegi_fastpred[j])-8)*Strip_Width;
	Temp_Pos_X = (gegi_channel[j]-16);
	//Temp_Pos_X = gegi_posx;
	}
	
      }
      else if (gegi_channel[j]<16){
      	Temp_Energy_Y += gegi_slow[j];

	if(gegi_eventcount<=2){
	//Temp_Pos_Y = (gegi_channel[j]+gegi_fastsucc[j]/(gegi_fastsucc[j]+gegi_fastpred[j])-8)*Strip_Width;
	Temp_Pos_Y = (gegi_channel[j]);
	//Temp_Pos_Y = gegi_posy;
	}
      }
      
	  GeGI_Image->Fill(Temp_Pos_X,Temp_Pos_Y);
    }
    
    GeGI_Total_X->Fill(Temp_Energy_X);
    GeGI_Total_Y->Fill(Temp_Energy_Y);
  }

  return;
}


void Fitting_GeGI_Efficiency_Calcs_X(){

  //Define the histogram and get data from the already made histogram
  //h2 = (TH1D*)fin->Get(Form(GeGI_Strip_Y_Name,GeGI_Strip_Y_Number)); 
  for(int GeGI_Strip_X_Number = 0;GeGI_Strip_X_Number<16;GeGI_Strip_X_Number++){
    // Setup up new function
    GeGI_f1[GeGI_Strip_X_Number] = new TF1(Form(GeGI_Strip_X_Name,GeGI_Strip_X_Number),"gaus(0)+pol2(3)",
					   Peak_Location-Peak_Width,Peak_Location+Peak_Width);
 
    // Set initial paramters to assist with initial conditions
    GeGI_f1[GeGI_Strip_X_Number]->SetParameter(0,100);
    GeGI_f1[GeGI_Strip_X_Number]->SetParameter(1,Peak_Location+1);
    GeGI_f1[GeGI_Strip_X_Number]->SetParameter(2,1);

    // Fit the histogram parameters
    GeGI_Histograms[GeGI_Strip_X_Number+16]->Fit(Form(GeGI_Strip_X_Name,GeGI_Strip_X_Number),"RMNQ");
  
    // Get gaussian amplitude
    GeGI_AmplitudeX[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParameter(0);
    GeGI_AmplitudeX_Error[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParError(0);
    
    // Get the gaussian offset
    GeGI_MeanX[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParameter(1);
    GeGI_MeanX_Error[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParError(1);
 
    // Get the gaussian standard deviation
    GeGI_SigmaX[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParameter(2);
    GeGI_SigmaX_Error[GeGI_Strip_X_Number] = GeGI_f1[GeGI_Strip_X_Number]->GetParError(2);

    // Calculate the area
    GeGI_AreaX[GeGI_Strip_X_Number] = sqrt(2*Pi)*GeGI_AmplitudeX[GeGI_Strip_X_Number]*GeGI_SigmaX[GeGI_Strip_X_Number]/Channel_Energy_Width;
    GeGI_AreaX_Error[GeGI_Strip_X_Number] = sqrt(2*Pi)*sqrt(pow((GeGI_SigmaX[GeGI_Strip_X_Number]*
								 GeGI_AmplitudeX_Error[GeGI_Strip_X_Number]),2)+
							    pow((GeGI_AmplitudeX[GeGI_Strip_X_Number]*
								 GeGI_SigmaX_Error[GeGI_Strip_X_Number]),2))/Channel_Energy_Width;
    
    //Calculate the resolution
    GeGI_ResolutionX[GeGI_Strip_X_Number] = 100.0*2.3548*GeGI_SigmaX[GeGI_Strip_X_Number]/GeGI_MeanX[GeGI_Strip_X_Number];
    GeGI_ResolutionX_Error[GeGI_Strip_X_Number] = 2.3545*sqrt(pow(GeGI_SigmaX_Error[GeGI_Strip_X_Number]/
								  GeGI_MeanX[GeGI_Strip_X_Number],2)+
							      pow(GeGI_MeanX_Error[GeGI_Strip_X_Number]/
								  pow(GeGI_MeanX[GeGI_Strip_X_Number],2),2));
    // Calculate the Absolute Efficiency
    GeGI_EfficiencyX[GeGI_Strip_X_Number] = 100*GeGI_AreaX[GeGI_Strip_X_Number]/Photon_Emission;
    GeGI_EfficiencyX_Error[GeGI_Strip_X_Number] = 100*GeGI_AreaX_Error[GeGI_Strip_X_Number]/Photon_Emission;
    
    cout << " " << GeGI_SigmaX[GeGI_Strip_X_Number] << "," << GeGI_SigmaX_Error[GeGI_Strip_X_Number] << endl;

    // MANUAL CALCULATION OF AREA/EFFICIENCY USING COVELL METHOD
    // int Area_Lower_Bound;
    // int Area_Upper_Bound;
    // int Number_of_Channels_Photopeak=20; // number of channels covering photopeak
    // int Number_of_Channels_Background=10; // number of channels covering upper and lower background regions
    // double Bin_Width = 0.5; // Bin Width in energy, this is defined by the total number of bins and the energy range covered

    
    // Area_Lower_Bound = Peak_Location/Bin_Width-Number_of_Channels_Photopeak/2;
    // Area_Upper_Bound = Peak_Location/Bin_Width+Number_of_Channels_Photopeak/2;
    
    // GeGI_AreaX_Manual[GeGI_Strip_X_Number] = (GeGI_Histograms[GeGI_Strip_X_Number+16]->Integral(Area_Lower_Bound,Area_Upper_Bound))-
    //   Number_of_Channels_Photopeak*(GeGI_Histograms[GeGI_Strip_X_Number+16]->Integral(Area_Lower_Bound-Number_of_Channels_Background,Area_Lower_Bound-1)+
    //   GeGI_Histograms[GeGI_Strip_X_Number+16]->Integral(Area_Upper_Bound+1,Area_Upper_Bound+Number_of_Channels_Background))/
    //   (2*Number_of_Channels_Background);

    // cout << GeGI_
    //  Strip_X_Number << "Area Manual: "  << GeGI_AreaX_Manual[GeGI_Strip_X_Number] << "Fitting: " << GeGI_AreaX[GeGI_Strip_X_Number] << " Difference: " <<  GeGI_AreaX_Manual[GeGI_Strip_X_Number]/ GeGI_AreaX[GeGI_Strip_X_Number] << endl;
    

  }
}


void Fitting_GeGI_Efficiency_Calcs_Y(){

  //Define the histogram and get data from the already made histogram
  //h2 = (TH1D*)fin->Get(Form(GeGI_Strip_Y_Name,GeGI_Strip_Y_Number)); 
  for(int GeGI_Strip_Y_Number = 0;GeGI_Strip_Y_Number<16;GeGI_Strip_Y_Number++){
    // Setup up new function
    GeGI_f2[GeGI_Strip_Y_Number] = new TF1(Form(GeGI_Strip_Y_Name,GeGI_Strip_Y_Number),"gaus(0)+pol2(3)",
					   Peak_Location-Peak_Width,Peak_Location+Peak_Width);
 
    // Set initial paramters to assist with initial conditions
    GeGI_f2[GeGI_Strip_Y_Number]->SetParameter(0,1);
    GeGI_f2[GeGI_Strip_Y_Number]->SetParameter(1,Peak_Location);
    GeGI_f2[GeGI_Strip_Y_Number]->SetParameter(2,5);

    // Fit the histogram parameters
    GeGI_Histograms[GeGI_Strip_Y_Number]->Fit(Form(GeGI_Strip_Y_Name,GeGI_Strip_Y_Number),"RMNQ");
  
    // Get gaussian amplitude
    GeGI_AmplitudeY[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParameter(0);
    GeGI_AmplitudeY_Error[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParError(0);
    
    // Get the gaussian offset
    GeGI_MeanY[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParameter(1);
    GeGI_MeanY_Error[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParError(1);
 
    // Get the gaussian standard deviation
    GeGI_SigmaY[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParameter(2);
    GeGI_SigmaY_Error[GeGI_Strip_Y_Number] = GeGI_f2[GeGI_Strip_Y_Number]->GetParError(2);

    // Calculate the area
    GeGI_AreaY[GeGI_Strip_Y_Number] = sqrt(2*Pi)*GeGI_AmplitudeY[GeGI_Strip_Y_Number]*GeGI_SigmaY[GeGI_Strip_Y_Number]/Channel_Energy_Width;
    GeGI_AreaY_Error[GeGI_Strip_Y_Number] = sqrt(2*Pi)*sqrt(pow((GeGI_SigmaY[GeGI_Strip_Y_Number]*
								 GeGI_AmplitudeY_Error[GeGI_Strip_Y_Number]),2)+
							    pow((GeGI_AmplitudeY[GeGI_Strip_Y_Number]*
								 GeGI_SigmaY_Error[GeGI_Strip_Y_Number]),2))/Channel_Energy_Width;
    //Calculate the resolution
    GeGI_ResolutionY[GeGI_Strip_Y_Number] = 100.0*2.3548*GeGI_SigmaY[GeGI_Strip_Y_Number]/GeGI_MeanX[GeGI_Strip_Y_Number];
    GeGI_ResolutionY_Error[GeGI_Strip_Y_Number] = 2.3545*sqrt(pow(GeGI_SigmaY_Error[GeGI_Strip_Y_Number]/
								  GeGI_MeanY[GeGI_Strip_Y_Number],2)+
							      pow(GeGI_MeanY_Error[GeGI_Strip_Y_Number]/
								  pow(GeGI_MeanY[GeGI_Strip_Y_Number],2),2));
    
    // Calculate the Absolute Efficiency
    GeGI_EfficiencyY[GeGI_Strip_Y_Number] = 100*GeGI_AreaY[GeGI_Strip_Y_Number]/Photon_Emission;
    GeGI_EfficiencyY_Error[GeGI_Strip_Y_Number] = 100*GeGI_AreaY_Error[GeGI_Strip_Y_Number]/Photon_Emission;

    cout << " " << GeGI_SigmaY[GeGI_Strip_Y_Number] << "," << GeGI_SigmaY_Error[GeGI_Strip_Y_Number] <<  endl;
  }
}

void Fitting_GeGI_Efficiency_Calcs_Total_X(){
  // Set up new function
  GeGI_f_Total_X = new TF1("GeGI_Total_Eff_X","gaus(0)",Peak_Location-Peak_Width,Peak_Location+Peak_Width);

  // Set initial paramters to assist with initial conditions
  GeGI_f_Total_X->SetParameter(0,1);
  GeGI_f_Total_X->SetParameter(1,Peak_Location);
  GeGI_f_Total_X->SetParameter(2,5);

  // Fit the histogram parameters
  //GeGI_Total_X->Fit("GeGI_Total_Eff_X","RMNQ");
  GeGI_Total_X->Fit("GeGI_Total_Eff_X","RMNQ");
  
  // Get gaussian amplitude
  GeGI_Total_Amplitude_X = GeGI_f_Total_X->GetParameter(0);
  GeGI_Total_Amplitude_X_Error = GeGI_f_Total_X->GetParError(0);

  // Get the gaussian offset
  GeGI_Total_Mean_X = GeGI_f_Total_X->GetParameter(1);
  GeGI_Total_Mean_X_Error = GeGI_f_Total_X->GetParError(1);

  // Get the gaussian standard deviation
  GeGI_Total_Sigma_X = GeGI_f_Total_X->GetParameter(2);
  GeGI_Total_Sigma_X_Error = GeGI_f_Total_X->GetParError(2);

  // Calculate the area
  GeGI_Total_Area_X = sqrt(2*Pi)*GeGI_Total_Amplitude_X*GeGI_Total_Sigma_X/Channel_Energy_Width;
  GeGI_Total_Area_X_Error = sqrt(2*Pi)*sqrt(pow((GeGI_Total_Sigma_X*GeGI_Total_Amplitude_X_Error),2)+
							    pow((GeGI_Total_Amplitude_X*GeGI_Total_Sigma_X_Error),2))/Channel_Energy_Width;

  //Calculate the resolution
  GeGI_Total_Resolution_X = 100.0*2.3548*GeGI_Total_Sigma_X/GeGI_Total_Mean_X;
  GeGI_Total_Resolution_X_Error = 2.3545*sqrt(pow(GeGI_Total_Sigma_X_Error/GeGI_Total_Mean_X,2)+pow(GeGI_Total_Mean_X_Error/
												    pow(GeGI_Total_Mean_X,2),2));

  // Calculate the Absolute Efficiency
  GeGI_Total_Efficiency_X = 100*GeGI_Total_Area_X/Photon_Emission;
  GeGI_Total_Efficiency_X_Error = 100*GeGI_Total_Area_X_Error/Photon_Emission;

  // Output results to command line
  cout << "GeGI Analysis for using Sum of X Channel Data" << endl;
  cout << "---------------------------------------------" << endl;
  cout << "GeGI Total Efficiency: " << GeGI_Total_Efficiency_X <<  endl;
  cout << "GeGI Total Resolution:  " <<  GeGI_Total_Resolution_X << endl << endl;

}

void Fitting_GeGI_Efficiency_Calcs_Total_Y(){
  /// Setup new function
  GeGI_f_Total_Y = new TF1("GeGI_Total_Eff_Y","gaus(0)",Peak_Location-Peak_Width,Peak_Location+Peak_Width);

  // Set initial paramters to assist with initial conditions
  GeGI_f_Total_Y->SetParameter(0,1);
  GeGI_f_Total_Y->SetParameter(1,Peak_Location);
  GeGI_f_Total_Y->SetParameter(2,5);

  // Fit the histogram parameters
  GeGI_Total_Y->Fit("GeGI_Total_Eff_Y","RMNQ");

  // Get gaussian amplitude
  GeGI_Total_Amplitude_Y = GeGI_f_Total_Y->GetParameter(0);
  GeGI_Total_Amplitude_Y_Error = GeGI_f_Total_Y->GetParError(0);

  // Get the gaussian offset
  GeGI_Total_Mean_Y = GeGI_f_Total_Y->GetParameter(1);
  GeGI_Total_Mean_Y_Error = GeGI_f_Total_Y->GetParError(1);

  // Get the gaussian standard deviation
  GeGI_Total_Sigma_Y = GeGI_f_Total_Y->GetParameter(2);
  GeGI_Total_Sigma_Y_Error = GeGI_f_Total_Y->GetParError(2);

  // Calculate the area
  GeGI_Total_Area_Y = sqrt(2*Pi)*GeGI_Total_Amplitude_Y*GeGI_Total_Sigma_Y/Channel_Energy_Width;
  GeGI_Total_Area_Y_Error = sqrt(2*Pi)*sqrt(pow((GeGI_Total_Sigma_Y*GeGI_Total_Amplitude_Y_Error),2)+
							    pow((GeGI_Total_Amplitude_Y*GeGI_Total_Sigma_Y_Error),2))/
                                                            Channel_Energy_Width;

  //Calculate the resolution
  GeGI_Total_Resolution_Y = 100.0*2.3548*GeGI_Total_Sigma_Y/GeGI_Total_Mean_Y;
  GeGI_Total_Resolution_Y_Error = 2.3545*sqrt(pow(GeGI_Total_Sigma_Y_Error/GeGI_Total_Mean_Y,2)+pow(GeGI_Total_Mean_Y_Error/
												    pow(GeGI_Total_Mean_Y,2),2));

  // Calculate the Absolute Efficiency
  GeGI_Total_Efficiency_Y = 100*GeGI_Total_Area_Y/Photon_Emission;
  GeGI_Total_Efficiency_Y_Error = 100*GeGI_Total_Area_Y_Error/Photon_Emission;

  // Output results to command line
  cout << "GeGI Analysis for using Sum of Y Channel Data" << endl;
  cout << "---------------------------------------------" << endl;
  cout << "GeGI Total Efficiency: " << GeGI_Total_Efficiency_Y << " " <<  endl;
  cout << "GeGI Total Resolution:  " <<  GeGI_Total_Resolution_Y << " " <<  endl << endl;



} 


void Fitting_GeGI_Efficiency_Manual_Calcs_X(){

  return;
}

void Fitting_GeGi_Efficiency_Manual_Calcs_Y(){

  return;
}




//------------------------------------------------------------------------------------------------------------
//
//                                   PLOTTING FUNCTIONS
//
//------------------------------------------------------------------------------------------------------------

//----------------------
// SIMULATED EFFICINCIES
//----------------------

void Plotting_Histograms(){

  // REMOVE COMMENTS TO INACT THAT PLOT
  //Efficiency_Plot_X(); // Plotting efficiency of simulated x channels
  //Efficiency_Plot_Y(); // Plotting efficiency of simulated y channels
  //Efficiency_Plot_XY(); // Plotting efficiency of simulated x and y strip efficiencies on single canvas
  //GeGI_Efficiency_Plot_X();
  //GeGI_Efficiency_Plot_Y();
  //GeGI_Efficiency_Plot_XY();
  //Comparison_Efficiency_Plot_X();
  //Comparison_Efficiency_Plot_Y();

  //Simulated_Multiplot_Strips_X();
  //Simulated_Miltiplot_Strips_Y();
  GeGI_Multiplot_Strips_X();
  GeGI_Multiplot_Strips_Y();

 
  //GeGI_Image->Draw("colz");

  return;
}


void Efficiency_Plot_X(){

  // Setup Canvas for Ploting Efficiencies
  auto Efficiency_Plot_X = new TCanvas("Efficiency_Plot_X","X Strip Efficiencies",1000,700);
  Efficiency_Plot_X->SetFillColor(0);
  Efficiency_Plot_X->SetGrid();
  Efficiency_Plot_X->GetFrame()->SetFillColor(21);
  Efficiency_Plot_X->GetFrame()->SetBorderSize(14);

  auto EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyX,0,EfficiencyX_Error);
  EfficienciesX->SetName("EfficienciesX");
  EfficienciesX->SetTitle("X Strip Absolute Efficiencies");
  EfficienciesX->SetMarkerColor(4);
  EfficienciesX->SetMarkerStyle(21);
  EfficienciesX->Draw("ALP");
  EfficienciesX->GetXaxis()->SetTitle("Channel Number");
  EfficienciesX->GetXaxis()->CenterTitle();
  EfficienciesX->GetYaxis()->SetTitle("Absolute Efficiency");
  EfficienciesX->GetYaxis()->CenterTitle();
  
  return;
}

void Efficiency_Plot_Y(){

  // Setup Canvas for Ploting Efficiencies
  auto Efficiency_Plot_Y = new TCanvas("Simulated Efficiency_Plot_Y","A Simple Graph with error bars",1024,768);
  Efficiency_Plot_Y->SetFillColor(0);
  Efficiency_Plot_Y->SetGrid();
  Efficiency_Plot_Y->GetFrame()->SetFillColor(21);
  Efficiency_Plot_Y->GetFrame()->SetBorderSize(14);

  auto EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyY,0,EfficiencyY_Error);
  EfficienciesY->SetName("EfficienciesY");
  EfficienciesY->SetTitle("Simulated GeGI Y Strip Absolute Efficiencies");
  EfficienciesY->SetMarkerStyle(33);
  EfficienciesY->SetMarkerColor(2);
  EfficienciesY->Draw("ALP");
  EfficienciesY->GetXaxis()->SetTitle("Channel Number");
  EfficienciesY->GetXaxis()->CenterTitle();
  EfficienciesY->GetYaxis()->SetTitle("Absolute Efficiency");
  EfficienciesY->GetYaxis()->CenterTitle();
  
  return;
}


void Efficiency_Plot_XY(){

  // Setup Canvas for Ploting Efficiencies
  auto c1 = new TCanvas("c1","Simulated GeGI Efficincies X and Y Strips",1000,700);
  c1->SetFillColor(0);
  c1->SetGrid();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(14);

  auto EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyX,0,EfficiencyX_Error);
  EfficienciesX->SetName("EfficienciesX");
  EfficienciesX->SetTitle("Simulated GeGI X and Y Strip Absolute Efficiencies");
  EfficienciesX->SetMarkerColor(4);
  EfficienciesX->SetMarkerStyle(21);
  EfficienciesX->Draw("ALP");
  auto EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyY,0,EfficiencyY_Error);
  EfficienciesX->SetName("EfficienciesY");
  EfficienciesY->SetMarkerStyle(33);
  EfficienciesY->SetMarkerColor(2);
  EfficienciesY->Draw("LP");
  EfficienciesX->GetXaxis()->SetTitle("Channel Number");
  EfficienciesX->GetXaxis()->CenterTitle();
  EfficienciesX->GetYaxis()->SetTitle("Absolute Efficiency");
  EfficienciesX->GetYaxis()->CenterTitle();
  EfficienciesX->GetYaxis()->SetRangeUser(0,1.5e-4);
  auto legend = new TLegend(0.75,0.75,0.88,0.88);
  legend->SetHeader("Legend","C"); // option "C" allows to center the header
  legend->AddEntry(EfficienciesX,"X Strips","ep");
  legend->AddEntry(EfficienciesY,"Y Strips","ep");
  legend->Draw();

  return;
}

// Plot histograms and fittings on 4x4 grid for X channels
void Simulated_Multiplot_Strips_X(){
  auto Sim_Multi_Strip_X=new TCanvas("Sim_Multi_Strip_X","Simulated X Histograms",1600,900);
  
  Sim_Multi_Strip_X->Divide(4,4);

  for(int i=0;i<16;i++){
    Sim_Multi_Strip_X->cd(i+1);
    h1 = (TH1D*)fin->Get(Form(Strip_X_Name,i));
    h1->GetXaxis()->SetRange(1280,1360);
    h1->GetXaxis()->SetTitle("Energy (keV)");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->CenterTitle();
    h1->Draw();
    f1[i]->Draw("same");
  }

  
  return;
}


// Plot histograms and fittings on 4x4 grid for Y channels
void Simulated_Miltiplot_Strips_Y(){
  auto Sim_Multi_Strip_Y=new TCanvas("Sim_Multi_Strip_Y","Simulated Y Histograms",1600,900);
  
  Sim_Multi_Strip_Y->Divide(4,4);

  for(int i=0;i<16;i++){
    Sim_Multi_Strip_Y->cd(i+1);
    h1 = (TH1D*)fin->Get(Form(Strip_Y_Name,i));
    h1->GetXaxis()->SetRange(1280,1360);
    h1->GetXaxis()->SetTitle("Energy (keV)");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->CenterTitle();
    h1->Draw();
    f2[i]->Draw("same");
  }
  return;
}


//----------------------
// MEASURED EFFICIENCIES
//----------------------

void GeGI_Efficiency_Plot_X(){

  // Setup Canvas for Ploting Efficiencies
  auto GeGI_Efficiency_Plot_X = new TCanvas("GeGI_Efficiency_Plot_X","GeGI X Strip Efficiencies",1000,700);
  GeGI_Efficiency_Plot_X->SetFillColor(0);
  GeGI_Efficiency_Plot_X->SetGrid();
  GeGI_Efficiency_Plot_X->GetFrame()->SetFillColor(21);
  GeGI_Efficiency_Plot_X->GetFrame()->SetBorderSize(14);

  auto GeGI_EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyX,0,GeGI_EfficiencyX_Error);
  GeGI_EfficienciesX->SetName("EfficienciesX");
  GeGI_EfficienciesX->SetTitle("GeGI X Strip Absolute Efficiencies");
  GeGI_EfficienciesX->SetMarkerColor(4);
  GeGI_EfficienciesX->SetMarkerStyle(21);
  GeGI_EfficienciesX->Draw("ALP");
  GeGI_EfficienciesX->GetXaxis()->SetTitle("Channel Number");
  GeGI_EfficienciesX->GetXaxis()->CenterTitle();
  GeGI_EfficienciesX->GetYaxis()->SetTitle("Absolute Efficiency");
  GeGI_EfficienciesX->GetYaxis()->CenterTitle();
  
  return;
}

void GeGI_Efficiency_Plot_Y(){

  // Setup Canvas for Ploting Efficiencies
  auto  GeGI_Efficiency_Plot_Y = new TCanvas(" GeGI_Efficiency_Plot_Y","GeGI Y Strip Efficiencies",1000,700);
  GeGI_Efficiency_Plot_Y->SetFillColor(0);
  GeGI_Efficiency_Plot_Y->SetGrid();
  GeGI_Efficiency_Plot_Y->GetFrame()->SetFillColor(21);
  GeGI_Efficiency_Plot_Y->GetFrame()->SetBorderSize(14);

  auto GeGI_EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyY,0,GeGI_EfficiencyY_Error);
  GeGI_EfficienciesY->SetName("EfficienciesY");
  GeGI_EfficienciesY->SetTitle("GeGI Y Strip Absolute Efficiencies");
  GeGI_EfficienciesY->SetMarkerStyle(33);
  GeGI_EfficienciesY->SetMarkerColor(2);
  GeGI_EfficienciesY->Draw("ALP");
  GeGI_EfficienciesY->GetXaxis()->SetTitle("Channel Number");
  GeGI_EfficienciesY->GetXaxis()->CenterTitle();
  GeGI_EfficienciesY->GetYaxis()->SetTitle("Absolute Efficiency");
  GeGI_EfficienciesY->GetYaxis()->CenterTitle();
  
  return;
}


void GeGI_Efficiency_Plot_XY(){

  // Setup Canvas for Ploting Efficiencies
  auto GeGI_Efficiencies_Plot_XY = new TCanvas("GeGI_XY","A Simple Graph with error bars",1000,700);
  GeGI_Efficiencies_Plot_XY->SetFillColor(0);
  GeGI_Efficiencies_Plot_XY->SetGrid();
  GeGI_Efficiencies_Plot_XY->GetFrame()->SetFillColor(21);
  GeGI_Efficiencies_Plot_XY->GetFrame()->SetBorderSize(14);

  auto GeGI_EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyX,0,GeGI_EfficiencyX_Error);
  GeGI_EfficienciesX->SetName("GeGI_EfficienciesX");
  GeGI_EfficienciesX->SetTitle("X and Y Strip Absolute Efficiencies");
  GeGI_EfficienciesX->SetMarkerColor(4);
  GeGI_EfficienciesX->SetMarkerStyle(21);
  GeGI_EfficienciesX->Draw("ALP");
  auto GeGI_EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyY,0,GeGI_EfficiencyY_Error);
  GeGI_EfficienciesX->SetName("GeGI_EfficienciesY");
  GeGI_EfficienciesY->SetMarkerStyle(33);
  GeGI_EfficienciesY->SetMarkerColor(2);
  GeGI_EfficienciesY->Draw("LP");
  GeGI_EfficienciesX->GetXaxis()->SetTitle("Channel Number");
  GeGI_EfficienciesX->GetXaxis()->CenterTitle();
  GeGI_EfficienciesX->GetYaxis()->SetTitle("Absolute Efficiency");
  GeGI_EfficienciesX->GetYaxis()->CenterTitle();
  GeGI_EfficienciesX->GetYaxis()->SetRangeUser(0,3e-4);
  auto legend = new TLegend(0.75,0.75,0.88,0.88);
  legend->SetHeader("Legend","C"); // option "C" allows to center the header
  legend->AddEntry(GeGI_EfficienciesX,"X Strips","ep");
  legend->AddEntry(GeGI_EfficienciesY,"Y Strips","ep");
  legend->Draw();

  return;
}


// Plot histograms and fittings on 4x4 grid for MEASURED X channels
void GeGI_Multiplot_Strips_X(){
  auto GeGI_Multi_Strip_X=new TCanvas("GeGI_Multi_Strip_X","Measured X Histograms",1600,900);
  
  GeGI_Multi_Strip_X->Divide(4,4);

  for(int i=0;i<16;i++){
    GeGI_Multi_Strip_X->cd(i+1);
    h1 = GeGI_Histograms[i+16];
    h1->GetXaxis()->SetRange((Peak_Location-20)/Channel_Energy_Width,(Peak_Location+20)/Channel_Energy_Width);
    h1->GetXaxis()->SetTitle("Energy (keV)");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->CenterTitle();
    h1->Draw();
    GeGI_f1[i]->Draw("same");
  }

  
  return;
}


// Plot histograms and fittings on 4x4 grid for MEASURED Y channels
void GeGI_Multiplot_Strips_Y(){
  auto GeGI_Multi_Strip_Y=new TCanvas("GeGI_Multi_Strip_Y","Measured Y Histograms",1600,900);
  
  GeGI_Multi_Strip_Y->Divide(4,4);

  for(int i=0;i<16;i++){
    GeGI_Multi_Strip_Y->cd(i+1);
    h1 = GeGI_Histograms[i];
    h1->GetXaxis()->SetRange((Peak_Location-20)/Channel_Energy_Width,(Peak_Location+20)/Channel_Energy_Width);
    h1->GetXaxis()->SetTitle("Energy (keV)");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->CenterTitle();
    h1->Draw();
    GeGI_f2[i]->Draw("same");
  }

  
  return;
}
//--------------------------------
// COMPARING MEASURED VS SIMULATED
//--------------------------------

void Comparison_Efficiency_Plot_X(){

  // Setup canvas for plotting efficiencies
  auto Comp_Eff_Plot_X = new TCanvas("Comp_Eff_X","Efficiency Comparision Between Measured and Simulated X Channels");
  Comp_Eff_Plot_X->SetFillColor(0);
  Comp_Eff_Plot_X->SetGrid();
  Comp_Eff_Plot_X->GetFrame()->SetFillColor(21);
  Comp_Eff_Plot_X->GetFrame()->SetBorderSize(14);

  auto Comp_EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyX,0,EfficiencyX_Error);
  Comp_EfficienciesX->SetName("GeGI_EfficienciesX");
  Comp_EfficienciesX->SetTitle("GeGI Measured and Simulated X Channel Efficiencies");
  Comp_EfficienciesX->SetMarkerColor(4);
  Comp_EfficienciesX->SetMarkerStyle(21);
  Comp_EfficienciesX->Draw("ALP");
  auto Comp_GeGI_EfficienciesX = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyX,0,GeGI_EfficiencyX_Error);
  Comp_GeGI_EfficienciesX->SetName("GeGI_EfficienciesX");
  Comp_GeGI_EfficienciesX->SetMarkerStyle(33);
  Comp_GeGI_EfficienciesX->SetMarkerColor(2);
  Comp_GeGI_EfficienciesX->Draw("LP");
  Comp_EfficienciesX->GetXaxis()->SetTitle("Channel Number");
  Comp_EfficienciesX->GetXaxis()->CenterTitle();
  Comp_EfficienciesX->GetYaxis()->SetTitle("Absolute Efficiency");
  Comp_EfficienciesX->GetYaxis()->CenterTitle();
  Comp_EfficienciesX->GetYaxis()->SetRangeUser(0,3e-4);
  auto legend = new TLegend(0.75,0.75,0.88,0.88);
  legend->SetHeader("Legend","C"); // option "C" allows to center the header
  legend->AddEntry(Comp_EfficienciesX,"Simulated Strips","ep");
  legend->AddEntry(Comp_GeGI_EfficienciesX,"Measured Strips","ep");
  legend->Draw();

  return;


}


void Comparison_Efficiency_Plot_Y(){

  // Setup canvas for plotting efficiencies
  auto Comp_Eff_Plot_Y = new TCanvas("Comp_Eff_Y","Efficiency Comparision Between Measured and Simulated Y Channels");
  Comp_Eff_Plot_Y->SetFillColor(0);
  Comp_Eff_Plot_Y->SetGrid();
  Comp_Eff_Plot_Y->GetFrame()->SetFillColor(21);
  Comp_Eff_Plot_Y->GetFrame()->SetBorderSize(14);

  auto Comp_EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,EfficiencyY,0,EfficiencyY_Error);
  Comp_EfficienciesY->SetName("GeGI_EfficienciesX");
  Comp_EfficienciesY->SetTitle("GeGI Measured and Simulated Y Channel Efficiencies");
  Comp_EfficienciesY->SetMarkerColor(4);
  Comp_EfficienciesY->SetMarkerStyle(21);
  Comp_EfficienciesY->Draw("ALP");
  auto Comp_GeGI_EfficienciesY = new TGraphErrors(Number_Of_Channels,Strip_Number,GeGI_EfficiencyY,0,GeGI_EfficiencyY_Error);
  Comp_GeGI_EfficienciesY->SetName("GeGI_EfficienciesY");
  Comp_GeGI_EfficienciesY->SetMarkerStyle(33);
  Comp_GeGI_EfficienciesY->SetMarkerColor(2);
  Comp_GeGI_EfficienciesY->Draw("LP");
  Comp_EfficienciesY->GetXaxis()->SetTitle("Channel Number");
  Comp_EfficienciesY->GetXaxis()->CenterTitle();
  Comp_EfficienciesY->GetYaxis()->SetTitle("Absolute Efficiency");
  Comp_EfficienciesY->GetYaxis()->CenterTitle();
  Comp_EfficienciesY->GetYaxis()->SetRangeUser(0,3e-4);
  auto legend = new TLegend(0.75,0.75,0.88,0.88);
  legend->SetHeader("Legend","C"); // option "C" allows to center the header
  legend->AddEntry(Comp_EfficienciesY,"Simulated Strips","ep");
  legend->AddEntry(Comp_GeGI_EfficienciesY,"Measured Strips","ep");
  legend->Draw();

  return;


}


// MISC Functions

void SRM_Input(){

  ifstream source_data;
  

  source_data.open(SRMfilename);
  
  source_data>>SRM_Num_of_Isotopes;

  cout << "  Number of Isotopes in SRM: " << SRM_Num_of_Isotopes << endl << endl; 

  cout << "Source Number, Half-Life (days), Activity (Bq)" << endl;

  for(int i=0; i<SRM_Num_of_Isotopes;i++){
  
    source_data>>SRM_Source_HL_Activity[i][0];
    source_data>>SRM_Source_HL_Activity[i][1];

    cout << "Source " << i << ", " << SRM_Source_HL_Activity[i][0] << ", " << SRM_Source_HL_Activity[i][1] << endl ;

  }
  
  
  source_data>>SRM_Date.tm_year;
  source_data>>SRM_Date.tm_mon;
  source_data>>SRM_Date.tm_mday;
  source_data>>SRM_Date.tm_hour;
  SRM_Date.tm_min=0;
  SRM_Date.tm_sec=0;
  
  cout << endl;
  cout << "Reference Date" << endl;
  cout << "Year: " << SRM_Date.tm_year << " Month: " << SRM_Date.tm_mon << " Day: " << SRM_Date.tm_mday << "  Hour: "
       << SRM_Date.tm_hour << endl << endl;
      
  Run_Date.tm_year=2020;
  Run_Date.tm_mon=5;
  Run_Date.tm_mday=0;
  Run_Date.tm_hour=0;
  Run_Date.tm_min=0;
  Run_Date.tm_sec=0;

  cout << "Run Date:" << endl;
  cout << "Year: " << Run_Date.tm_year << " Month: " << Run_Date.tm_mon << " Day: " << Run_Date.tm_mday << "  Hour: "
       << Run_Date.tm_hour << endl << endl;

  Time_Difference = difftime(mktime(&Run_Date),mktime(&SRM_Date))/86400; // in days since halflife is in days

  cout << "Time between refence and run dates: " << Time_Difference << " days" <<  endl << endl;

  //for(int isotope_num=0;isotope_num<SRM_Num_of_Isotopes;isotope_num++){
    


    for(int i=0;i<SRM_Num_of_Isotopes;i++){

      int Number_Of_Emissions=-1; 
      source_data>>Number_Of_Emissions; // Get number of emissions

      cout << "Isotope Number: " << i << " Number of Emissions: " << Number_Of_Emissions << endl;
      cout << "Emission Energy, Probability, Probabilty Error" << endl;
      
      for(int j=0;j<Number_Of_Emissions;j++){

	source_data>>SRM_Emission[i][j][0]; // Emision Eenergy
	source_data>>SRM_Emission[i][j][1]; // Emission probability
	source_data>>SRM_Emission[i][j][2]; // Emission probability error

	cout << SRM_Emission[i][j][0] << ", " << SRM_Emission[i][j][1] << ", " << SRM_Emission[i][j][2] << endl;
      
      }
      
      cout << endl;
    }
    //}					     
}

