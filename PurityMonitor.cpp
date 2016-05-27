// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <sstream>

// C headers
#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <string.h>
#include <thread>
#include <pthread.h>
#include <unistd.h>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TApplication.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TFile.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TView3D.h"
#include "TTree.h"

// Drift parameter for E-field to v-Field conversion
// const float DriftPar_0 = 1.31195;
// const float DriftPar_1 = 1.87803e-7;
// const float DriftPar_2 = 10.9707;
// const float DriftPar_3 = 0.216244;

int main();
std::vector<TH2F*> ReadDataFile(std::string,unsigned);
std::vector<TH2F*> SimpleFiled(unsigned);
void DrawHistos(std::vector<TH2F*>);
std::vector<float> VelocityFunction(std::vector<float>);
double Simulation(std::vector<TH2F*>, std::vector<double>, std::vector<float>, double, float, unsigned, unsigned);
void SimpleSimulation(std::vector<TH2F*>, std::vector<double>, std::vector<float>, double, float, unsigned, unsigned);

int main()
{
    std::vector<double> ParticlePosition = {0.1,-1};
    std::vector<float> MonitorBoundary = {9.5,4.7};
    std::vector<TH2F*> FieldHistograms;

    // Make Arguments
    unsigned NumberOfSteps = 200;
    double Frequency = 5.16816465e4; //5.1681646e4
    float PhaseShift = 0.0;

//   DriftVelocity = -2e4*DriftPar_0*pow(DriftField/100,DriftPar_1) /( 1+exp(DriftPar_2/pow(DriftField/100,DriftPar_3)) );

    unsigned NumberOfFiles = 11;
    std::string FileName = "FieldMaps/FieldMap_500_250_50_sinus_";

//   FieldHistograms = ReadDataFile(FileName,NumberOfFiles);
    FieldHistograms = SimpleFiled(NumberOfFiles);
    DrawHistos(FieldHistograms);

    std::cout << Simulation(FieldHistograms, ParticlePosition, MonitorBoundary, Frequency, PhaseShift, NumberOfFiles, NumberOfSteps) << std::endl;

    if(0)
    {
        double NewFrequency;
        std::vector<double> FrequencySweep = {1e4,1e5};
        std::vector<double> NumberOfPeriodes;

        for(unsigned range = 0; range < 2; range++) NumberOfPeriodes.push_back( Simulation(FieldHistograms, ParticlePosition, MonitorBoundary, FrequencySweep[range], PhaseShift, NumberOfFiles, NumberOfSteps) );

        for(unsigned sweep = 0; sweep < 10; sweep++)
        {
            NewFrequency = (FrequencySweep[1]+FrequencySweep[0])/2.0;
            double NewNoPeriods = Simulation(FieldHistograms, ParticlePosition, MonitorBoundary, NewFrequency, PhaseShift, NumberOfFiles, NumberOfSteps);

            std::cout << FrequencySweep[0] << " " << FrequencySweep[1] << " " << NewFrequency << " | " << NumberOfPeriodes[0] << " " << NumberOfPeriodes[1] << " " << NewNoPeriods << std::endl;

            if(NewNoPeriods > NumberOfPeriodes[0] || NewNoPeriods > NumberOfPeriodes[1])
            {
                if(NumberOfPeriodes[0] > NumberOfPeriodes[1])
                {
                    NumberOfPeriodes[1] = NewNoPeriods;
                    FrequencySweep[1] = NewFrequency;
                }
                else if(NumberOfPeriodes[0] != NumberOfPeriodes[1])
                {
                    NumberOfPeriodes[0] = NewNoPeriods;
                    FrequencySweep[0] = NewFrequency;
                }
                else
                {
                    std::cout << "Minimum found after " << sweep << " steps!" << std::endl;
                    break;
                }
            }
        }

        std::cout.precision(17);
        std::cout << "Ideal frequency is " << NewFrequency << std::endl;
    }

//   SectorList.resize(4*NumberOfFiles-2);


//   for(std::vector<unsigned>::iterator iter = SectorList.begin(); iter != SectorList.end(); iter++) std::cout << *iter << std::endl;

//   for(unsigned iter = 0; iter <10; iter++)
//   {
//     FieldHistograms.push_back(ReadDataFile(FileName+std::to_string(iter)+".txt"));
//   }

//   VerticalField -> SetMinimum(-8e4);

    return 0;
}

std::vector<TH2F*> ReadDataFile(std::string FileName, unsigned NumberOfFiles)
{
    bool r_BinCountFlag;

    std::string Line;
    std::string Cell;

    std::vector<std::string> HistoName = {"",""};

    std::vector<float> CoordMinMax;
    std::vector<float> CoordBinSize;
    std::vector<unsigned> NumberOfBins;

    std::vector<TH2F*> FieldHistograms;

    // Loop over files
    for(unsigned file = 0; file < NumberOfFiles; file++)
    {
        std::vector<std::vector<float>> RawFieldData;
        RawFieldData.resize(4); // Resize by the number of columns

        r_BinCountFlag = true;

        // Initialize file property vectors
        CoordMinMax = { (float)0xDEADBEEF , 0.0 , (float)0xDEADBEEF , 0.0 }; // Minimum and Maximum of both coordinates (r_min , r_max , z_min , z_max)
        CoordBinSize = { 0.0 , 0.0 }; // (delta_r , delta_z)
        NumberOfBins = { 0 , 0 }; // (#r_bins , #z_bins)

        std::ostringstream StringStreamFile;
        StringStreamFile << file;

        // Open field map file
        ifstream FieldFile (FileName + StringStreamFile.str() + ".txt");
        if(FieldFile.bad())
        {
            std::cout << "No such file or directory: " << (FileName + StringStreamFile.str() + ".txt") << std::endl;
            continue;
        }

        // Loop over lines until files end
        while(std::getline(FieldFile,Line))
        {
            if(Line[0] != '%')
            {
                // Loop over all columns
                for(unsigned column = 0; column < RawFieldData.size(); column++)
                {
                    // Only read data if data stream works
                    if(FieldFile >> Cell)
                    {
                        // Replace NaN with DEADBEEF
                        if(Cell != "NaN") RawFieldData[column].push_back(atof(Cell.c_str()));
                        else	      RawFieldData[column].push_back((float)0xDEADBEEF);

                        // Find minimum and maximum in data
                        if(column<2 && RawFieldData[column].back()<CoordMinMax[2*column]) 	     CoordMinMax[2*column] = RawFieldData[column].back();
                        else if(column<2 && RawFieldData[column].back()>CoordMinMax[2*column+1]) CoordMinMax[2*column+1] = RawFieldData[column].back();
                    }
                } // End of column loop

                // Rais number of bins for z if the r value repeats
                if(RawFieldData[0].back()==RawFieldData[0].front())
                {
                    NumberOfBins[1]++;
                    // Set r bin count flag false after the second repetition of the r value
                    if(r_BinCountFlag && RawFieldData[0].size()>1) r_BinCountFlag = false;
                }
                if(r_BinCountFlag) NumberOfBins[0]++;
            }
        } // End of line loop
        FieldFile.close();

        // Calculate bin size for both coordinates
        CoordBinSize[0] = (CoordMinMax[1]-CoordMinMax[0])/((float)NumberOfBins[0]-1.0);
        CoordBinSize[1] = (CoordMinMax[3]-CoordMinMax[2])/((float)NumberOfBins[1]-1.0);

        // Set histogram names
        HistoName = {"Radial_Field_"+StringStreamFile.str(),"Vertical_Field_"+StringStreamFile.str()};

        // Initialize and reset histograms for both coordinates
        FieldHistograms.push_back(new TH2F(HistoName[0].c_str(),HistoName[0].c_str(),NumberOfBins[0]+1,CoordMinMax[0]-CoordBinSize[0],CoordMinMax[1]+CoordBinSize[0],NumberOfBins[1],CoordMinMax[2],CoordMinMax[3]+CoordBinSize[1]));
        FieldHistograms.back() -> Reset();

        FieldHistograms.push_back(new TH2F(HistoName[1].c_str(),HistoName[1].c_str(),NumberOfBins[0]+1,CoordMinMax[0]-CoordBinSize[0],CoordMinMax[1]+CoordBinSize[0],NumberOfBins[1],CoordMinMax[2],CoordMinMax[3]+CoordBinSize[1]));
        FieldHistograms.back() -> Reset();

        // Loop to fill histograms
        for(unsigned z_bin = 0; z_bin < NumberOfBins[1]; z_bin++) for(unsigned r_bin = 0; r_bin < NumberOfBins[0]+1; r_bin++)
            {
                if(r_bin != 0) // Fill usual bins
                {
                    FieldHistograms[2*file] -> SetBinContent(r_bin+1,z_bin+1,RawFieldData[2][z_bin*NumberOfBins[0]+(r_bin-1)]);
                    FieldHistograms[2*file+1] -> SetBinContent(r_bin+1,z_bin+1,RawFieldData[3][z_bin*NumberOfBins[0]+(r_bin-1)]);
                }
                else // Fill first radial bin again with E_r[0] and E_z[0], in order to include r=0 in the histogram
                {
                    FieldHistograms[2*file] -> SetBinContent(r_bin+1,z_bin+1,RawFieldData[2][z_bin*NumberOfBins[0]+r_bin]);
                    FieldHistograms[2*file+1] -> SetBinContent(r_bin+1,z_bin+1,RawFieldData[3][z_bin*NumberOfBins[0]+r_bin]);
                }
            } // end histogram fill loop

        // Debug
//     for(unsigned iter = 0; iter < CoordMinMax.size(); iter++) std::cout << CoordMinMax[iter] << " ";
//     std::cout << CoordBinSize[0] << " " << CoordBinSize[1] << " " << NumberOfBins[0] << " " << NumberOfBins[1] << std::endl;

    } // End of loop over files

    // Return Histograms
    return FieldHistograms;
}

std::vector< TH2F* > SimpleFiled(unsigned NumberOfFiles)
{
    float E_0 = -70000;

    std::vector<TH2F*> FieldHistograms;

    std::vector<float> CoordMinMax = {0.0,10.0,-5.0,5.0};
    std::vector<float> CoordBinSize = {0.0,0.0};
    std::vector<unsigned> NumberOfBins = {100,100};

    std::vector<std::string> HistoName = {"",""};

    for(unsigned file = 0; file < NumberOfFiles; file++)
    {
        std::ostringstream StringStreamFile;
        StringStreamFile << file;

        HistoName = {"Radial_Field_"+StringStreamFile.str(),"Vertical_Field_"+StringStreamFile.str()};

        CoordBinSize[0] = (CoordMinMax[1]-CoordMinMax[0])/((float)NumberOfBins[0]-1.0);
        CoordBinSize[1] = (CoordMinMax[3]-CoordMinMax[2])/((float)NumberOfBins[1]-1.0);

        FieldHistograms.push_back(new TH2F(HistoName[0].c_str(),HistoName[0].c_str(),NumberOfBins[0]+1,CoordMinMax[0]-CoordBinSize[0],CoordMinMax[1]+CoordBinSize[0],NumberOfBins[1],CoordMinMax[2],CoordMinMax[3]+CoordBinSize[1]));
        FieldHistograms.back() -> Reset();

        FieldHistograms.push_back(new TH2F(HistoName[1].c_str(),HistoName[1].c_str(),NumberOfBins[0]+1,CoordMinMax[0]-CoordBinSize[0],CoordMinMax[1]+CoordBinSize[0],NumberOfBins[1],CoordMinMax[2],CoordMinMax[3]+CoordBinSize[1]));
        FieldHistograms.back() -> Reset();

        // Loop to fill histograms
        for(unsigned z_bin = 0; z_bin < NumberOfBins[1]; z_bin++) for(unsigned r_bin = 0; r_bin < NumberOfBins[0]+1; r_bin++)
            {
                float E_z;
                if(z_bin < (unsigned)(NumberOfBins[1]/2)) E_z = E_0 + z_bin*500;
                else E_z = 2.0*(E_0 + (unsigned)(NumberOfBins[1]/2*500)) - (z_bin-(unsigned)(NumberOfBins[1]/2))*500;

                FieldHistograms[2*file] -> SetBinContent(r_bin+1,z_bin+1,0.0);
                FieldHistograms[2*file+1] -> SetBinContent(r_bin+1,z_bin+1,E_z);
            } // end histogram fill loop

    }

    return FieldHistograms;
}

void DrawHistos(std::vector<TH2F*> FieldHistograms)
{
    std::string PicFormat = ".gif";
    std::string AnimationTime = "+20";
    std::vector<std::string> PicNames;

    PicNames.push_back("Radial_Field");
    PicNames.push_back("Vertical_Field");

    // Deleting files with the same name in this directory
    for(unsigned pics = 0; pics < PicNames.size(); pics++)
    {
        PicNames[pics] += PicFormat;

        if(std::remove( PicNames[pics].c_str() ) != 0)
            std::cout << "Creating new file: " << PicNames[pics] << std::endl;
        else
            std::cout << "Recreating file: " << PicNames[pics] << std::endl;

        if(PicFormat==".gif")
            PicNames[pics] += AnimationTime;
    }

    // Loop over all Histograms with even numbers
    for(unsigned entry = 0; entry < FieldHistograms.size(); entry+=2)
    {
        FieldHistograms[entry] -> SetStats(0);
        FieldHistograms[entry+1] -> SetStats(0);

        FieldHistograms[entry] -> SetMaximum(4e5);
//     FieldHistograms[entry+1] -> SetMaximum(4e5);


        FieldHistograms[entry] -> SetMinimum(-4e5);
//     FieldHistograms[entry+1] -> SetMinimum(-4e5);

        TCanvas * C1 = new TCanvas("E_r","E_r",800,500);
        FieldHistograms[entry] -> Draw("colz");
        C1 -> Print(PicNames[0].c_str(), PicFormat.c_str());

        TCanvas * C2 = new TCanvas("E_z","E_z",800,500);
        FieldHistograms[entry+1] -> Draw("colz");
        C2 -> Print(PicNames[1].c_str(), PicFormat.c_str());

        delete C1;
        delete C2;
    }
}

std::vector<float> VelocityFunction(std::vector<float> ElectricField)
{
    const float ChargeMobility = 475.; // cm^2/(V s)
    float ConversionFactor = 0.01; // Conversion from V/m to V/cm
    std::vector<float> DriftVelocity;

    for(int coord = 0; coord < 2; coord++) DriftVelocity.push_back(ChargeMobility*ConversionFactor*ElectricField[coord]);

    return DriftVelocity;
}

double Simulation(std::vector<TH2F*> FieldHistograms, std::vector<double> ParticlePosition, std::vector<float> MonitorBoundary, double Frequency, float PhaseShift, unsigned NumberOfFiles, unsigned NumberOfSteps)
{
    TApplication *App = new TApplication("TheApp",0,0);

    float InterpolationFactor;
    double Delta_t = 1.0/(Frequency*(double)NumberOfSteps);
    double Time = 0.0;
    unsigned NumberOfSectors = 4;

    std::vector<unsigned> SectorList;
    std::vector<float> DriftVelocity;
    std::vector<float> EField;

    EField.resize(3);

    std::vector<TCanvas*> PositionCanvas;
    PositionCanvas.push_back(new TCanvas("r(t)","r(t)",1000,500));
    PositionCanvas.push_back(new TCanvas("z(t)","z(t)",1000,500));

    std::vector<TGraph> PositionInTime;
    PositionInTime.resize(2);

    float FileStep;

    // Make a order of file calls in a sinus
    for(unsigned sector = 0; sector < NumberOfSectors; sector++) for(unsigned files = 0; files < NumberOfFiles-1; files++)
        {
            SectorList.push_back( files + (sector%2)*(NumberOfFiles-1-2*files) );
        }

//   SectorList = {10,10,10,10};

    unsigned GraphEntry = 0;

    while(ParticlePosition[0] < MonitorBoundary[0] && fabs(ParticlePosition[1]) < MonitorBoundary[1])
    {
        //Calculate Steps
        for(unsigned step = 0; step < NumberOfSteps; step++)
        {
            FileStep = (float)step*(float)SectorList.size()/(float)NumberOfSteps + (float)SectorList.size()*(float)PhaseShift/2.0;

            InterpolationFactor = FileStep - (int)FileStep;
            EField[2] = 0.0;

            if(ParticlePosition[0] < MonitorBoundary[0] && fabs(ParticlePosition[1]) < MonitorBoundary[1])
            {
                for(unsigned coord = 0; coord < 2; coord++)
                {
                    EField[coord] = 0.0;
                    if(fmod(FileStep,SectorList.size()) > (float)NumberOfSectors*((float)NumberOfFiles-1.0)/2.0)
                    {
                        EField[coord] -= pow(-1,coord)*(1.0-InterpolationFactor)*( FieldHistograms[2*SectorList[(int)FileStep%(SectorList.size())]+coord]->Interpolate(ParticlePosition[0],-ParticlePosition[1]) );
                        EField[coord] -= pow(-1,coord)*InterpolationFactor*( FieldHistograms[2*SectorList[((int)FileStep+1)%(SectorList.size())]+coord]->Interpolate(ParticlePosition[0],-ParticlePosition[1]) );
                    }
                    else
                    {
                        EField[coord] -= (1.0-InterpolationFactor)*( FieldHistograms[2*SectorList[(int)FileStep%(SectorList.size())]+coord]->Interpolate(ParticlePosition[0],ParticlePosition[1]) );
                        EField[coord] -= InterpolationFactor*( FieldHistograms[2*SectorList[((int)FileStep+1)%(SectorList.size())]+coord]->Interpolate(ParticlePosition[0],ParticlePosition[1]) );
                    }
                    EField[2] += pow(EField[coord],2);
                }
                EField[2] = sqrt(EField[2]);

                DriftVelocity = VelocityFunction(EField);

                for(unsigned coord = 0; coord < ParticlePosition.size(); coord++) ParticlePosition[coord] += DriftVelocity[coord]*Delta_t;

                for(unsigned coord = 0; coord < PositionInTime.size(); coord++)
                {
                    PositionInTime[coord].SetPoint(GraphEntry,Time,ParticlePosition[coord]);
                    PositionCanvas[coord] -> cd();
                    PositionInTime[coord].Draw("APL");
                    PositionCanvas[coord] -> Update();
                }
                Time += Delta_t;
                GraphEntry++;
            }
        }
    }

    //Run TApplication
//   App->Run();
    delete App;

    return Time*Frequency;
}

