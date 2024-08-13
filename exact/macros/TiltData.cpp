R__LOAD_LIBRARY(libExACT.so)
#include <TH1.h>
#include <TGraph.h>
#include <vector>
#include <fstream>

class Tilting
{
public:
    Tilting(const std::string &filename);
    ~Tilting();
    std::vector<double> times;
    std::vector<double> tilts;
};

Tilting::Tilting(const std::string &filename)
{
    cout << "Loading file: " << filename << endl;

    std::ifstream file(filename);
    std::string line;

    // Skip the header line
    std::getline(file, line);

    // Read data lines and populate the vectors
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double time, tilt;
        iss >> time >> tilt;
        times.push_back(time);
        tilts.push_back(tilt);
    }
};

Tilting::~Tilting(){};

void TiltData()
{
    // Declare Variables
    const std::string FileName = "../macros/Tilt_Sensor_Data.csv";
    const int offsetrange = 938;
    double offset[offsetrange] = {0};
    double corrected_tilt[offsetrange];

    // Load tilt data
    Tilting tilt(FileName);

    // Convert tilt into easier to handle arrays
    double time[tilt.times.size()];
    copy(tilt.times.begin(), tilt.times.end(), time);

    double til[tilt.tilts.size()];
    copy(tilt.tilts.begin(), tilt.tilts.end(), til);

    // Apply corrections to the offset until it forms a smooth curve
    for (int i = 0; i < offsetrange; i++)
    {
        if (i < 15)
        {
            offset[i] = 0 - til[i] + til[100] - 0.2;
        }
        else if (15 <= i && i < 233)
        {
            offset[i] = -til[i] + til[100];
        }
        else if (233 <= i && i < 246)
        {
            offset[i] = +0.7 - til[i] + til[100] - 0.4;
        }
        else if (246 <= i && i < 340)
        {
            offset[i] = -3.5 - til[i] + til[100];
        }
        else if (340 <= i && i < 345)
        {
            offset[i] = -2.05 - til[i] + til[100] - 0.4;
        }
        else if (345 <= i && i < 358)
        {
            offset[i] = -3.0 - til[i] + til[100] - 0.4;
        }
        else if (358 <= i && i < 383)
        {
            offset[i] = +4.0 - til[i] + til[100] - 0.4;
        }
        else if (382 <= i && i < 626)
        {
            offset[i] = -2.8 - til[i] + til[100];
        }
        else if (626 <= i && i < 668)
        {
            offset[i] = +0.3 - til[i] + til[100];
        }
        else if (668 <= i && i < 676)
        {
            offset[i] = +3.3 - til[i] + til[100];
        }
        else if (676 <= i && i < 725)
        {
            offset[i] = -2.5 - til[i] + til[100];
        }
        else if (725 <= i && i < 759)
        {
            offset[i] = 4.4 - til[i] + til[100];
        }
        else if (759 <= i && i < 950)
        {
            offset[i] = -1.5 - til[i] + til[100];
        }
        else
        {
            offset[i] = -til[i] + til[100];
        }
    }

    // Apply correction to the tilt data
    for (int i = 0; i < offsetrange; i++)
    {
        corrected_tilt[i] = til[i] + offset[i];
    }

    // Save results to file
    ofstream corrTiltFile;
    corrTiltFile.open("corrected_tilt.csv");
    corrTiltFile << "Time in sec"
                 << ", "
                 << "Raw Tilt in deg"
                 << ", "
                 << "Offset in deg"
                 << ", "
                 << "Corrected tilt in deg" << endl;
    ;
    for (int i = 0; i < offsetrange; i++)
    {
        corrTiltFile << time[i] << ", " << til[i] << ", " << offset[i] << ", " << corrected_tilt[i] << endl;
    }
    corrTiltFile.close();

    // Create a plot of the results
    TCanvas *c1 = new TCanvas("c1", "Exclusion graphs examples", 200, 10, 600, 400);
    c1->SetGrid();

    TMultiGraph *mg = new TMultiGraph();
    TGraph *g1 = new TGraph(offsetrange, time, til);
    TGraph *g2 = new TGraph(offsetrange, time, offset);
    TGraph *g3 = new TGraph(offsetrange, time, corrected_tilt);
    g2->SetLineColor(4);
    g3->SetLineColor(2);

    mg->Add(g1);
    mg->Add(g2);
    mg->Add(g3);
    mg->Draw("AL");

    mg->GetHistogram()->GetYaxis()->SetTitle("Angle in deg");
    mg->GetHistogram()->GetXaxis()->SetTitle("Time in sec since start of the night on the 13th");

    auto legend = new TLegend(0.1, 0.75, 0.3, 0.9);
    legend->AddEntry(g1, "Original data", "l");
    legend->AddEntry(g2, "Correction factor", "l");
    legend->AddEntry(g3, "Corrected data", "l");
    legend->Draw();
}