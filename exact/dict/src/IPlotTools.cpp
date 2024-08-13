#include "IPlotTools.h"
#include "TBox.h"

IPlotTools::IPlotTools(){

}
IPlotTools::~IPlotTools(){
    
}
void IPlotTools::FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8; 
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}
int IPlotTools::FindPixel(int nx, int ny)
{
    int pixelID;
    pixelID = (nx*4+(ny%4)+(ny/4)*128);
    return pixelID;
}
int IPlotTools::GetMUSICCursorID(int x, int y)
{
    int MUSIC_column = x/2;
    int MUSIC_row = y/4;
    int MUSIC_ID = MUSIC_column+MUSIC_row*16;

    return MUSIC_ID;  
}
void IPlotTools::DrawMUSICBoundaries() {
    TBox *b = new TBox(-0.5, -0.5, 1.5, 3.5);
    b->SetFillStyle(0);
    b->SetLineColor(kRed);
    b->Draw();

    for (int i = 0; i < 512 / 8; i++) {
        TBox *bn = (TBox *)b->Clone();
        bn->SetX1((i % 16) * 2 - 0.5);
        bn->SetX2((i % 16) * 2 + 1.5);
        bn->SetY1((i / 16) * 4 - 0.5);
        bn->SetY2((i / 16) * 4 + 3.5);
        bn->Draw();
    }
}
