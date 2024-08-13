#include<iostream>


class IPlotTools{
	public:
		IPlotTools();
		~IPlotTools();
		static void FindBin(int pixelID, int *nx, int *ny);
		static int FindPixel(int nx, int ny);
		static int GetMUSICCursorID(int x, int y);
		static void DrawMUSICBoundaries(); // Function declaration added here

};
