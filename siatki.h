#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

struct _Node
{
    int id;
    double x, y, t;
    bool bc;
};

struct _Element
{
    int ID[4];
};

class _Grid
{
    int varTab[10];
    // Zawartosc tablicy varTab:
    // 0 - SimulationTime
    // 1 - SimulationStepTime
    // 2 - Conductivity
    // 3 - Alfa
    // 4 - Tot
    // 5 - InitialTemp
    // 6 - Density
    // 7 - SpecificHeat
    // 8 - number of Nodes
    // 9 - number of Elements
    _Node* listaN;
    _Element* listaE;

public:
    _Grid(string);

    void wypisz();

    int getVar(int n)
    {
        return varTab[n];
    }

    double getNx(int i)
    {
        _Node* wskPom = &(listaN[i]);
        return wskPom->x;
    }

    double getNy(int i)
    {
        _Node* wskPom = &(listaN[i]);
        return wskPom->y;
    }

    int getEl(int n, int i)
    {
        _Element* wskEl = &(listaE[n]);
        return wskEl->ID[i];
    }

    int getNbc(int i)
    {
        _Node* wskPom = &(listaN[i]);
        if(wskPom->bc)
            return 1;
        else
            return 0;
    }
};
