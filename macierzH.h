#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

class _matrixH
{
    double ** macH;
    double ** macC;
    double * vecP;
    double * vecTemp;

    double ** macRoz;

    double ** niestH; // Niestacionarna H

    double ** pomC;
    double ** pomRoz;   // Pomocnicza macierz rozszrzona
    double * pomP;      // Pomocniczy wektor P
    double * pomTemp;   // pomocniczy wektor temperatur

    int rozmH;

public:
    _matrixH(int noN);

    void wypisz();

    void wypiszP();

    void wypiszTemp();

    void wypiszC();

    void obliczH(double dane[2][4], int nodeId[4], double w1, double w2, int con, int dens, int speh);

    void obliczHbc(double dane[2][4], int nodeId[4], int nodeBC[4], double w1, double w2, int alf);

    void obliczP(double dane[2][4], int nodeId[4], int nodeBC[4], double wag1, double wag2, int alf, int tot);

    int elimGauss();

    void stanNiestacionarny(int krok, int initemp);

    int krokSymulacji(int, double*);

    void obliczH_3P(double dane[2][4], int nodeId[4], double w1, double w2, int con, int dens, int speh);
};
