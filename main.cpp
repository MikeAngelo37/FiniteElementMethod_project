#include <iostream>
#include <math.h>
#include "siatki.h"
#include "jakobian.h"
#include "macierzH.h"

using namespace std;

int main()
{
// Ladowane danych do programu
    _Grid Tablica("Test1_4_4.txt");

// Wypisywanie zaladowanych danych
    Tablica.wypisz();

// Tworzenie macierzy H
    int n = Tablica.getVar(9);
    int m = Tablica.getVar(8);
    _matrixH macH(m);
    //macH.wypisz();  // Debug

// Obliczanie czastkowych H
    for(int i=0;i<n;i++)
    {
        double element[2][4];
        int idNode[4];
        int cond = Tablica.getVar(2);
        int dence = Tablica.getVar(6);
        int specheat = Tablica.getVar(7);

        for(int j=0;j<4;j++)
        {
            idNode[j] = Tablica.getEl(i, j);
            element[0][j] = Tablica.getNx(idNode[j]-1);
            element[1][j] = Tablica.getNy(idNode[j]-1);

            //cout << "Node " << idNode[j] << endl;
            //cout << " " << element[0][j] << " | " << element[1][j] << endl;
        }
        cout << "Obliczanie dla elementu " << i+1 << "..." << endl;
        macH.obliczH(element, idNode, 1.0, 1.0, cond, dence, specheat);
    }

// Wypisywanie macierzy H
    macH.wypisz();


// Obliczanie czastkowych Hbc i wektora P
    for(int i=0;i<n;i++)
    {
        double element[2][4];
        int idNode[4];
        int nodeBc[4];
        int alfa = Tablica.getVar(3);
        int tot = Tablica.getVar(4);

        for(int j=0;j<4;j++)
        {
            idNode[j] = Tablica.getEl(i, j);
            element[0][j] = Tablica.getNx(idNode[j]-1);
            element[1][j] = Tablica.getNy(idNode[j]-1);
            nodeBc[j] = Tablica.getNbc(idNode[j]-1);

            //cout << "Node " << idNode[j];
            //if(nodeBc[j])
            //    cout << "\t(BC)";
            //cout << endl;
            //cout << " " << element[0][j] << " | " << element[1][j] << endl;
        }
        cout << "Obliczanie dla elementu " << i+1 << "..." << endl;
        macH.obliczHbc(element, idNode, nodeBc, 1.0, 1.0, alfa);
        macH.obliczP(element, idNode, nodeBc, 1.0, 1.0, alfa, tot);
    }

// Wypisywanie macierzy H + Hbc
    macH.wypisz();

// Wypisywanie wektora P
    macH.wypiszP();

// [H] = [H] + [C]/dT
// {P} = {P}+{[C]/dT}*{T0}
// Wyliczyæ to przed wylizaniem temperature

// Obliczanie i wypisywanie wektora temperatur
    if(macH.elimGauss()==-1)
        cout << "Error!" << endl;
    else
        macH.wypiszTemp();

// Wypisywaie C
    macH.wypiszC();

// Obliczanie stau niestacionarnego
    cout << "Obliczanie wektora T dla kazdego kroku czasowego" << endl;
    int krokCzasowy = Tablica.getVar(1);
    int tempStartowa = Tablica.getVar(5);
    //macH.stanNiestacionarny(krokCzasowy, tempStartowa);
    n = Tablica.getVar(8);
    double * vecTemp = new double[n];
    for(int i=0;i<n;i++)
    {
        vecTemp[i] = Tablica.getVar(5);
    }

    cout << "--Wektor temperatur - wartosci poczatkowe--" << endl;
    //cout << endl;

    for(int i=0;i<n;i++)
    {
        cout << " | " << vecTemp[i];
    }
    cout << endl;
    //cout << "------------" << endl;
    //cout << endl;


    int czasSym = Tablica.getVar(0);
    for(int i=krokCzasowy;i<=czasSym;i=i+krokCzasowy)
    {
        //cout << "--------------" << endl;
        //cout << "T = " << i << endl;

        macH.krokSymulacji(krokCzasowy, vecTemp);

        //cout << "--Wektor temperatur--" << endl;
        for(int i=0;i<n;i++)
        {
            //cout << vecTemp[i] << " | ";
        }
        //cout << endl;
        cout << "------------" << endl;
        // Najmniejsza i najwieksza temperatura
        double TempMin = Tablica.getVar(4) * 1.0;   // Temperatura otoczenia (najwieksza mozliwa)
        double TempMax = Tablica.getVar(5) * 1.0;   // Temperatura obiektu (najmniejsza mozliwa)
        for(int i=0;i<n;i++)
        {
            if(TempMin>vecTemp[i])
                TempMin=vecTemp[i];

            if(TempMax<vecTemp[i])
                TempMax=vecTemp[i];
        }
        cout << " T = " << i << " | Temp. Min. =  " << TempMin << " | Temp. Max. =  " << TempMax <<  endl;
        //cout << endl;
        //cout << " ------------" << endl;
    }

    return 0;
}
