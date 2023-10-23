#include <iostream>
#include "jakobian.h"

using namespace std;

// Funkcje ksztaltu

//  Wzgledem ksi

double pFksiN_1(double eta)
{
    return -0.25*(1-eta);
}

double pFksiN_2(double eta)
{
    return 0.25*(1-eta);
}

double pFksiN_3(double eta)
{
    return 0.25*(1+eta);
}

double pFksiN_4(double eta)
{
    return -0.25*(1+eta);
}

//  Wzgledem eta

double pFetaN_1(double ksi)
{
    return -0.25*(1-ksi);
}

double pFetaN_2(double ksi)
{
    return -0.25*(1+ksi);
}

double pFetaN_3(double ksi)
{
    return 0.25*(1+ksi);
}

double pFetaN_4(double ksi)
{
    return 0.25*(1-ksi);
}

// Jakobian przeksztalcenia

double JakobianEta(double A, double B, double C, double D, double eta)
{
    return pFksiN_1(eta)*A + pFksiN_2(eta)*B + pFksiN_3(eta)*C + pFksiN_4(eta)*D;
}

double JakobianKsi(double A, double B, double C, double D, double ksi)
{
    return pFetaN_1(ksi)*A + pFetaN_2(ksi)*B + pFetaN_3(ksi)*C + pFetaN_4(ksi)*D;
}

// Dla macierzy Hbc

double funKsz1(double ksi, double eta)
{
    return 0.25*(1.0-ksi)*(1.0-eta);
}

double funKsz2(double ksi, double eta)
{
    return 0.25*(1.0+ksi)*(1.0-eta);
}

double funKsz3(double ksi, double eta)
{
    return 0.25*(1.0+ksi)*(1.0+eta);
}

double funKsz4(double ksi, double eta)
{
    return 0.25*(1.0-ksi)*(1.0+eta);
}
/*
double funKsz_w(double ksi, double eta, int w)
{
    switch(w)
    {
    // Dol ( N1 i N2 )
    case 0:
        return 0.25*(1.0-ksi)*(1.0-eta);
        break;
    // Lewo ( N2 i N3 )
    case 1:
        return 0.25*(1.0+ksi)*(1.0-eta);
        break;
    // Gora ( N3 i N4 )
    case 2:
        return 0.25*(1.0+ksi)*(1.0+eta);
        break;
    // Prawo ( N4 i N1 )
    case 3:
        return 0.25*(1.0-ksi)*(1.0+eta);
        break;
    // Blad
    default:
        return -1.0;
        break;
    }
}
*/
double funKsz_w(double ksi, double eta, int w)
{
    switch(w)
    {
    // Dol ( N1 i N2 )
    case 0:
        return 0.25*(1.0+ksi)*(1.0-eta);
        break;
    // Prawo ( N2 i N3 )
    case 1:
        return 0.25*(1.0+ksi)*(1.0+eta);
        break;
    // Gora ( N3 i N4 )
    case 2:
        return 0.25*(1.0-ksi)*(1.0+eta);
        break;
    // Lewo ( N4 i N1 )
    case 3:
        return 0.25*(1.0-ksi)*(1.0-eta);
        break;
    // Blad
    default:
        return -1.0;
        break;
    }
}









