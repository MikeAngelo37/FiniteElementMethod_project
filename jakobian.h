#include <iostream>

using namespace std;

// Funkcje ksztaltu

//  Wzgledem ksi

double pFksiN_1(double eta);
double pFksiN_2(double eta);
double pFksiN_3(double eta);
double pFksiN_4(double eta);

//  Wzgledem eta

double pFetaN_1(double ksi);
double pFetaN_2(double ksi);
double pFetaN_3(double ksi);
double pFetaN_4(double ksi);

// Jakobian przeksztalcenia

double JakobianEta(double A, double B, double C, double D, double eta);

double JakobianKsi(double A, double B, double C, double D, double ksi);

double funKsz1(double ksi, double eta);
double funKsz2(double ksi, double eta);
double funKsz3(double ksi, double eta);
double funKsz4(double ksi, double eta);

double funKsz_w(double ksi, double eta, int w);
