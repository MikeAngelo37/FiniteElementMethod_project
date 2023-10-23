#include <iostream>
#include <math.h>
#include <vector>
#include "macierzH.h"
#include "jakobian.h"

using namespace std;

const double eps = 1e-12; // sta³a przyblizenia zera

_matrixH::_matrixH(int noN) // number of Nodes
{
    // Macierz H + Hbc
    macH = new double *[noN];
    for(int i=0;i<noN;i++)
    {
        macH[i] = new double[noN];
        for(int j=0;j<noN;j++)
        {
            macH[i][j] = 0.0;
        }
    }

    // Macierz C
    macC = new double *[noN];
    for(int i=0;i<noN;i++)
    {
        macC[i] = new double[noN];
        for(int j=0;j<noN;j++)
        {
            macC[i][j] = 0.0;
        }
    }

    // Wektor P i wektor Temperatur
    vecP = new double [noN];
    vecTemp = new double [noN];
    for(int i=0;i<noN;i++)
    {
        vecP[i] = 0.0;
        vecTemp[i] = 0.0;
    }

    // Macierz roszerzona
    macRoz = new double *[noN];
    for(int i=0;i<noN;i++)
    {
        macRoz[i] = new double[noN+1];
        for(int j=0;j<(noN+1);j++)
        {
            macRoz[i][j] = 0.0;
        }
    }

    rozmH = noN;
}

void _matrixH::wypisz()
{
    cout << "---Macierz_H---" << endl;
    for(int i=0;i<rozmH;i++)
    {
        for(int j=0;j<rozmH;j++)
        {
            cout << " ";
            cout.precision(6);
            cout.width(8);
            cout.fill();
            cout << macH[i][j] << " |";
            //cout << "\t|";
        }
        cout << endl << endl;
    }
    cout << "---------------" << endl;
}


void _matrixH::wypiszC()
{
    cout << "---Macierz_C---" << endl;
    for(int i=0;i<rozmH;i++)
    {
        for(int j=0;j<rozmH;j++)
        {
            cout << " ";
            cout.precision(6);
            cout.width(8);
            cout.fill();
            cout << macC[i][j] << " |";
            //cout << "\t|";
        }
        cout << endl << endl;
    }
    cout << endl;
    cout << "---------------" << endl;
}


void _matrixH::wypiszP()
{
    cout << "--Wektor_P--" << endl;
    for(int i=0;i<rozmH;i++)
    {
        cout << vecP[i] << " | ";
    }
    cout << endl;
    cout << "------------" << endl;
    cout << endl;
}


void _matrixH::wypiszTemp()
{
    cout << "--Wektor_T--" << endl;
    for(int i=0;i<rozmH;i++)
    {
        cout << vecTemp[i] << " | ";
    }
    cout << endl;
    cout << "------------" << endl;
    cout << endl;
}


void _matrixH::obliczH(double dane[2][4], int nodeId[4], double wag1, double wag2, int con, int dens, int speh)
{
    double conductivity = 1.0 * con;
    double density = 1.0 * dens;
    double specheat = 1.0 * speh;
    //double conductivity = 25.0;

    double nody[2][4]; // 0 - X, 1 - Y
    for(int i=0;i<4;i++)
    {
        nody[0][i] = dane[0][i]; // X
        nody[1][i] = dane[1][i]; // Y
    }
/*
    cout << "id | X | Y" << endl;
    for(int i=0;i<4;i++)
    {
        cout << nodeId[i] << " | " << nody[0][i] << " | " << nody[1][i] << endl;
    }
*/
    //cout << "Macierz H" << endl;

    double ksi[4] = {-1.0, 1.0, -1.0, 1.0};
    double eta[4] = {-1.0, -1.0, 1.0, 1.0};
    for(int i=0;i<4;i++)
    {
        ksi[i] = ksi[i] * 1.0/sqrt(3);
        eta[i] = eta[i] * 1.0/sqrt(3);
    }

    //cout << " ksi | eta" << endl;
    //for(int i=0;i<4;i++)
    //{
    //    cout << " " << ksi[i] << " | " << eta[i] << endl;
    //}

    double poKsi[4][4];
    double poEta[4][4];

    for(int i=0;i<4;i++)
    {
        poKsi[i][0] = pFksiN_1(eta[i]);
        poKsi[i][1] = pFksiN_2(eta[i]);
        poKsi[i][2] = pFksiN_3(eta[i]);
        poKsi[i][3] = pFksiN_4(eta[i]);
    }

    for(int i=0;i<4;i++)
    {
        poEta[i][0] = pFetaN_1(ksi[i]);
        poEta[i][1] = pFetaN_2(ksi[i]);
        poEta[i][2] = pFetaN_3(ksi[i]);
        poEta[i][3] = pFetaN_4(ksi[i]);
    }

    //for(int i=0;i<4;i++)
    //{
    //    for(int j=0;j<4;j++)
    //    {
    //        cout << " " << poEta[i][j] << " |";
    //    }
    //    cout << endl;
    //}

    //double nody[2][4] = {// 0 - X, 1 - Y
    //    {0.0, 0.025, 0.025, 0.0},
    //    {0.0, 0.0, 0.025, 0.025}
    //};

    double mJakob[2][2] = {0.0};

    // dX / dKsi
    mJakob[0][0] = (nody[0][0]*poKsi[0][0]) + (nody[0][1]*poKsi[0][1]) + (nody[0][2]*poKsi[0][2]) + (nody[0][3]*poKsi[0][3]);
    // dY / dKsi
    mJakob[0][1] = (nody[1][0]*poKsi[0][0]) + (nody[1][1]*poKsi[0][1]) + (nody[1][2]*poKsi[0][2]) + (nody[1][3]*poKsi[0][3]);
    // dX / dEta
    mJakob[1][0] = (nody[0][0]*poEta[0][0]) + (nody[0][1]*poEta[0][1]) + (nody[0][2]*poEta[0][2]) + (nody[0][3]*poEta[0][3]);
    // dY / dEta
    mJakob[1][1] = (nody[1][0]*poEta[0][0]) + (nody[1][1]*poEta[0][1]) + (nody[1][2]*poEta[0][2]) + (nody[1][3]*poEta[0][3]);

    //cout << mJakob[0][0] << " | " << mJakob[0][1] << endl;
    //cout << mJakob[1][0] << " | " << mJakob[1][1] << endl;

    double detJakob = (mJakob[0][0]*mJakob[1][1]) - (mJakob[0][1]*mJakob[1][0]);
    //cout << "detJ = " << detJakob << endl;

    double odwJakob[2][2];
    odwJakob[0][0] =  mJakob[1][1] * (1/detJakob);
    odwJakob[0][1] = -mJakob[0][1] * (1/detJakob);
    odwJakob[1][0] = -mJakob[1][0] * (1/detJakob);
    odwJakob[1][1] =  mJakob[0][0] * (1/detJakob);

    double pcX[4][4]; // pochodne funkcji ksztaltu po X
    double pcY[4][4]; // pochodne funkcji ksztaltu po Y

    for(int i=0;i<4;i++) // pcX
    {
        for(int j=0;j<4;j++)
        {
            pcX[i][j] = odwJakob[0][0]*poKsi[i][j] + odwJakob[0][1]*poEta[i][j];
        }
    }

    for(int i=0;i<4;i++) // pcY
    {
        for(int j=0;j<4;j++)
        {
            pcY[i][j] = odwJakob[1][0]*poKsi[i][j] + odwJakob[1][1]*poEta[i][j];
        }
    }
/*
    cout << "pcX" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << pcX[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    cout << "pcY" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << pcY[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;
*/
    double localH[4][4] = {0.0};
    double localC[4][4] = {0.0};
    double tempH[4][4] = {0.0};
    double tempC[4][4] = {0.0};
    double tempCrazyC[4][4] = {0.0};
    //double w1 = 1.0, w2 = 1.0;

    double ksiC[4] = {-1.0, -1.0, 1.0, 1.0};
    double etaC[4] = {-1.0, 1.0, -1.0, 1.0};
    for(int i=0;i<4;i++)
    {
        ksiC[i] = ksiC[i] * 1.0/sqrt(3);
        etaC[i] = etaC[i] * 1.0/sqrt(3);
    }

    for(int i=0;i<4;i++)
    {
        tempC[i][0] = funKsz1(ksiC[i], etaC[i]);
        tempC[i][1] = funKsz2(ksiC[i], etaC[i]);
        tempC[i][2] = funKsz3(ksiC[i], etaC[i]);
        tempC[i][3] = funKsz4(ksiC[i], etaC[i]);
    }

    //cout << " {C} Funkcji ksztaltu: " << endl;
    //for(int i=0;i<4;i++)
    //{
    //    cout << tempC[i][0] << " | " << tempC[i][1] << " | " << tempC[i][2] << " | " << tempC[i][3] << endl;
    //}

    //for(int i=0;i<4;i++)
    //{
    //    for(int j=0;j<4;j++)
    //    {
    //        localC[i][j] = tempC[i][j] * tempC[j][i] * detJakob * density * specheat;
    //    }
    //}

    double pcXrazyT[4][4]; // iloczyn transponowanej i pcX
    double pcYrazyT[4][4]; // iloczyn transponowanej i pcY

    for(int n=0;n<4;n++)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {   // Mno¿enie macierzy przez nie same transponowane
                pcXrazyT[i][j] = pcX[n][i]*pcX[n][j];
                pcYrazyT[i][j] = pcY[n][i]*pcY[n][j];
                tempCrazyC[i][j] = tempC[n][i] * tempC[n][j];
            }
        }


        for(int i=0;i<4;i++)
        {

            for(int j=0;j<4;j++)
            {
                tempH[i][j] = (pcXrazyT[i][j] + pcYrazyT[i][j]) * conductivity * detJakob;
                localH[i][j] = localH[i][j] + tempH[i][j];
                //localH[i][j] = localH[i][j] * wag1 * wag2;

                // Macierz C
                //tempCrazyC[i][j] = tempC[i][j] * tempC[j][i]; // Razy transponowana
                tempCrazyC[i][j] = tempCrazyC[i][j] * detJakob * density * specheat;
                localC[i][j] = localC[i][j] + tempCrazyC[i][j];
                //tempC[i][j] = (pcXrazyT[i][j] + pcYrazyT[i][j]) * conductivity * detJakob;
                //localC[i][j] = localC[i][j] + tempC[i][j];
            }
        }

        //for(int i=0;i<4;i++)
        //{
        //    for(int j=0;j<4;j++)
        //    {
        //         cout << " " << tempH[i][j] << " |";
        //    }
        //    cout << endl;
        //}
        //cout << endl;
    }
/*
    cout << "DEBUG: Lokalna H" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << localH[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    cout << "DEBUG: Lokalna C" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << localC[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;
*/
    // Agregacja
    int Agr_x[4][4];
    int Agr_y[4][4];

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            Agr_x[i][j] = nodeId[i];
            Agr_y[i][j] = nodeId[j];
        }
    }

    //cout << "Macierz agregacji:" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            int pomX = Agr_x[i][j];// - 1;
            int pomY = Agr_y[i][j];// - 1;
            //cout << " " << Agr_x[i][j] << "," << Agr_y[i][j] << " |";

            macH[pomX-1][pomY-1] =  macH[pomX-1][pomY-1] + localH[i][j];    // Macierz H
            macC[pomX-1][pomY-1] =  macC[pomX-1][pomY-1] + localC[i][j];    // Macierz C
        }
        //cout << endl;
    }
    //cout << endl << endl;
}


void _matrixH::obliczHbc(double dane[2][4], int nodeId[4], int nodeBC[4], double wag1, double wag2, int alf)
{
    double alpha = 1.0 * alf;
    double sumaHbc[4][4] = {0.0};

    double nody[2][4]; // 0 - X, 1 - Y

    for(int i=0;i<4;i++)
    {
        nody[0][i] = dane[0][i]; // X
        nody[1][i] = dane[1][i]; // Y
    }
/*
    cout << "id | X | Y" << endl;
    for(int i=0;i<4;i++)
    {
        cout << nodeId[i];
        if(nodeBC[i])
            cout << " | BC";
        else
            cout << " |   ";

        cout << " | " << nody[0][i] << " | " << nody[1][i] << endl;
    }

    double ksi1[4] = {-1.0/sqrt(3), 1.0, -1.0/sqrt(3), -1.0};
    double ksi2[4] = { 1.0/sqrt(3), 1.0,  1.0/sqrt(3), -1.0};

    double eta1[4] = {-1.0, -1.0/sqrt(3),  1.0, -1.0/sqrt(3)};
    double eta2[4] = {-1.0,  1.0/sqrt(3),  1.0,  1.0/sqrt(3)};
*/
    double ksi1[4] = {-1.0/sqrt(3), 1.0, 1.0/sqrt(3), -1.0};
    double ksi2[4] = { 1.0/sqrt(3), 1.0,  -1.0/sqrt(3), -1.0};

    double eta1[4] = {-1.0, -1.0/sqrt(3),  1.0, 1.0/sqrt(3)};
    double eta2[4] = {-1.0,  1.0/sqrt(3),  1.0,  -1.0/sqrt(3)};

    //cout << "Obliczanie lokalnej Hbc:" << endl;
    for(int i=0;i<4;i++)
    {
        int pom1 = i;
        int pom2 = (i+1)%4; // Zabespieczenie, by indeksy
                            // nie wyszly poza tablice
        // Liczenie boku
        //if(1) // DEBUG
        if(nodeBC[pom1] && nodeBC[pom2])
        {
            //cout << "Bok " << 1+i << ": " << nodeId[pom1] << " | " << nodeId[pom2] << endl;

            double pomX1, pomX2, pomY1, pomY2;

            pomX1 = nody[0][pom1];
            pomX2 = nody[0][pom2];

            pomY1 = nody[1][pom1];
            pomY2 = nody[1][pom2];

            // Jakobian
            double detJAB = (sqrt(pow(pomX2-pomX1, 2)+pow(pomY2-pomY1, 2)))/2.0;

            double H_pc1[4][4] = { 0 };
            double H_pc2[4][4] = { 0 };
            double H_bc[4][4] = { 0 };

            //ksi[i*2]        eta[i*2]      // Pierwszy punkt
            //ksi[(i*2)+1]    eta[(i*2)+1]  // Drugi punktt

            double dane_tab[2][4] = { 0 };

            for(int j=0;j<4;j++)
            {
                double pomKsi, pomEta;
                //dane_tab[0][4] // tu ma byc petla po czterech
                //dane_tab[1][4] // funkcjach ksztaltu

                pomKsi = ksi1[j];
                pomEta = eta1[j];
                dane_tab[0][j] = funKsz_w(pomKsi, pomEta, i);

                pomKsi = ksi2[j];
                pomEta = eta2[j];
                dane_tab[1][j] = funKsz_w(pomKsi, pomEta, i);
            }
/*
            cout << "Dane: P1: "<< nodeId[0] <<" | P2: "<< nodeId[1] <<" | P3: "<< nodeId[2] <<" | P4: "<< nodeId[3] <<" |" << endl;
            for(int i=0;i<2;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    cout << " " << dane_tab[i][j] << "\t|";
                }
                cout << endl;
            }
            cout << endl;
*/
        // Mnozenie macierzy pc1
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    H_pc1[i][j] = dane_tab[0][i] * dane_tab[0][j];
                    //cout << "| " << dane_tab[0][3+i] << " * " << dane_tab[0][3+j] << "\t";
                    H_pc1[i][j] = H_pc1[i][j] * wag1 * alpha;
                }
                //cout << "|" << endl;
            }
/*
            cout << "Macierz H_pc1:" << endl;
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    cout << "| " << H_pc1[i][j] << "\t";
                }
                cout << "|" << endl;
            }
            cout << endl;
*/
        // Mnozenie macierzy pc2
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    H_pc2[i][j] = dane_tab[1][i] * dane_tab[1][j];
                    //cout << "| " << dane_tab[0][3+i] << " * " << dane_tab[0][3+j] << "\t";
                    H_pc2[i][j] = H_pc2[i][j] * wag2 * alpha;
                }
                //cout << "|" << endl;
            }
            //cout << endl;
/*
            cout << "Macierz H_pc2:" << endl;
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    cout << "| " << H_pc2[i][j] << "\t";
                }
                cout << "|" << endl;
            }
            cout << endl;
*/
        // Obliczanie macierzy H_bc
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    //H_pc2[i][j] = dane_tab[1][3+i] * dane_tab[1][3+j];
                    //cout << "| " << dane_tab[0][3+i] << " * " << dane_tab[0][3+j] << "\t";
                    //H_bc[i][j] = H_pc2[i][j] * w1 * alpha;
                    H_bc[i][j] = (H_pc1[i][j] + H_pc2[i][j]) * detJAB; // Jakobian = 0.0125
                    sumaHbc[i][j] = sumaHbc[i][j] + H_bc[i][j];
                }
                //cout << "|" << endl;
            }
            //cout << endl;
/*
            cout << "Macierz H_bc:" << endl;
            for(int i=0;i<4;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    cout << "| " << H_bc[i][j] << "\t";
                }
                cout << "|" << endl;
            }
            cout << endl;
*/
        }
/*
        cout << "sumaHbc:" << endl;
        for(int i=0;i<4;i++) // i - wiersz
        {
            for(int j=0;j<4;j++) // j - kolumna
            {
                cout << "| " << sumaHbc[i][j] << "\t";
            }
            cout << "|" << endl;
        }
        cout << endl;
*/
    }
    // Agregacja
    int Agr_x[4][4];
    int Agr_y[4][4];

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            Agr_x[i][j] = nodeId[i];
            Agr_y[i][j] = nodeId[j];
        }
    }

    //cout << "Macierz agregacji:" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            int pomX = Agr_x[i][j];// - 1;
            int pomY = Agr_y[i][j];// - 1;
            //cout << " " << Agr_x[i][j] << "," << Agr_y[i][j] << " |";

            macH[pomX-1][pomY-1] =  macH[pomX-1][pomY-1] + sumaHbc[i][j];
        }
        //cout << endl;
    }
    //cout << endl << endl;

}


void _matrixH::obliczP(double dane[2][4], int nodeId[4], int nodeBC[4], double wag1, double wag2, int alf, int tot)
{
    double alpha = 1.0 * alf;
    double temperOt = 1.0 * tot;

    double vectorP[4] = {0.0};

    double nody[2][4]; // 0 - X, 1 - Y

    for(int i=0;i<4;i++)
    {
        nody[0][i] = dane[0][i]; // X
        nody[1][i] = dane[1][i]; // Y
    }
/*
    cout << "id | X | Y" << endl;
    for(int i=0;i<4;i++)
    {
        cout << nodeId[i];
        if(nodeBC[i])
            cout << " | BC";
        else
            cout << " |   ";

        cout << " | " << nody[0][i] << " | " << nody[1][i] << endl;
    }
*/
    double ksi1[4] = {-1.0/sqrt(3), 1.0, 1.0/sqrt(3), -1.0};
    double ksi2[4] = { 1.0/sqrt(3), 1.0,  -1.0/sqrt(3), -1.0};

    double eta1[4] = {-1.0, -1.0/sqrt(3),  1.0, 1.0/sqrt(3)};
    double eta2[4] = {-1.0,  1.0/sqrt(3),  1.0,  -1.0/sqrt(3)};

    //cout << "Obliczanie wektora {P}:" << endl;
    for(int i=0;i<4;i++)
    {
        int pom1 = i;
        int pom2 = (i+1)%4; // Zabespieczenie, by indeksy
                            // nie wyszly poza tablice
        // Liczenie boku
        //if(1) // DEBUG
        if(nodeBC[pom1] && nodeBC[pom2])
        {
            //cout << "Bok " << 1+i << ": " << nodeId[pom1] << " | " << nodeId[pom2] << endl;

            double pomX1, pomX2, pomY1, pomY2;

            pomX1 = nody[0][pom1];
            pomX2 = nody[0][pom2];

            pomY1 = nody[1][pom1];
            pomY2 = nody[1][pom2];

            // Jakobian
            double detJAB = (sqrt(pow(pomX2-pomX1, 2)+pow(pomY2-pomY1, 2)))/2.0;

            double P_pc1[4] = { 0.0 };
            double P_pc2[4] = { 0.0 };
            double P_sum[4] = { 0.0 };

            //ksi[i*2]        eta[i*2]      // Pierwszy punkt
            //ksi[(i*2)+1]    eta[(i*2)+1]  // Drugi punktt

            double dane_tab[2][4] = { 0 };

            for(int j=0;j<4;j++)
            {
                double pomKsi, pomEta;
                //dane_tab[0][4] // tu ma byc petla po czterech
                //dane_tab[1][4] // funkcjach ksztaltu

                pomKsi = ksi1[j];
                pomEta = eta1[j];
                dane_tab[0][j] = funKsz_w(pomKsi, pomEta, i);

                pomKsi = ksi2[j];
                pomEta = eta2[j];
                dane_tab[1][j] = funKsz_w(pomKsi, pomEta, i);
            }
/*
            cout << "Dane: P1: "<< nodeId[0] <<" | P2: "<< nodeId[1] <<" | P3: "<< nodeId[2] <<" | P4: "<< nodeId[3] <<" |" << endl;
            for(int i=0;i<2;i++) // i - wiersz
            {
                for(int j=0;j<4;j++) // j - kolumna
                {
                    cout << " " << dane_tab[i][j] << "\t|";
                }
                cout << endl;
            }
            cout << endl;
*/
        // Mnozenie macierzy pc1
            for(int i=0;i<4;i++) // i - wiersz
            {
                P_pc1[i] = dane_tab[0][i] * wag1 * alpha * temperOt;
            }

        // Mnozenie macierzy pc2
            for(int i=0;i<4;i++) // i - wiersz
            {
                P_pc2[i] = dane_tab[1][i] * wag2 * alpha * temperOt;
            }

        // Obliczanie macierzy H_bc
            for(int i=0;i<4;i++) // i - wiersz
            {
                P_sum[i] = (P_pc1[i] + P_pc2[i]) * detJAB;
                vectorP[i] = vectorP[i] + P_sum[i];
            }
/*
            cout << "Wetor {P} lokalny:" << endl;
            for(int i=0;i<4;i++) // i - wiersz
            {
                cout << "| " << P_sum[i] << "\t";
            }
            cout << "|" << endl;
*/
        }
    }
/*
    cout << "Zsumowany wektor P:" << endl;
    for(int i=0;i<4;i++) // i - wiersz
    {
        cout << "| " << vectorP[i] << "\t";
    }
    cout << endl;
*/
    // Agregacja
    int Agreg[4];

    for(int i=0;i<4;i++)
    {
        Agreg[i] = nodeId[i];
    }

    //cout << "Macierz agregacji:" << endl;
    for(int i=0;i<4;i++)
    {
        int pomA = Agreg[i];// - 1;
        //cout << " " << Agreg[i] << " |";

        vecP[pomA-1] = vecP[pomA-1] + vectorP[i];

    }
    //cout << endl << endl;
}


int _matrixH::elimGauss()
{
    //int n, double ** AB, double * X
    //int rozmH, double ** macRoz, double * vecTemp
    int n = rozmH;

    int i, j, k;
    double m, s;

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            macRoz[i][j] = macH[i][j];
        }
    }

    for(int i=0;i<rozmH;i++)
    {
        macRoz[i][rozmH] = vecP[i];
    }

// eliminacja wspó³czynników

    for( i = 0; i < n - 1; i++ )
    {
        for( j = i + 1; j < n; j++ )
        {
            if( fabs ( macRoz[i][i] ) < eps )
            {
                return -1; // eliminacja nie powiodla sie
            }
            m = -macRoz[j][i] / macRoz[i][i];
            for( k = i + 1; k <= n; k++ )
            {
                macRoz[j][k] += m * macRoz[i][k];
            }
        }
    }

// wyliczanie niewiadomych

    for( i = n - 1; i >= 0; i-- )
    {
        s = macRoz[i][n];
        for( j = n - 1; j >= i + 1; j-- )
        s -= macRoz[i][j] * vecTemp[j];
        if( fabs ( macRoz[i][i] ) < eps )
        {
            return -1; // eliminacja nie powiodla sie
        }
        vecTemp[i] = s / macRoz[i][i];
    }
    return 0; // eliminacja powiodla sie
}


void _matrixH::stanNiestacionarny(int krok, int initemp)
{
    double stepTime = 1.0 * krok;
    double tempStart = 1.0 * initemp;
    int n = rozmH;

    // macierz C przez krok czasowy
    double ** pomC = new double *[n];
    for(int i=0;i<n;i++)
    {
        pomC[i] = new double[n];
        for(int j=0;j<n;j++)
        {
            pomC[i][j] = (macC[i][j])/stepTime;
        }
    }

    // Obliczona H + macierz C przez krok czasowy
    niestH = new double *[n];
    for(int i=0;i<n;i++)
    {
        niestH[i] = new double[n];
        for(int j=0;j<n;j++)
        {
            niestH[i][j] = macH[i][j] + pomC[i][j];
        }
    }

    // Inicjowanie pomocniczych T i P
    pomP = new double[n];
    pomTemp = new double[n];
    for(int i=0;i<n;i++)
    {
        pomP[i] = vecP[i];
        pomTemp[i] = tempStart;
    }

    // Wypisywanie wyniku
    cout << "---Stan_niestacionarny---" << endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout << " ";
            cout.precision(6);
            cout.width(8);
            cout.fill();
            cout << niestH[i][j] << " |";
            //cout << "\t|";
        }
        cout << endl << endl;
    }
    cout << "-------------------------" << endl;
}


int _matrixH::krokSymulacji(int krok, double * Temperatury)
{   // [H] + ([C]/krokCzas) * wecTempK = ( ([C]/krokCzas) * wecTempP ) - {P}
    // [A] * {X} = {B}

    int n = rozmH;
    double ** C_przez_krok; // = [C]/krokCzas
    double * C_razy_Temp;   // {B} = ( ([C]/krokCzas) * wecTempP ) - {P}
    double ** H_plus_C;     // [A] = [H] + ([C]/krokCzas)
    double * Temp_K;        // {X} = wecTempK

    // Wyliczanie [C]/krokCzas
    C_przez_krok = new double *[n];
    for(int i=0;i<n;i++)
    {
        C_przez_krok[i] = new double[n];
        for(int j=0;j<n;j++)
        {
            C_przez_krok[i][j] =  macC[i][j] / krok;
            //cout << " " << C_przez_krok[i][j];
        }
        //cout << endl;
    }

    // Obliczanie ( ([C]/krokCzas) * wecTempP ) + {P}
    C_razy_Temp = new double[n];
    for(int i=0;i<n;i++)
    {
        C_razy_Temp[i] = 0.0;
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            C_razy_Temp[i] = C_razy_Temp[i] + C_przez_krok[i][j] * Temperatury[j];
        }
    }

    for(int i=0;i<n;i++)
    {
        C_razy_Temp[i] = C_razy_Temp[i] + vecP[i];
        //cout << " " << C_razy_Temp[i];
    }
    //cout << endl;

    // Wyliczanie [H] + ([C]/krokCzas)
    H_plus_C = new double *[n];
    for(int i=0;i<n;i++)
    {
        H_plus_C[i] = new double[n];
        for(int j=0;j<n;j++)
        {
            H_plus_C[i][j] = macH[i][j] + C_przez_krok[i][j];
            //cout << " " << H_plus_C[i][j];
        }
        //cout << endl;
    }

    // Obliczanie {X} w rownaniu [A] * {X} = {B}
    Temp_K = new double[n];
    for(int i=0;i<n;i++)
    {
        Temp_K[i] = 0.0;
    }

    // Pomocnicza macierz rozszeerzona
    double ** pomMR = new double *[n];
    for(int i=0;i<n;i++)
    {
        pomMR[i] = new double[n+1];
        for(int j=0;j<n+1;j++)
        {
            pomMR[i][j] = 0.0;
        }
    }

    // Wczytywanie [A]
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            pomMR[i][j] = H_plus_C[i][j];
        }
    }

    // Wczytywanie {B}
    for(int i=0;i<n;i++)
    {
        pomMR[i][n] = C_razy_Temp[i];
    }

    // Elimiacja Gaussa
      // {B} = C_razy_Temp
      // [A] = H_plus_C
      // {X} = Temp_K
    double m, s;
    double eps = 1e-12;

  // eliminacja współczynników

    for(int i=0;i<n-1;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            if( fabs( pomMR[i][i] ) < eps )
                return false;
            m = -pomMR[j][i] / pomMR[i][i];
            for( int k=i+1;k<=n;k++)
            pomMR[j][k] += m * pomMR[i][k];
        }
    }

    // wyliczanie niewiadomych

    for(int i=n-1;i>=0;i--)
    {
        s = pomMR[i][n];
        for(int j=n-1;j>=i+1;j--)
        s -= pomMR[i][j] * Temp_K[j];
        if( fabs( pomMR [i][i] ) < eps )
            return false;
        Temp_K[i] = s / pomMR [i][i];
    }

    for(int i=0;i<n;i++)
    {
        Temperatury[i] = Temp_K[i];
        //cout << " " << Temp_K[i];
    }
    //cout << endl;

    // Wyszukiwanie najmniejszej i najwiekszej temperatury


    return 0;
}


/*
void _matrixH::obliczH_3P(double dane[2][4], int nodeId[4], double wag1, double wag2, int con, int dens, int speh)
{
    //Punkty calkowania:
    // Ksi: 0.0, -1.0*sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), sqrt(3.0/5.0), sqrt(3.0/5.0), 0.0, -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0)
    // Eta: 0.0, -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), sqrt(3.0/5.0), sqrt(3.0/5.0), 0.0

    double conductivity = 1.0 * con;
    double density = 1.0 * dens;
    double specheat = 1.0 * speh;
    //double conductivity = 25.0;

    double nody[2][4]; // 0 - X, 1 - Y
    for(int i=0;i<4;i++)
    {
        nody[0][i] = dane[0][i]; // X
        nody[1][i] = dane[1][i]; // Y
    }
/*
    cout << "id | X | Y" << endl;
    for(int i=0;i<4;i++)
    {
        cout << nodeId[i] << " | " << nody[0][i] << " | " << nody[1][i] << endl;
    }

    //cout << "Macierz H" << endl;

    double ksi[9] = {0.0, -1.0*sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), sqrt(3.0/5.0), sqrt(3.0/5.0), 0.0, -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0)};
    double eta[9] = {0.0, -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0), -1.0*sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), sqrt(3.0/5.0), sqrt(3.0/5.0), 0.0};

    //cout << " ksi | eta" << endl;
    //for(int i=0;i<4;i++)
    //{
    //    cout << " " << ksi[i] << " | " << eta[i] << endl;
    //}

    double poKsi[4][4];
    double poEta[4][4];

    for(int i=0;i<4;i++)
    {
        poKsi[i][0] = pFksiN_1(eta[i]);
        poKsi[i][1] = pFksiN_2(eta[i]);
        poKsi[i][2] = pFksiN_3(eta[i]);
        poKsi[i][3] = pFksiN_4(eta[i]);
    }

    for(int i=0;i<4;i++)
    {
        poEta[i][0] = pFetaN_1(ksi[i]);
        poEta[i][1] = pFetaN_2(ksi[i]);
        poEta[i][2] = pFetaN_3(ksi[i]);
        poEta[i][3] = pFetaN_4(ksi[i]);
    }

    //for(int i=0;i<4;i++)
    //{
    //    for(int j=0;j<4;j++)
    //    {
    //        cout << " " << poEta[i][j] << " |";
    //    }
    //    cout << endl;
    //}

    //double nody[2][4] = {// 0 - X, 1 - Y
    //    {0.0, 0.025, 0.025, 0.0},
    //    {0.0, 0.0, 0.025, 0.025}
    //};

    double mJakob[2][2] = {0.0};

    // dX / dKsi
    mJakob[0][0] = (nody[0][0]*poKsi[0][0]) + (nody[0][1]*poKsi[0][1]) + (nody[0][2]*poKsi[0][2]) + (nody[0][3]*poKsi[0][3]);
    // dY / dKsi
    mJakob[0][1] = (nody[1][0]*poKsi[0][0]) + (nody[1][1]*poKsi[0][1]) + (nody[1][2]*poKsi[0][2]) + (nody[1][3]*poKsi[0][3]);
    // dX / dEta
    mJakob[1][0] = (nody[0][0]*poEta[0][0]) + (nody[0][1]*poEta[0][1]) + (nody[0][2]*poEta[0][2]) + (nody[0][3]*poEta[0][3]);
    // dY / dEta
    mJakob[1][1] = (nody[1][0]*poEta[0][0]) + (nody[1][1]*poEta[0][1]) + (nody[1][2]*poEta[0][2]) + (nody[1][3]*poEta[0][3]);

    //cout << mJakob[0][0] << " | " << mJakob[0][1] << endl;
    //cout << mJakob[1][0] << " | " << mJakob[1][1] << endl;

    double detJakob = (mJakob[0][0]*mJakob[1][1]) - (mJakob[0][1]*mJakob[1][0]);
    //cout << "detJ = " << detJakob << endl;

    double odwJakob[2][2];
    odwJakob[0][0] =  mJakob[1][1] * (1/detJakob);
    odwJakob[0][1] = -mJakob[0][1] * (1/detJakob);
    odwJakob[1][0] = -mJakob[1][0] * (1/detJakob);
    odwJakob[1][1] =  mJakob[0][0] * (1/detJakob);

    double pcX[4][4]; // pochodne funkcji ksztaltu po X
    double pcY[4][4]; // pochodne funkcji ksztaltu po Y

    for(int i=0;i<4;i++) // pcX
    {
        for(int j=0;j<4;j++)
        {
            pcX[i][j] = odwJakob[0][0]*poKsi[i][j] + odwJakob[0][1]*poEta[i][j];
        }
    }

    for(int i=0;i<4;i++) // pcY
    {
        for(int j=0;j<4;j++)
        {
            pcY[i][j] = odwJakob[1][0]*poKsi[i][j] + odwJakob[1][1]*poEta[i][j];
        }
    }
/*
    cout << "pcX" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << pcX[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    cout << "pcY" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << pcY[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    double localH[4][4] = {0.0};
    double localC[4][4] = {0.0};
    double tempH[4][4] = {0.0};
    double tempC[4][4] = {0.0};
    double tempCrazyC[4][4] = {0.0};
    //double w1 = 1.0, w2 = 1.0;

    double ksiC[4] = {-1.0, -1.0, 1.0, 1.0};
    double etaC[4] = {-1.0, 1.0, -1.0, 1.0};
    for(int i=0;i<4;i++)
    {
        ksiC[i] = ksiC[i] * 1.0/sqrt(3);
        etaC[i] = etaC[i] * 1.0/sqrt(3);
    }

    for(int i=0;i<4;i++)
    {
        tempC[i][0] = funKsz1(ksiC[i], etaC[i]);
        tempC[i][1] = funKsz2(ksiC[i], etaC[i]);
        tempC[i][2] = funKsz3(ksiC[i], etaC[i]);
        tempC[i][3] = funKsz4(ksiC[i], etaC[i]);
    }

    //cout << " {C} Funkcji ksztaltu: " << endl;
    //for(int i=0;i<4;i++)
    //{
    //    cout << tempC[i][0] << " | " << tempC[i][1] << " | " << tempC[i][2] << " | " << tempC[i][3] << endl;
    //}

    //for(int i=0;i<4;i++)
    //{
    //    for(int j=0;j<4;j++)
    //    {
    //        localC[i][j] = tempC[i][j] * tempC[j][i] * detJakob * density * specheat;
    //    }
    //}

    double pcXrazyT[4][4]; // iloczyn transponowanej i pcX
    double pcYrazyT[4][4]; // iloczyn transponowanej i pcY

    for(int n=0;n<4;n++)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {   // Mno¿enie macierzy przez nie same transponowane
                pcXrazyT[i][j] = pcX[n][i]*pcX[n][j];
                pcYrazyT[i][j] = pcY[n][i]*pcY[n][j];
                tempCrazyC[i][j] = tempC[n][i] * tempC[n][j];
            }
        }


        for(int i=0;i<4;i++)
        {

            for(int j=0;j<4;j++)
            {
                tempH[i][j] = (pcXrazyT[i][j] + pcYrazyT[i][j]) * conductivity * detJakob;
                localH[i][j] = localH[i][j] + tempH[i][j];
                //localH[i][j] = localH[i][j] * wag1 * wag2;

                // Macierz C
                //tempCrazyC[i][j] = tempC[i][j] * tempC[j][i]; // Razy transponowana
                tempCrazyC[i][j] = tempCrazyC[i][j] * detJakob * density * specheat;
                localC[i][j] = localC[i][j] + tempCrazyC[i][j];
                //tempC[i][j] = (pcXrazyT[i][j] + pcYrazyT[i][j]) * conductivity * detJakob;
                //localC[i][j] = localC[i][j] + tempC[i][j];
            }
        }

        //for(int i=0;i<4;i++)
        //{
        //    for(int j=0;j<4;j++)
        //    {
        //         cout << " " << tempH[i][j] << " |";
        //    }
        //    cout << endl;
        //}
        //cout << endl;
    }
/*
    cout << "DEBUG: Lokalna H" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << localH[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    cout << "DEBUG: Lokalna C" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout << " " << localC[i][j] << " |";
        }
        cout << endl;
    }
    cout << endl;

    // Agregacja
    int Agr_x[4][4];
    int Agr_y[4][4];

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            Agr_x[i][j] = nodeId[i];
            Agr_y[i][j] = nodeId[j];
        }
    }

    //cout << "Macierz agregacji:" << endl;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            int pomX = Agr_x[i][j];// - 1;
            int pomY = Agr_y[i][j];// - 1;
            //cout << " " << Agr_x[i][j] << "," << Agr_y[i][j] << " |";

            macH[pomX-1][pomY-1] =  macH[pomX-1][pomY-1] + localH[i][j];    // Macierz H
            macC[pomX-1][pomY-1] =  macC[pomX-1][pomY-1] + localC[i][j];    // Macierz C
        }
        //cout << endl;
    }
    //cout << endl << endl;
}
*/
