#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "siatki.h"

using namespace std;

_Grid::_Grid(string adres)
{
    fstream plik;
    plik.open( adres, ios::in );

    string dane;
    for( int i=0; i<10; i++ )
    {
        int j, sp;
        getline( plik, dane );
        //cout << "[" << dane << "]" << endl;   // Debugowanie

        int n = dane.size();
        //cout << "[" << n << "]" << endl;    // Debugowanie

        for( j=n; j>0; j-- )
        {
            //cout << " >" << j << endl;      // Debugowanie

            sp = dane[j-1];
            //cout << " [" << dane[j-1] << "][" << sp << "]" << endl;     // Debugowanie

            if(sp==32) // Szukamy " "
                break;
        }
        string liczba = dane.substr(j);
        //cout << "[" << liczba << "]" << endl;   // Debugowanie

        varTab[i] = atoi(liczba.c_str());
        //cout << varTab[i] << endl;      // Debugowanie
    }

    //for( int i=0; i<10; i++ )       //
    //{                               // Debugowanie
    //    cout << varTab[i] << endl;  //
    //}                               //


// Wczytywanie danych do _Node
    getline( plik, dane );
    //cout << "[" << dane << "]" << endl;   // Debugowanie

    int m = varTab[8];
    listaN = new _Node[m];
    _Node* wskPom;
    for(int i=0;i<m;i++)
    {
        getline( plik, dane );
        //cout << "[" << dane << "]" << endl;   // Debugowanie

        int n, j, sp;
        string liczba;
        wskPom = &(listaN[i]);

    // Wczytywanie Y
        n = dane.size();
        for(j=n;j>0;j--)
        {
            sp = dane[j-1];
            if(sp==44)  // Szukamy ","
                break;
        }
        liczba = dane.substr(j);
        wskPom->y = stod(liczba);
        //cout << "Y = " << wskPom->y << endl;    // Debugowanie

        dane.erase(j-1);
        //cout << "[" << dane << "]" << endl;     // Debugowanie

    // Wczytywanie X
        n = dane.size();
        for(j=n;j>0;j--)
        {
            sp = dane[j-1];
            if(sp==44)  // Szukamy ","
                break;
        }
        liczba = dane.substr(j);
        wskPom->x = stod(liczba);
        //cout << "X = " << wskPom->x << endl;    // Debugowanie

        dane.erase(j-1);
        //cout << "[" << dane << "]" << endl;     // Debugowanie

    // Wczytywanie ID                         // Do naprawienia
        wskPom->id = atoi(dane.c_str());      // Juz naprawione
        //cout << "ID = " << wskPom->id << endl;    // Debugowanie
    }

    //for(int i=0;i<m;i++)            //
    //{                               // Debugowanie
    //    wskPom = &(listaN[i]);      //
    //    cout << wskPom->id << ":\t";//
    //    cout << wskPom->x << ", ";  //
    //    cout << wskPom->y << endl;  //
    //}                               //

// Wczytywanie danych do _Element
    getline(plik, dane);
    //cout << "[" << dane << "]" << endl;   // Debugowanie

    m = varTab[9];
    listaE = new _Element[m];
    _Element* wskEl;

    for(int i=0;i<m;i++)
    {
        getline( plik, dane );
        //cout << "[" << dane << "]" << endl;   // Debugowanie


        int n, j, sp;
        string liczba;
        wskEl = &(listaE[i]);

    // Wczytywanie ID[k]
        for(int k=3;k>=0;k--)
        {
            n = dane.size();
            for(j=n;j>0;j--)
            {
                sp = dane[j-1];
                if(sp==44)  // Szukamy ","
                    break;
            }
            liczba = dane.substr(j);
            wskEl->ID[k] = atoi(liczba.c_str());
            //cout << "ID[" << k << "] = " << wskEl->ID[k] << endl;    // Debugowanie
            dane.erase(j-1);
        }
    }

    //for(int i=0;i<m;i++)                //
    //{                                   //
    //    wskEl = &(listaE[i]);           //
    //    cout << i+1 << ": ";              // Debugowanie
    //    cout << wskEl->ID[0] << ", ";   //
    //    cout << wskEl->ID[1] << ", ";   //
    //    cout << wskEl->ID[2] << ", ";   //
    //    cout << wskEl->ID[3] << endl;   //
    //}                                   //

//Wczytywanie do BC
    getline( plik, dane );
    getline( plik, dane );
    //cout << "[" << dane << "]" << endl;   // Debugowanie

    string liczba;

    int i, sp, n = dane.size();
    bool koniec = true;


    for(i=n;i>0;i--)
    {
        sp = dane[i-1];
        //cout << i << ": " << sp << endl;    //debugowanie
        if(sp==44)  // Szukamy ","
        {
            liczba = dane.substr(i+1);
            int nodenum = atoi(liczba.c_str());
            wskPom = &(listaN[nodenum-1]);
            wskPom->bc = true;

            //cout << "[ ! ]:" << liczba << endl; // Debugwanie
            //cout << "[ i ]:" << nodenum << endl; // Debugwanie
            dane.erase(i-1);
        }
    }
    liczba = dane.substr(i);
    int nodenum = atoi(liczba.c_str());
    wskPom = &(listaN[nodenum-1]);
    wskPom->bc = true;

    //cout << "[ ! ]:" << liczba << endl; // Debugwanie
    //cout << "[ i ]:" << nodenum << endl; // Debugwanie

    //for(i=0;i<varTab[8];i++)                //
    //{                               // Debugowanie
    //    wskPom = &(listaN[i]);      //
    //    cout << i << ": ";          //
    //    if(wskPom->bc)              //
    //        cout << "true" << endl; //
    //    else                        //
    //        cout << "false" << endl;//
    //}                               //

    plik.close();
}

void _Grid::wypisz()
{
    _Node* wskPom;
    _Element* wskEl;
    int m = varTab[8];
    int n = varTab[9];

    cout << "Dane:" << endl;
    for( int i=0; i<10; i++ )
    {
        cout << varTab[i] << endl;
    }
    cout << endl;

    cout << "Nody: " << varTab[8] << endl;
    for(int i=0;i<m;i++)
    {
        wskPom = &(listaN[i]);
        cout << wskPom->id << ":\t";
        cout << wskPom->x << ", ";
        cout << wskPom->y << endl;
    }
    cout << endl;

    cout << "Elementy: " << varTab[9] << endl;
    for(int i=0;i<n;i++)
    {
        wskEl = &(listaE[i]);
        cout << i+1 << ": ";
        cout << wskEl->ID[0] << ", ";
        cout << wskEl->ID[1] << ", ";
        cout << wskEl->ID[2] << ", ";
        cout << wskEl->ID[3] << endl;
    }
    cout << endl;

    cout << "BC:" << endl;
    for(int i=0;i<m;i++)
    {
        wskPom = &(listaN[i]);
        if(wskPom->bc)
            cout << " " << i+1;
    }
    cout << endl << endl;
}
