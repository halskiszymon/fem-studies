//
// Created by Szymon Halski on 23/10/2024.
//

#include "structures.h"
#include <vector>
#include <cmath>

using namespace std;

#ifndef JAKOBIAN_H
#define JAKOBIAN_H

struct ElementJakobian {
    vector<Node> nodes;

    // lab 10/11 - macierz C lokalna
    vector<vector<double>> C_local;

    // lab 7 - macierz Hbc
    vector<vector<double>> Hbc_local;

    double dlugoscBoku(int sideIndex) {
        int node1, node2;
        if (sideIndex == 0) { node1 = 0; node2 = 1; }
        else if (sideIndex == 1) { node1 = 1; node2 = 2; }
        else if (sideIndex == 2) { node1 = 2; node2 = 3; }
        else if (sideIndex == 3) { node1 = 3; node2 = 0; }
        else {
            cerr << "Nieprawidlowy index boku: " << sideIndex << endl;
            return 0.0;
        }

        double dx = nodes[node2].x - nodes[node1].x;
        double dy = nodes[node2].y - nodes[node1].y;

        return sqrt(dx * dx + dy * dy);
    }

    // lab 4 - macierz H lokalna
    vector<vector<double>> H_local;

    ElementJakobian() {
        H_local = vector<vector<double>>(4, vector<double>(4, 0.0));
        Hbc_local = vector<vector<double>>(4, vector<double>(4, 0.0));
        C_local = vector<vector<double>>(4, vector<double>(4, 0.0));
    }

    double dN_dx[4], dN_dy[4];

    // pochodne funkcji kształtu w punkcie Gaussa
    void pochodneFunkcjiKsztaltu(double ksi, double eta, double N[4], double dN_dKsi[4], double dN_dEta[4]) {
        N[0] = 0.25 * (1 - ksi) * (1 - eta);
        N[1] = 0.25 * (1 + ksi) * (1 - eta);
        N[2] = 0.25 * (1 + ksi) * (1 + eta);
        N[3] = 0.25 * (1 - ksi) * (1 + eta);

        if (dN_dKsi != nullptr) {
            dN_dKsi[0] = -0.25 * (1 - eta);
            dN_dKsi[1] =  0.25 * (1 - eta);
            dN_dKsi[2] =  0.25 * (1 + eta);
            dN_dKsi[3] = -0.25 * (1 + eta);
        }

        if (dN_dEta != nullptr) {
            dN_dEta[0] = -0.25 * (1 - ksi);
            dN_dEta[1] = -0.25 * (1 + ksi);
            dN_dEta[2] =  0.25 * (1 + ksi);
            dN_dEta[3] =  0.25 * (1 - ksi);
        }
    }

    // obliczanie macierzy jakobianu
    vector<vector<double>> jakobian(double dN_dKsi[4], double dN_dEta[4]) {
        vector<vector<double>> J(2, vector<double>(2, 0.0));

        for (int i = 0; i < 4; i++) {
            J[0][0] += dN_dKsi[i] * nodes[i].x;
            J[0][1] += dN_dKsi[i] * nodes[i].y;
            J[1][0] += dN_dEta[i] * nodes[i].x;
            J[1][1] += dN_dEta[i] * nodes[i].y;
        }

        return J;
    }

    // obliczanie wyznacznika macierzy jakobianu
    double wyznacznik(const vector<vector<double>>& J) {
        return J[0][0] * J[1][1] - J[0][1] * J[1][0];
    }

    // obliczanie odwrotnosci macierzy jakobianu
    vector<vector<double>> odwrotnoscJakobian(const vector<vector<double>>& J, double detJ) {
        vector<vector<double>> invJ(2, vector<double>(2, 0.0));

        invJ[0][0] = J[1][1] / detJ;
        invJ[0][1] = -J[0][1] / detJ;
        invJ[1][0] = -J[1][0] / detJ;
        invJ[1][1] = J[0][0] / detJ;

        return invJ;
    }

    // obliczanie pochodnych funkcji ksztaltu dla odwrotnosci Jakobiana
    void pochodneFunkcjiKsztaltuOdwrotnosc(vector<vector<double>> invJ, double dN_dKsi[4], double dN_dEta[4], bool wypiszWynik = false) {
        for (int k = 0; k < 4; k++) {
            dN_dx[k] = invJ[0][0] * dN_dKsi[k] + invJ[0][1] * dN_dEta[k];
            dN_dy[k] = invJ[1][0] * dN_dKsi[k] + invJ[1][1] * dN_dEta[k];
        }

        if (wypiszWynik) {
            cout << "dN/dx:" << endl;
            cout << "    ";
            for (int k = 0; k < 4; k++) {
                cout << dN_dx[k] << " ";
            }
            cout << endl;

            cout << "dN/dy:" << endl;
            cout << "    ";
            for (int k = 0; k < 4; k++) {
                cout << dN_dy[k] << " ";
            }
            cout << endl << endl;
        }
    }
};

#endif //JAKOBIAN_H
