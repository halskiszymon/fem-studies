//
// Created by Szymon Halski on 06/11/2024.
//
#pragma once

#include <iostream>
#include <string>
#include <vector>

#ifndef POMOCNICZE_H
#define POMOCNICZE_H

using namespace std;

// obliczanie macierzy C lokalnej
void obliczMacierzC(ElementJakobian& elem, const double points[], const double weights[], int numPoints, double density, double specificHeat) {
    double N[4];

    for (int i = 0; i < numPoints; i++) {
        // pętla po ksi
        double ksi = points[i];
        double weightKsi = weights[i];

        for (int j = 0; j < numPoints; j++) {
            // pętla po eta
            double eta = points[j];
            double weightEta = weights[j];

            // wyznacznik wagi (2D) - iloczyn wag Gaussa
            double weight = weightKsi * weightEta;

            double dN_dKsi[4], dN_dEta[4];
            elem.pochodneFunkcjiKsztaltu(ksi, eta, N, dN_dKsi, dN_dEta);

            vector<vector<double>> J = elem.jakobian(dN_dKsi, dN_dEta);
            double detJ = elem.wyznacznik(J);

            for (int m = 0; m < 4; m++) {
                for (int n = 0; n < 4; n++) {
                    elem.C_local[m][n] += N[m] * N[n] * density * specificHeat * detJ * weight;
                }
            }
        }
    }
}

// funkcja rozwiązująca układ równań liniowych Ax = b gauusem bez czesciowego wyboru elementow
vector<double> eliminacjaGaussa(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // eliminacja w przod
    for (int k = 0; k < n; k++) {
        // normalizacja elementu glownego
        double pivot = A[k][k];
        for (int j = k; j < n; j++) {
            A[k][j] /= pivot;
        }
        b[k] /= pivot;

        // zerowanie elementow
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // podstawiamy wstecz
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }

    return x;
}

// obliczanie macierzy Hbc
void obliczHbcDlaBoku(ElementJakobian& elem, int sideIndex, const double points[], const double weights[], int numPoints, double alfa) {
    if (elem.Hbc_local.empty()) {
        elem.Hbc_local = vector<vector<double>>(4, vector<double>(4, 0.0));
    }

    // dlugosc boku
    double detJ = elem.dlugoscBoku(sideIndex) / 2.0;

    for (int i = 0; i < numPoints; i++) {
        double ksi, eta;
        double weight = weights[i];

        // okreslamy ksi i eta dla kazdego boku
        if (sideIndex == 0) {
            eta = -1.0;
            ksi = points[i];
        } else if (sideIndex == 1) {
            ksi = 1.0;
            eta = points[i];
        } else if (sideIndex == 2) {
            eta = 1.0;
            ksi = -points[i];
        } else if (sideIndex == 3) {
            ksi = -1.0;
            eta = -points[i];
        } else {
            cerr << "Nieprawidlowy indeks: " << sideIndex << endl;
            return;
        }

        // okreslamy pochodne funkcji ksztaltu w ksi i eta
        double N[4];
        elem.pochodneFunkcjiKsztaltu(ksi, eta, N, nullptr, nullptr);

        // obliczamy hbc
        for (int m = 0; m < 4; m++) {
            for (int n = 0; n < 4; n++) {
                elem.Hbc_local[m][n] += N[m] * N[n] * alfa * weight * detJ;
            }
        }
    }
}

// obliczanie macierzy H lokalnej
void obliczMacierzH(ElementJakobian& elem, const double points[], const double weights[], int numPoints, double conductivity) {
    elem.H_local = vector<vector<double>>(4, vector<double>(4, 0.0));

    double N[4];

    for (int i = 0; i < numPoints; i++) { // ksi
        double ksi = points[i];
        double weightKsi = weights[i];
        for (int j = 0; j < numPoints; j++) { // eta
            double eta = points[j];
            double weightEta = weights[j];
            double weight = weightKsi * weightEta;

            double dN_dKsi[4], dN_dEta[4];
            elem.pochodneFunkcjiKsztaltu(ksi, eta, N, dN_dKsi, dN_dEta);

            vector<vector<double>> J = elem.jakobian(dN_dKsi, dN_dEta);
            double detJ = elem.wyznacznik(J);
            vector<vector<double>> invJ = elem.odwrotnoscJakobian(J, detJ);

            elem.pochodneFunkcjiKsztaltuOdwrotnosc(invJ, dN_dKsi, dN_dEta);

            for (int m = 0; m < 4; m++) {
                for (int n = 0; n < 4; n++) {
                    elem.H_local[m][n] += (elem.dN_dx[m] * elem.dN_dx[n] + elem.dN_dy[m] * elem.dN_dy[n]) * detJ * conductivity * weight;
                }
            }
        }
    }
}

void obliczWypiszJakobian(ElementJakobian& elem, const double points[], int numPoints, const string& elementName) {
    int pc = 0;
    cout << endl << elementName << endl;

    for (int i = 0; i < numPoints; i++) { // ksi
        double ksi = points[i];
        for (int j = 0; j < numPoints; j++) { // eta
            double eta = points[j];
            pc++;
            cout << "Punkt całkowania nr. " << pc << endl;

            double N[4];
            double dN_dKsi[4], dN_dEta[4];
            elem.pochodneFunkcjiKsztaltu(ksi, eta, N, dN_dKsi, dN_dEta);

            vector<vector<double>> J = elem.jakobian(dN_dKsi, dN_dEta);
            double detJ = elem.wyznacznik(J);
            vector<vector<double>> invJ = elem.odwrotnoscJakobian(J, detJ);

            cout << "Punkt Gaussa (ksi, eta): (" << ksi << ", " << eta << ")" << endl;
            cout << "Macierz Jakobiana:" << endl;
            for (const auto& row : J) {
                cout << "    ";
                for (double val : row) {
                    cout << val << " ";
                }
                cout << endl;
            }

            cout << "Wyznacznik Jakobiana: " << detJ << endl;

            cout << "Odwrotność Jakobiana:" << endl;
            for (const auto& row : invJ) {
                cout << "    ";
                for (double val : row) {
                    cout << val << " ";
                }
                cout << endl;
            }

            // Obliczanie dN/dx i dN/dy
            elem.pochodneFunkcjiKsztaltuOdwrotnosc(invJ, dN_dKsi, dN_dEta, true);
        }
    }
}

#endif //POMOCNICZE_H
