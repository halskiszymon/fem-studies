//
// Created by Szymon Halski on 03/10/2024.
//

#ifndef STRUCTURES_H
#define STRUCTURES_H

struct Node {
    double x;
    double y;
    bool isBoundary; // do macierzy Hbc - czy jest na granicy

    Node(double x1 = 0.0, double y1 = 0.0) : x(x1), y(y1), isBoundary(false) {};
};

struct Element {
    // ID[1x4] - identyfikatory tablicy np {1, 2, 6, 5}
    int ID[4];
};

struct Grid {
    int nN = 0; // liczba wezlow
    int nE = 0; // liczba elementow

    Node* nodes; // wezly
    Element* elements; // elementy
};

// do odczytu z pliku
struct GlobalData {
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
};

#endif //STRUCTURES_H
