//
// Created by Szymon Halski on 10/10/2024.
//

#ifndef METODA_GAUSSA_H
#define METODA_GAUSSA_H

#include <iostream>
#include <cmath>

using namespace std;

struct PointsWeightTable {
	double* point;
	double* weight;

	PointsWeightTable(int size) {
		point = new double[size];
		weight = new double[size];
	}
};

double metodaGaussa2D(int points, double (*func)(double)) {
	PointsWeightTable pwt = PointsWeightTable(points);

	if (points == 2) {
		pwt.point[0] = -sqrt(1.0 / 3.0);
		pwt.point[1] = sqrt(1.0 / 3.0);

		pwt.weight[0] = 1.0;
		pwt.weight[1] = 1.0;
	} else if (points == 3) {
		pwt.point[0] = -sqrt(3.0 / 5.0);
		pwt.point[1] = 0.0;
		pwt.point[2] = sqrt(3.0 / 5.0);

		pwt.weight[0] = 5.0 / 9.0;
		pwt.weight[1] = 8.0 / 9.0;
		pwt.weight[2] = 5.0 / 9.0;
	} else if (points == 4) {
		pwt.point[0] = -0.8611363115940526;
		pwt.point[1] = -0.3399810435848563;
		pwt.point[2] = 0.3399810435848563;
		pwt.point[3] = 0.8611363115940526;

		pwt.weight[0] = 0.3478548451374539;
		pwt.weight[1] = 0.6521451548625461;
		pwt.weight[2] = 0.6521451548625461;
		pwt.weight[3] = 0.3478548451374539;
	} else {
		cerr << "Nieobsługiwana liczba punktów: " << points << endl;
		return 0.0;
	}

	double suma = 0;

	for (int i = 0; i < points; i++) {
		suma += pwt.weight[i] * func(pwt.point[i]);
	}

	return suma;
}

double metodaGaussa3D(int points, double (*func)(double, double)) {
	PointsWeightTable pwt = PointsWeightTable(points);

	if (points == 2) {
		pwt.point[0] = -sqrt(1.0 / 3.0);
		pwt.point[1] = sqrt(1.0 / 3.0);

		pwt.weight[0] = 1.0;
		pwt.weight[1] = 1.0;
	} else if (points == 3) {
		pwt.point[0] = -sqrt(3.0 / 5.0);
		pwt.point[1] = 0.0;
		pwt.point[2] = sqrt(3.0 / 5.0);

		pwt.weight[0] = 5.0 / 9.0;
		pwt.weight[1] = 8.0 / 9.0;
		pwt.weight[2] = 5.0 / 9.0;
	} else if (points == 4) {
		pwt.point[0] = -0.8611363115940526;
		pwt.point[1] = -0.3399810435848563;
		pwt.point[2] = 0.3399810435848563;
		pwt.point[3] = 0.8611363115940526;

		pwt.weight[0] = 0.3478548451374539;
		pwt.weight[1] = 0.6521451548625461;
		pwt.weight[2] = 0.6521451548625461;
		pwt.weight[3] = 0.3478548451374539;
	} else {
		cerr << "Nieobsługiwana liczba punktów: " << points << endl;
		return 0.0;
	}

	double suma = 0;

	for (int i = 0; i < points; i++) {
		for (int j = 0; j < points; j++) {
			suma += pwt.weight[i] * pwt.weight[j] * func(pwt.point[i], pwt.point[j]);
		}
	}

	return suma;
}

#endif //METODA_GAUSSA_H