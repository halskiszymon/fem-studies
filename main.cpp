#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include "metodaGaussa.h"
#include "structures.h"
#include "jakobian.h"
#include "pomocnicze.h"

using namespace std;

// defs
int getDataFromFile(const string &fileName, GlobalData &globalData, Grid &grid);
void displayGlobalData(const GlobalData &globalData);
void displayGrid(Grid grid);

// lab 2
double f1(double x) {
	return 2.0 * pow(x, 2) + 0.1 * x + 3.0;
}

double f2(double x, double y) {
	return -2.0 * pow(x, 2) * y + 2.0 * x * y + 4.0;
}

double f3(double x, double y) {
	return -5.0 * pow(x, 2) * y + 2.0 * x * y + 10.0;
}

// main
int main() {
	cout << fixed << setprecision(8);

	// lab 1
	// struktura wezla i siatki
	// wczytywanie z pliku siatki
	// tworzenie siatki z wczytanego pliku
	cout << "--- lab 1 ---" << endl;
	cout << "Aktualny dir: "  << filesystem::current_path() << endl;

	GlobalData globalData = GlobalData();
	Grid grid = Grid();

	// int result = getDataFromFile("/Users/szymonhalski/CLionProjects/mes/Test1_4_4.txt", globalData, grid);
	// int result = getDataFromFile("/Users/szymonhalski/CLionProjects/mes/Test2_4_4_MixGrid.txt", globalData, grid);
	int result = getDataFromFile("/Users/szymonhalski/CLionProjects/mes/Test3_31_31_kwadrat.txt", globalData, grid);

	if (result == 0) {
		cerr << "Wystapil blad podczas importowania danych pliku." << endl;
		return 0;
	}

	cout << "Import zakonczony." << endl << endl;

	cout << "Global data:" << endl;
	displayGlobalData(globalData);

	cout << endl << "Siatka:" << endl;
	displayGrid(grid);

	// lab 2
	// funkcja z metoda gauusa ktora w ukladzie -1,1 realizuje calkowanie 2 i 3 punktowe
	// cout << endl << "--- lab 2 ---" << endl;
	// cout << "f1 - Wynik działania dla 2d - 2 punktów: " << metodaGaussa2D(2, f1) << endl;
	// cout << "f1 - Wynik działania dla 2d - 3 punktów: " << metodaGaussa2D(3, f1) << endl;
	//
	// cout << "f2 - Wynik działania dla 3d - 2 punktów: " << metodaGaussa3D(2, f2) << endl;
	// cout << "f2 - Wynik działania dla 3d - 3 punktów: " << metodaGaussa3D(3, f2) << endl;
	//
	// cout << "f3 - Wynik działania dla 3d - 2 punktów: " << metodaGaussa3D(2, f3) << endl;
	// cout << "f3 - Wynik działania dla 3d - 3 punktów: " << metodaGaussa3D(3, f3) << endl;

	// lab 3
	// macierz jakobiana
	// obliczanie wyznacznika macierzy jakobiana oraz jej odwrotnosci
	// cout << endl << "--- lab 3 ---" << endl;
	// ElementJakobian elem1, elem2, elem3;
	// elem1.nodes = {{0, 0}, {0.025, 0}, {0.025, 0.025}, {0, 0.025}};
 //    elem2.nodes = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
 //    elem3.nodes = {{0, 0}, {4, 0}, {4, 4}, {0, 5}};

	// pochodne funkcji kształtu w punkcie
	// maja byc obliczane przy uzyciu gaussa
	// double dN_dKsi[4] = {-0.25, 0.25, 0.25, -0.25};
	// double dN_dEta[4] = {-0.25, -0.25, 0.25, 0.25};

	// punkty i wagi dla 2 punktowego calkowania
	double points2[2] = {-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};
	double weights2[2] = {1.0, 1.0};

	// punkty i wagi dla 3 punktowego calkowania
	double points3[3] = {-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};
	double weights3[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

	// lab 7
	// calkowanie 4 punktowe - modyfikacja metodaGauusa
	// dodanie macierzy Hbc
	// punkty i wagi dla 4 punktowego calkowania
	double points4[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
	double weights4[4] = {0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539};

	// Element 1
	// obliczWypiszJakobian(elem1, points2, 2, "Element 1");
	//
	// Element 2
	// obliczWypiszJakobian(elem2, points2, 2, "Element 2");
	//
	// Element 3
	// obliczWypiszJakobian(elem3, points2, 2, "Element 3");

	// lab 4
	// cout << endl << "--- lab 4 ---" << endl;
	// cout << "--- Obliczanie macierzy H dla przykladow z pdf ---" << endl;

	// przyklad z pdf nr. 6 - przyklad nr 1
	// ElementJakobian elem6_1;
	// elem6_1.nodes = {{0, 0}, {0.025, 0}, {0.025, 0.025}, {0, 0.025}};
	//
	// obliczWypiszJakobian(elem6_1, points3, 3, "Element pdf 6. nr 1");
	//
	// obliczMacierzH(elem6_1, points3, weights3, 3, 30);
	//
	// cout << "Macierz H dla elementu pdf 6. nr 1:" << endl;
	// for (const auto& row : elem6_1.H_local) {
	// 	for (double val : row) {
	// 		cout << val << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << endl;

	// przyklad z pdf nr. 6 - przyklad nr 2
	// ElementJakobian elem6_2;
	// elem6_2.nodes = {{0.01, -0.01}, {0.025, 0}, {0.025, 0.025}, {0, 0.025}};
	//
	// obliczWypiszJakobian(elem6_2, points2, 2, "Element pdf 6. nr 2");
	//
	// obliczMacierzH(elem6_2, points2, weights2, 2, 30);
	//
	// cout << "Macierz H dla elementu pdf 6. nr 2:" << endl;
	// for (const auto& row : elem6_2.H_local) {
	// 	for (double val : row) {
	// 		cout << val << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << endl;

	cout << endl << "--- lab 4.5 ---" << endl;
	cout << "--- Obliczanie Jakobianów dla elementów z pliku ---" << endl;

	// punkty i wagi dla 2 punktowego całkowania
	// double points2[2] = {-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};

	// macierz H globalna lab6
	vector<vector<double>> H_global(grid.nN, vector<double>(grid.nN, 0.0));

	// lab 8/9 - wektor P
	vector<double> P_global(grid.nN, 0.0);

	// lab 10/11 - macierz C
	vector<vector<double>> C_global(grid.nN, vector<double>(grid.nN, 0.0));

	for (int e = 0; e < grid.nE; e++) {
		Element currentElement = grid.elements[e];

		// pobierz wspolrzedne wezlow dla aktualnego elementu
		ElementJakobian elemJakobian;
		elemJakobian.nodes.clear();

		for (int i = 0; i < 4; i++) {
			int nodeIndex = currentElement.ID[i] - 1; // indexy w pliku zaczynają się od 1, więc odejmujemy 1
			Node node = grid.nodes[nodeIndex];
			elemJakobian.nodes.push_back(node);
		}

		string elementName = "Element " + to_string(e + 1);
		obliczWypiszJakobian(elemJakobian, points2, 2, elementName);

		obliczMacierzH(elemJakobian, points2, weights2, 2, globalData.conductivity);

		// wektor P lokalny
		vector<double> P_local(4, 0.0);

		// lab 7 - macierz Hbc
		for (int side = 0; side < 4; side++) {
			int node1Index = currentElement.ID[side] - 1;
			int node2Index = currentElement.ID[(side + 1) % 4] - 1;

			Node node1 = grid.nodes[node1Index];
			Node node2 = grid.nodes[node2Index];

			if (node1.isBoundary && node2.isBoundary) {
				// oblicz macierz Hbc dla tego boku
				obliczHbcDlaBoku(elemJakobian, side, points2, weights2, 2, globalData.alfa);

				// lab 8/9
				// obliczamy wektor P lokalny dla tego boku
				// N[m] * alfa * Tot

				double detJ_boku = elemJakobian.dlugoscBoku(side) / 2.0;

				for (int ip = 0; ip < 2; ip++) {
					double ksi, eta;
					double weight = weights2[ip];

					if (side == 0) {
						eta = -1.0;
						ksi = points2[ip];
					} else if (side == 1) {
						ksi = 1.0;
						eta = points2[ip];
					} else if (side == 2) {
						eta = 1.0;
						ksi = -points2[ip];
					} else if (side == 3) {
						ksi = -1.0;
						eta = -points2[ip];
					}

					double N[4];
					elemJakobian.pochodneFunkcjiKsztaltu(ksi, eta, N, nullptr, nullptr);

					for (int m = 0; m < 4; m++) {
						P_local[m] += N[m] * globalData.alfa * globalData.tot * weight * detJ_boku;
					}
				}
			}
		}

		cout << "Macierz Hbc dla elementu " << e + 1 << ":" << endl;
		for (const auto& row : elemJakobian.Hbc_local) {
			for (double val : row) {
				cout << val << " ";
			}
			cout << endl;
		}
		cout << endl;

		// dodajemy Hbc do H lokalnej
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				elemJakobian.H_local[i][j] += elemJakobian.Hbc_local[i][j];
			}
		}

		cout << "Wektor P_local dla elementu " << e + 1 << ":" << endl;
		for (int i = 0; i < 4; i++) {
			cout << P_local[i] << " ";
		}
		cout << endl << endl;

		// dodajemy P lokalne do P globalnego
		for (int i = 0; i < 4; i++) {
			int global_i = currentElement.ID[i] - 1;
			P_global[global_i] += P_local[i];
		}

		cout << "Macierz H dla elementu " << e + 1 << ":" << endl;
		for (const auto& row : elemJakobian.H_local) {
			for (double val : row) {
				cout << val << " ";
			}
			cout << endl;
		}
		cout << endl;

		// lab 6 - obliczanie macierzy H globalnej
		for (int i = 0; i < 4; i++) {
			int global_i = currentElement.ID[i] - 1;
			for (int j = 0; j < 4; j++) {
				int global_j = currentElement.ID[j] - 1;
				H_global[global_i][global_j] += elemJakobian.H_local[i][j];
			}
		}
		// end lab 6


		// lab 10/11 - macierz C globalna
		obliczMacierzC(elemJakobian, points2, weights2, 2,
					   globalData.density, globalData.specificHeat);

		cout << "Macierz C (lokalna) dla elementu " << e + 1 << ":" << endl;
		for (auto &row : elemJakobian.C_local) {
			for (double val : row) {
				cout << val << " ";
			}
			cout << endl;
		}
		cout << endl;

		// dodajemy C lokalne do C globalnej
		for (int i = 0; i < 4; i++) {
			int global_i = currentElement.ID[i] - 1;
			for (int j = 0; j < 4; j++) {
				int global_j = currentElement.ID[j] - 1;
				C_global[global_i][global_j] += elemJakobian.C_local[i][j];
			}
		}
	}

	// lab 5
	// nadrobienie kodu

	// lab 6
	// stworzenie macierzy H globalnej

	cout << "Wektor P globalny:" << endl;
	for (int i = 0; i < grid.nN; i++) {
		cout << P_global[i] << endl;
	}
	cout << endl;

	cout << "Macierz H globalna:" << endl;
	for (const auto& row : H_global) {
		for (double val : row) {
			cout << val << " ";
		}
		cout << endl;
	}

	// lab 10/11
	cout << endl << "Macierz C globalna:" << endl;
	for (auto &row : C_global) {
		for (double val : row) {
			cout << val << " ";
		}
		cout << endl;
	}
	cout << endl;

	// lab 10
	// obliczanie tempratury rozwiazujac uklad rownan
	vector<double> t(grid.nN, 0.0);

	t = eliminacjaGaussa(H_global, P_global);

	cout << endl;

	cout << "Temperatury w węzłach:" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << "T[" << i + 1 << "] = " << t[i] << endl;
	}

	// lab 11/12
	// final program code
	// solver dla zagadnienia nieustalonego wielokrokowego (transient)
	// wykorzystujemy H_global, C_global, P_global oraz parametry: simulationTime, simulationStepTime

	// liczba kroków czasowych
    int nSteps = (int) (globalData.simulationTime / globalData.simulationStepTime);

    // wektor temperatur T początkowych
    vector<double> T(grid.nN, globalData.initialTemp);

    // A = ( C_global/dt + H_global )
    // b = ( C_global/dt * T(n) + P_global )
    vector<vector<double>> A(grid.nN, vector<double>(grid.nN, 0.0));
    vector<double>         b(grid.nN, 0.0);

    double dt = globalData.simulationStepTime;

    cout << endl << "---- lab 12 - final program ----" << endl << endl;

    for(int step = 1; step <= nSteps; step++) {
        double currentTime = step * dt;
        cout << "Krok czasowy: " << step << ", czas = " << currentTime << endl;

        // A = C_global/dt + H_global
        for(int i = 0; i < grid.nN; i++){
            for(int j = 0; j < grid.nN; j++){
                A[i][j] = H_global[i][j] + (C_global[i][j] / dt);
            }
        }

        // b = (C_global/dt)*T(n) + P_global
        //   b[i] = sum( C_global[i][j]/dt * T[j] ) + P_global[i]
        for(int i = 0; i < grid.nN; i++){
            double sumCT = 0.0;
            for(int j = 0; j < grid.nN; j++){
                sumCT += (C_global[i][j] / dt) * T[j];
            }
            b[i] = sumCT + P_global[i];
        }

        // uklad A * T(n+1) = b - rozwiazujemy metoda eliminacji Gaussa
        vector<double> Tnew = eliminacjaGaussa(A, b);

        T = Tnew;

        double minT = T[0], maxT = T[0];
        for (int i = 1; i < grid.nN; i++) {
            if (T[i] < minT) minT = T[i];
            if (T[i] > maxT) maxT = T[i];
        }
        cout << "   -> Min T = " << minT << ", Max T = " << maxT << endl;
    }

    cout << endl << "Temperatury w wezlach:" << endl;
    for (int i = 0; i < T.size(); i++) {
        cout << "T[" << i+1 << "] = " << T[i] << endl;
    }

    return 0;
}

// funcs
int getDataFromFile(const string &fileName, GlobalData &globalData, Grid &grid) {
	ifstream inputFile(fileName);

	if (!inputFile.is_open()) {
		cerr << "Blad podczas otwierania pliku." << endl;
		return 0;
	}

	string line;

	while (getline(inputFile, line)) {
		istringstream iss(line);

		string keyword;
		iss >> keyword;

		if (keyword == "SimulationTime") {
			iss >> globalData.simulationTime;
		} else if (keyword == "SimulationStepTime") {
			iss >> globalData.simulationStepTime;
		} else if (keyword == "Conductivity") {
			iss >> globalData.conductivity;
		} else if (keyword == "Alfa") {
			iss >> globalData.alfa;
		} else if (keyword == "Tot") {
			iss >> globalData.tot;
		} else if (keyword == "InitialTemp") {
			iss >> globalData.initialTemp;
		} else if (keyword == "Density") {
			iss >> globalData.density;
		} else if (keyword == "SpecificHeat") {
			iss >> globalData.specificHeat;
		} else if (keyword == "Nodes") {
			string x;

			iss >> x >> grid.nN;
		} else if (keyword == "Elements") {
			string x;
			iss >> x >> grid.nE;
		} else if (keyword == "*Node") {
			string nvm;

			grid.nodes = new Node[grid.nN];

			for (int i = 0; i < grid.nN; i++) {
				string x, y;

				inputFile >> nvm >> x >> y;

				x.pop_back();
				y.pop_back();

				grid.nodes[i].x = stod(x);
				grid.nodes[i].y = stod(y);
			}
		} else if (keyword == "*Element,") {
			string x;

			grid.elements = new Element[grid.nE];

			for (int i = 0; i < grid.nE; i++) {
				string one, two, three, four;

				inputFile >> x >> one >> two >> three >> four;

				one.pop_back();
				two.pop_back();
				three.pop_back();

				grid.elements[i].ID[0] = stoi(one);
				grid.elements[i].ID[1] = stoi(two);
				grid.elements[i].ID[2] = stoi(three);
				grid.elements[i].ID[3] = stoi(four);
			}
		}
	}

	// do macierzy Hbc - lab 7
	// oblicz min oraz max z koordynatow
	double minX = grid.nodes[0].x, maxX = grid.nodes[0].x;
	double minY = grid.nodes[0].y, maxY = grid.nodes[0].y;

	for (int i = 1; i < grid.nN; i++) {
		if (grid.nodes[i].x < minX) minX = grid.nodes[i].x;
		if (grid.nodes[i].x > maxX) maxX = grid.nodes[i].x;
		if (grid.nodes[i].y < minY) minY = grid.nodes[i].y;
		if (grid.nodes[i].y > maxY) maxY = grid.nodes[i].y;
	}

	// sprawdz czy jest na granicy i ustaw flage
	for (int i = 0; i < grid.nN; i++) {
		if (grid.nodes[i].x == minX || grid.nodes[i].x == maxX || grid.nodes[i].y == minY || grid.nodes[i].y == maxY) {
			grid.nodes[i].isBoundary = true;
		}
	}

	return 1;
}

void displayGlobalData(const GlobalData &globalData) {
	cout << "SimulationTime: " << globalData.simulationTime << endl;
	cout << "SimulationStepTime: " << globalData.simulationStepTime << endl;
	cout << "Conductivity: " << globalData.conductivity << endl;
	cout << "Alfa: " << globalData.alfa << endl;
	cout << "Tot: " << globalData.tot << endl;
	cout << "InitialTemp: " << globalData.initialTemp << endl;
	cout << "Density: " << globalData.density << endl;
	cout << "SpecificHeat: " << globalData.specificHeat << endl;
}

void displayGrid(Grid grid) {
	cout << "Wezly: " << endl;
	for (int i = 0; i < grid.nN; i++) {
		cout << "x: " << grid.nodes[i].x << " " << "y: " << grid.nodes[i].y << endl;
	}

	cout << endl;

	cout << "Elementy: " << endl;
	for (int i = 0; i < grid.nE; i++) {
		cout << grid.elements[i].ID[0] << " " << grid.elements[i].ID[1] << " " << grid.elements[i].ID[2] << " " << grid.elements[i].ID[3] << endl;
	}
}
