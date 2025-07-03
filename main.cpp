#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

float f1(float x) {
    return pow(x, 2) * pow(sin(x), 3);  
}

float f2(float x) {
    return exp(x * x) * (1 - x);  
}

float transform_x(float t, float a, float b) {
    return (b - a) * t / 2.0 + (a + b) / 2.0;
}

float jacobian(float a, float b) {
    return (b - a) / 2.0;
}

float gauss_legendre_single(float (*func)(float), float a, float b,
    vector<float>& nodes, vector<float>& weights) {
    float result = 0.0;
    float jac = jacobian(a, b);

    for (int i = 0; i < nodes.size(); i++) {
        float x = transform_x(nodes[i], a, b);
        result += weights[i] * func(x);
    }

    return result * jac;
}

float gauss_legendre_composite(float (*func)(float), float a, float b,
    vector<float>& nodes, vector<float>& weights, int n_subdivisions) {

    float h = (b - a) / n_subdivisions;  
    float total_result = 0.0;

    for (int i = 0; i < n_subdivisions; i++) {
        float a_sub = a + i * h;
        float b_sub = a + (i + 1) * h;
        total_result += gauss_legendre_single(func, a_sub, b_sub, nodes, weights);
    }

    return total_result;
}

// Funkcja obliczajaca blad wzgledny
float calculate_error(float numerical, float analytical) {
    if (analytical != 0) {
        return abs((numerical - analytical) / analytical) * 100.0;
    }
    return abs(numerical - analytical);
}

int main() {
    // Wezly i wagi dla roznych rzedow kwadratury

    // 2 wezly
    vector<float> nodes_2 = { -0.5773502692, 0.5773502692 };
    vector<float> weights_2 = { 1.0, 1.0 };

    // 3 wezly
    vector<float> nodes_3 = { -0.7745966692, 0.0, 0.7745966692 };
    vector<float> weights_3 = { 0.5555555556, 0.8888888889, 0.5555555556 };

    // 4 wezly
    vector<float> nodes_4 = { -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 };
    vector<float> weights_4 = { 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 };

    // Rozne liczby podzialow do testowania
    vector<int> subdivisions = { 1, 2, 4, 8, 16, 32 };

    cout << fixed << setprecision(8);

    // CALKA 1: f1(x) = x^2 * sin^3(x) w przedziale [1.0, 4.8]
    float a1 = 1.0, b1 = 4.8;
    float val1 = -10.9001;  // Wartosc analityczna dla calki 1
    cout << setprecision(2) << endl;
    cout << "============================================================================" << endl;
    cout << "CALKA 1: f(x) = x^2 * sin^3(x) w przedziale [" << a1 << ", " << b1 << "]" << endl;
    cout << "Wartosc analityczna: " << val1 << endl;
    cout << "============================================================================" << endl;
    cout << setprecision(6) << endl;
    // Nag³ówek tabeli
    cout << "Podzialy | 2 wezly        | 3 wezly        | 4 wezly        |" << endl;
    cout << "---------|----------------|----------------|----------------|" << endl;

    for (int n : subdivisions) {
        float result1_2 = gauss_legendre_composite(f1, a1, b1, nodes_2, weights_2, n);
        float result1_3 = gauss_legendre_composite(f1, a1, b1, nodes_3, weights_3, n);
        float result1_4 = gauss_legendre_composite(f1, a1, b1, nodes_4, weights_4, n);

        cout << setw(8) << n << " | "
            << setw(14) << result1_2 << " | "
            << setw(14) << result1_3 << " | "
            << setw(14) << result1_4 << " |" << endl;
    }

    cout << "\nBledy wzgledne dla Calki 1 [%]:" << endl;
    cout << "Podzialy | 2 wezly    | 3 wezly    | 4 wezly    |" << endl;
    cout << "---------|------------|------------|------------|" << endl;

    for (int n : subdivisions) {
        float result1_2 = gauss_legendre_composite(f1, a1, b1, nodes_2, weights_2, n);
        float result1_3 = gauss_legendre_composite(f1, a1, b1, nodes_3, weights_3, n);
        float result1_4 = gauss_legendre_composite(f1, a1, b1, nodes_4, weights_4, n);

        cout << setw(8) << n << " | "
            << setw(10) << calculate_error(result1_2, val1) << " | "
            << setw(10) << calculate_error(result1_3, val1) << " | "
            << setw(10) << calculate_error(result1_4, val1) << " |" << endl;
    }

    cout << "\n" << endl;

    // CALKA 2: f2(x) = e^(x^2) * (1-x) w przedziale [-1.5, 3.2]
    float a2 = -1.5, b2 = 3.2;
    float val2 = -9364.62;  // Wartosc analityczna dla calki 2
    cout << setprecision(2) << endl;
    cout << "============================================================================" << endl;
    cout << "CALKA 2: f(x) = e^(x^2) * (1-x) w przedziale [" << a2 << ", " << b2  << "]" << endl;
    cout << "Wartosc analityczna: " << val2 << endl;
    cout << "============================================================================" << endl;
    cout << setprecision(6) << endl;
    // Naglowek tabeli
    cout << "Podzialy | 2 wezly        | 3 wezly        | 4 wezly        |" << endl;
    cout << "---------|----------------|----------------|----------------|" << endl;

    for (int n : subdivisions) {
        float result2_2 = gauss_legendre_composite(f2, a2, b2, nodes_2, weights_2, n);
        float result2_3 = gauss_legendre_composite(f2, a2, b2, nodes_3, weights_3, n);
        float result2_4 = gauss_legendre_composite(f2, a2, b2, nodes_4, weights_4, n);

        cout << setw(8) << n << " | "
            << setw(14) << result2_2 << " | "
            << setw(14) << result2_3 << " | "
            << setw(14) << result2_4 << " |" << endl;
    }

    cout << "\nBledy wzgledne dla Calki 2 [%]:" << endl;
    cout << "Podzialy | 2 wezly    | 3 wezly    | 4 wezly    |" << endl;
    cout << "---------|------------|------------|------------|" << endl;

    for (int n : subdivisions) {
        float result2_2 = gauss_legendre_composite(f2, a2, b2, nodes_2, weights_2, n);
        float result2_3 = gauss_legendre_composite(f2, a2, b2, nodes_3, weights_3, n);
        float result2_4 = gauss_legendre_composite(f2, a2, b2, nodes_4, weights_4, n);

        cout << setw(8) << n << " | "
            << setw(10) << calculate_error(result2_2, val2) << " | "
            << setw(10) << calculate_error(result2_3, val2) << " | "
            << setw(10) << calculate_error(result2_4, val2) << " |" << endl;
    }

    // Analiza zbieznosci - porownanie kolejnych podzialow
    cout << "\n============================================================================" << endl;
    cout << "ANALIZA ZBIEZNOSCI - roznice miedzy kolejnymi podzialami" << endl;
    cout << "============================================================================" << endl;

    cout << "\nCALKA 1:" << endl;
    for (int i = 1; i < subdivisions.size(); i++) {
        int n_prev = subdivisions[i - 1];
        int n_curr = subdivisions[i];

        float result1_2_prev = gauss_legendre_composite(f1, a1, b1, nodes_2, weights_2, n_prev);
        float result1_2_curr = gauss_legendre_composite(f1, a1, b1, nodes_2, weights_2, n_curr);
        float result1_3_prev = gauss_legendre_composite(f1, a1, b1, nodes_3, weights_3, n_prev);
        float result1_3_curr = gauss_legendre_composite(f1, a1, b1, nodes_3, weights_3, n_curr);
        float result1_4_prev = gauss_legendre_composite(f1, a1, b1, nodes_4, weights_4, n_prev);
        float result1_4_curr = gauss_legendre_composite(f1, a1, b1, nodes_4, weights_4, n_curr);

        cout << "Podzialy " << n_prev << " -> " << n_curr << ":" << endl;
        cout << "  2 wezly: |" << result1_2_curr << " - " << result1_2_prev << "| = "
            << abs(result1_2_curr - result1_2_prev) << endl;
        cout << "  3 wezly: |" << result1_3_curr << " - " << result1_3_prev << "| = "
            << abs(result1_3_curr - result1_3_prev) << endl;
        cout << "  4 wezly: |" << result1_4_curr << " - " << result1_4_prev << "| = "
            << abs(result1_4_curr - result1_4_prev) << endl;
    }

    cout << "\nCALKA 2:" << endl;
    for (int i = 1; i < subdivisions.size(); i++) {
        int n_prev = subdivisions[i - 1];
        int n_curr = subdivisions[i];

        float result2_2_prev = gauss_legendre_composite(f2, a2, b2, nodes_2, weights_2, n_prev);
        float result2_2_curr = gauss_legendre_composite(f2, a2, b2, nodes_2, weights_2, n_curr);
        float result2_3_prev = gauss_legendre_composite(f2, a2, b2, nodes_3, weights_3, n_prev);
        float result2_3_curr = gauss_legendre_composite(f2, a2, b2, nodes_3, weights_3, n_curr);
        float result2_4_prev = gauss_legendre_composite(f2, a2, b2, nodes_4, weights_4, n_prev);
        float result2_4_curr = gauss_legendre_composite(f2, a2, b2, nodes_4, weights_4, n_curr);

        cout << "Podzialy " << n_prev << " -> " << n_curr << ":" << endl;
        cout << "  2 wezly: |" << result2_2_curr << " - " << result2_2_prev << "| = "
            << abs(result2_2_curr - result2_2_prev) << endl;
        cout << "  3 wezly: |" << result2_3_curr << " - " << result2_3_prev << "| = "
            << abs(result2_3_curr - result2_3_prev) << endl;
        cout << "  4 wezly: |" << result2_4_curr << " - " << result2_4_prev << "| = "
            << abs(result2_4_curr - result2_4_prev) << endl;
    }

    // Zapis do pliku dla wykresu zbieznosci
    ofstream file("convergence_composite.txt");
    if (file.is_open()) {
        file << "# Wyniki dla roznych podzialow - Calka 1: x^2 * sin^3(x)" << endl;
        file << "# Podzialy 2_wezly 3_wezly 4_wezly" << endl;
        for (int n : subdivisions) {
            float result1_2 = gauss_legendre_composite(f1, a1, b1, nodes_2, weights_2, n);
            float result1_3 = gauss_legendre_composite(f1, a1, b1, nodes_3, weights_3, n);
            float result1_4 = gauss_legendre_composite(f1, a1, b1, nodes_4, weights_4, n);
            file << n << " " << result1_2 << " " << result1_3 << " " << result1_4 << endl;
        }

        file << "# Wyniki dla roznych podzialow - Calka 2: e^(x^2) * (1-x)" << endl;
        file << "# Podzialy 2_wezly 3_wezly 4_wezly" << endl;
        for (int n : subdivisions) {
            float result2_2 = gauss_legendre_composite(f2, a2, b2, nodes_2, weights_2, n);
            float result2_3 = gauss_legendre_composite(f2, a2, b2, nodes_3, weights_3, n);
            float result2_4 = gauss_legendre_composite(f2, a2, b2, nodes_4, weights_4, n);
            file << n << " " << result2_2 << " " << result2_3 << " " << result2_4 << endl;
        }
        file.close();
        cout << "\nDane zapisane do pliku convergence_composite.txt" << endl;
    }

    return 0;
}