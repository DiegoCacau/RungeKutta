#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

struct State {
    double S;
    double I;
    double R;
    double D;
};


double beta = 0.5;   // taxa de transmissão
double gamma = 0.3;  // taxa de recuperação
int N;            // população total

State derivate(const State& y) {
    State dydt;
    dydt.S = -beta * y.S * y.I / (double)N;
    dydt.I = (beta * y.S * y.I / (double)N) - (gamma * y.I);
    dydt.R = gamma * y.I;
    return dydt;
}

// Runge-Kutta 
State rk4(const State& y, int h) {
    State k1 = derivate(y);

    State y2 = { y.S + 0.5*h*k1.S, y.I + 0.5*h*k1.I, y.R + 0.5*h*k1.R };
    State k2 = derivate(y2);

    State y3 = { y.S + 0.5*h*k2.S, y.I + 0.5*h*k2.I, y.R + 0.5*h*k2.R };
    State k3 = derivate(y3);

    State y4 = { y.S + h*k3.S, y.I + h*k3.I, y.R + h*k3.R };
    State k4 = derivate(y4);

    State yn;
    yn.S = y.S + (h/6.0)*(k1.S + 2*k2.S + 2*k3.S + k4.S);
    yn.I = y.I + (h/6.0)*(k1.I + 2*k2.I + 2*k3.I + k4.I);
    yn.R = y.R + (h/6.0)*(k1.R + 2*k2.R + 2*k3.R + k4.R);
    return yn;
}

int main() {
    // Condições iniciais
    State y0;
    N = 90000;
    y0.I = 200;    // infectados
    y0.S = N - y0.I;   // suscetíveis
    y0.R = 0;     // recuperados
    y0.D = 0;     // mortos


    std::cout << "gama: " << gamma << "\n beta: " << beta << "\n\n";


    // parâmetros da simulação
    int T = 90;   // tempo total
    int h = 1;     // passo de integração
    int steps = (T/h);

    std::ofstream file("sir_rk4.csv");

    State y = y0;
    int t = 0;

    for (int i = 0; i <= steps; i++) {
        file << std::fixed 
             << t << "," << int(y.S) << "," << int(y.I) << "," << int(y.R) << "\n";

        y = rk4(y, h);
        t += h;
    }

    file.close();

    return 0;
}
