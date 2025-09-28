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


double beta = 3.5;   // taxa de transmissão
double gamma = 0.6;  // taxa de recuperação
double mi = 0.1;     // morte natural
double d = 0.3;       // morte por doença
int N;            // população total

State derivate(const State& y) {
    State dydt;
    dydt.S = (mi * N) + (d * y.I) - (mi * y.S) - (beta * y.S * y.I / (double)N);
    dydt.I = (beta * y.S * y.I / (double)N) - (gamma * y.I) - (mi * y.I) - (d * y.I);
    dydt.R = (gamma * y.I) - (mi * y.R);
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
    N = 40000;
    y0.I = 200;    // infectados
    y0.S = N - y0.I;   // suscetíveis
    y0.R = 0;     // recuperados
    y0.D = 0;     // mortos


    std::cout << "gama: " << gamma << "\n beta: " << beta << "\n\n";


    // parâmetros da simulação
    int T = 100;   // tempo total
    int h = 1;     // passo de integração
    int steps = (T/h);

    std::ofstream file("sir_rk4.csv");

    State y = y0;
    int t = 0;

    for (int i = 0; i <= steps; i++) {
        file << std::fixed 
             << t << "," << float(y.S/N) << "," << float(y.I/N) << "," << float(y.R/N) << "\n";

        if((y.S + y.I + y.R) > N-1 && N+1 < (y.S + y.I + y.R)){
            std::cout << "erro no tempo: " << i << "\n";
            std::cout << y.S+y.I+y.R << "," << int(y.S) << "," << int(y.I) << "," << int(y.R) << "\n";
        }

        y = rk4(y, h);
        t += h;
    }

    file.close();

    return 0;
}
