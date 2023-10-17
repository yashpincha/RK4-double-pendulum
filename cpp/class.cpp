#include <iostream>
#include <cmath>

const double GRAVITY = -9.81;

class Pendulum {
private:
    int num;
    double t;
    double theta;
    double omega;
    double domega;
    double length;
    double mass;
    double x, y, dx, dy, ddx, ddy;
    double energy, Ep, Ek;

public:
    // Constructor
    Pendulum(int NUM, double TIME, double THETA, double OMEGA, double DOMEGA, double LENGTH, double MASS);

    // Member functions
    void init(int NUM, double TIME, double THETA, double OMEGA, double DOMEGA, double LENGTH, double MASS);
    void calculateCoordinates();
    void calculateEnergy();
    void calculateEp();
    void calculateEk();
};

// Constructor definition
Pendulum::Pendulum(int NUM, double TIME, double THETA, double OMEGA, double DOMEGA, double LENGTH, double MASS)
    : num(NUM), t(TIME), theta(THETA), omega(OMEGA), domega(DOMEGA), length(LENGTH), mass(MASS) {}

// Initialize function
void Pendulum::init(int NUM, double TIME, double THETA, double OMEGA, double DOMEGA, double LENGTH, double MASS) {
    num = NUM;
    t = TIME;
    theta = THETA;
    omega = OMEGA;
    domega = DOMEGA;
    length = LENGTH;
    mass = MASS;
}

// Calculate coordinates
void Pendulum::calculateCoordinates() {
    y = length * cos(theta);
    x = length * sin(theta);
    dy = -omega * length * sin(theta);
    dx = omega * length * cos(theta);
    ddy = -length * (domega * sin(theta) + pow(omega, 2) * cos(theta));
    ddx = length * (domega * cos(theta) - pow(omega, 2) * sin(theta));
}

// Calculate energy
void Pendulum::calculateEnergy() {
    calculateEp();
    calculateEk();
    energy = Ep + Ek;
}

// Calculate potential energy
void Pendulum::calculateEp() {
    Ep = -mass * y * GRAVITY;
}

// Calculate kinetic energy
void Pendulum::calculateEk() {
    double v2 = pow(dx, 2) + pow(dy, 2);
    Ek = 0.5 * mass * v2;
}

int main() {
    // Example usage
    Pendulum pendulumObj(1, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    pendulumObj.calculateCoordinates();
    pendulumObj.calculateEnergy();

    return 0;
}
