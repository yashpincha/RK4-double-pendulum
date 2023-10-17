#include <iostream>
#include <cmath>
#include <fstream>
#include "Header.h"

const double GRAVITY = -9.81;

void normalizeAngle(double &theta) {
    while (theta < -M_PI || theta > M_PI) {
        if (theta < -M_PI) {
            theta = 2 * M_PI + theta;
        }
        if (theta > M_PI) {
            theta = theta - 2 * M_PI;
        }
    }
}

double calculateFirstEquation(double omega1, double omega2, double theta1, double theta2,
                               double l1, double l2, double m1, double m2) {
    double dtheta = theta2 - theta1;
    double numerator = pow(omega1, 2) * m2 * l1 * cos(dtheta) * sin(dtheta) +
                       pow(omega2, 2) * m2 * l2 * sin(dtheta) -
                       (m1 + m2) * GRAVITY * sin(theta1) +
                       m2 * cos(dtheta) * GRAVITY * sin(theta2);
    double denominator = (m1 + m2) * l2 - m2 * l2 * pow(cos(dtheta), 2);
    return numerator / denominator;
}

double calculateSecondEquation(double omega1, double omega2, double theta1, double theta2,
                                double l1, double l2, double m1, double m2) {
    double dtheta = theta2 - theta1;
    double numerator = -pow(omega2, 2) * m2 * l2 * cos(dtheta) * sin(dtheta) +
                       (m1 + m2) * (GRAVITY * sin(theta1) * cos(dtheta) -
                                    l1 * pow(omega1, 2) * sin(dtheta) - GRAVITY * sin(theta2));
    double denominator = (m1 + m2) * l2 - m2 * l2 * pow(cos(dtheta), 2);
    return numerator / denominator;
}

double calculateG1(double omega) {
    return omega;
}

double calculateG2(double omega) {
    return omega;
}

void updateCoordinates(pendulum &pendulum1, pendulum &pendulum2) {
    pendulum1.coordinates();
    pendulum2.coordinates();
    pendulum2.AddX(pendulum1.GetX());
    pendulum2.AddY(pendulum1.GetY());
    pendulum2.AddDx(pendulum1.GetDx());
    pendulum2.AddDy(pendulum1.GetDy());
    pendulum2.AddDdx(pendulum1.GetDdx());
    pendulum2.AddDdy(pendulum1.GetDdy());

    pendulum1.ek();
    pendulum1.ep();
    pendulum2.ek();
    pendulum2.ep();
}

void rungeKutta(pendulum &pendulum1, pendulum &pendulum2, double h) {
    double t = pendulum1.GetT();
    double omega1 = pendulum1.GetOmega();
    double omega2 = pendulum2.GetOmega();
    double theta1 = pendulum1.GetTheta();
    double theta2 = pendulum2.GetTheta();
    double length1 = pendulum1.GetLength();
    double length2 = pendulum2.GetLength();
    double mass1 = pendulum1.GetMass();
    double mass2 = pendulum2.GetMass();

    double k1 = calculateFirstEquation(omega1, omega2, theta1, theta2, length1, length2, mass1, mass2);
    double j1 = calculateSecondEquation(omega1, omega2, theta1, theta2, length1, length2, mass1, mass2);
    double b1 = calculateG1(omega1);
    double p1 = calculateG2(omega2);

    double k2 = calculateFirstEquation(omega1 + h / 2. * k1, omega2 + h / 2. * j1,
                                        theta1 + h / 2. * b1, theta2 + h / 2. * p1,
                                        length1, length2, mass1, mass2);
    double j2 = calculateSecondEquation(omega1 + h / 2. * k1, omega2 + h / 2. * j1,
                                         theta1 + h / 2. * b1, theta2 + h / 2. * p1,
                                         length1, length2, mass1, mass2);
    double b2 = calculateG1(omega1 + h / 2. * b1);
    double p2 = calculateG2(omega2 + h / 2. * p1);

    double k3 = calculateFirstEquation(omega1 + h / 2. * k2, omega2 + h / 2. * j2,
                                        theta1 + h / 2. * b2, theta2 + h / 2. * p2,
                                        length1, length2, mass1, mass2);
    double j3 = calculateSecondEquation(omega1 + h / 2. * k2, omega2 + h / 2. * j2,
                                         theta1 + h / 2. * b2, theta2 + h / 2. * p2,
                                         length1, length2, mass1, mass2);
    double b3 = calculateG1(omega1 + h / 2. * b2);
    double p3 = calculateG2(omega2 + h / 2. * p2);

    double k4 = calculateFirstEquation(omega1 + h * k3, omega2 + h * j3,
                                        theta1 + h * b3, theta2 + h * p3,
                                        length1, length2, mass1, mass2);
    double j4 = calculateSecondEquation(omega1 + h * k3, omega2 + h * j3,
                                         theta1 + h * b3, theta2 + h * p3,
                                         length1, length2, mass1, mass2);
    double b4 = calculateG1(omega1 + h * b3);
    double p4 = calculateG2(omega2 + h * p3);

    double newOmega1 = omega1 + h / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
    double newOmega2 = omega2 + h / 6. * (j1 + 2 * j2 + 2 * j3 + j4);
    double newTheta1 = theta1 + h / 6. * (b1 + 2 * b2 + 2 * b3 + b4);
    double newTheta2 = theta2 + h / 6. * (p1 + 2 * p2 + 2 * p3 + p4);

    double newDomega1 = calculateFirstEquation(newOmega1, newOmega2, newTheta1, newTheta2, mass1, mass2, length1, length2);
    double newDomega2 = calculateSecondEquation(newOmega1, newOmega2, newTheta1, newTheta2, mass1, mass2, length1, length2);

    normalizeAngle(newTheta1);
    normalizeAngle(newTheta2);

    pendulum1.SetOmega(newOmega1);
    pendulum2.SetOmega(newOmega2);
    pendulum1.SetTheta(newTheta1);
    pendulum2.SetTheta(newTheta2);
    pendulum1.SetT(t + h);
    pendulum2.SetT(t + h);
    pendulum1.SetDomega(newDomega1);
    pendulum2.SetDomega(newDomega2);
    updateCoordinates(pendulum1, pendulum2);
}

int main() {
    pendulum pendulum1;
    pendulum pendulum2;
    ifstream initialConditions("initialisation.txt");
    ofstream coordinates("coordinates");
    ofstream altCoord("coordinatescolon");
    ofstream thetaValues("theta");
    ofstream energyValues("energy");
    ofstream phaseValues("phase");
    int num1, num2;
    double length1, length2, mass1, mass2, theta1, theta2, omega1, omega2, domega1, domega2;
    double tMin, tMax, h;

    initialConditions >> num1 >> length1 >> mass1 >> theta1 >> omega1 >> domega1;
    initialConditions >> num2 >> length2 >> mass2 >> theta2 >> omega2 >> domega2;
    initialConditions >> tMin >> tMax >> h;

    pendulum1.init(num1, tMin, theta1, omega1, domega1, length1, mass1);
    pendulum2.init(num2, tMin, theta2, omega2, domega2, length2, mass2);
    updateCoordinates(pendulum1, pendulum2);

    altCoord << "#" << "time" << " " << "X1" << " " << "Y1" << " " << "X2" << " " << "Y2" << endl;
    thetaValues << "#" << "time" << " " << "theta1" << " " << "theta2" << endl;
    phaseValues << "#" << "time" << " " << "theta1" << " " << "omega1" << " " << "theta2" << " " << "omega2" << endl;
    energyValues << "#" << "time" << " " << "Etot" << " " << "Ek" << " " << "Ep" << endl;

    int i = 0;
    double t = tMin;
    while (t < tMax) {
        if (i % 100 == 0) {
            altCoord << pendulum1.GetT() << " " << pendulum1.GetX() << " " << pendulum1.GetY() << " "
                     << pendulum2.GetX() << " " << pendulum2.GetY() << endl;

            coordinates << pendulum1.GetT() << " " << "0" << " " << "0" << endl;
            coordinates << pendulum1.GetT() << " " << pendulum1.GetX() << " " << pendulum1.GetY() << endl;
            coordinates << pendulum2.GetT() << " " << pendulum2.GetX() << " " << pendulum2.GetY() << endl;
            coordinates << " " << endl;
            coordinates << " " << endl;

            thetaValues << pendulum1.GetT() << " " << pendulum1.GetTheta() << " " << pendulum2.GetTheta() << endl;

            phaseValues << pendulum1.GetT() << " " << pendulum1.GetTheta() << " " << pendulum1.GetOmega() << " "
                        << pendulum2.GetTheta() << " " << pendulum2.GetOmega() << endl;

            double Ek = pendulum1.GetEk() + pendulum2.GetEk();
            double Ep = pendulum1.GetEp() + pendulum2.GetEp();
            double E = pendulum1.GetEnergy() + pendulum2.GetEnergy();

            energyValues << pendulum1.GetT() << " " << E << " " << Ek << " " << Ep << endl;
        }

        rungeKutta(pendulum1, pendulum2, h);
        t += h;
        i += 1;
    }

    ofstream traceEnergy("traceenergy.gnu");
    traceEnergy << "set title \"Energy traced as a function of time\"" << endl;
    traceEnergy << "set xlabel \"time t\"" << endl;
    traceEnergy << "set ylabel \"energy E\"" << endl;
    traceEnergy << "set style fill transparent solid 0.2 noborder" << endl;
    traceEnergy << "plot \"energy\" using 1:2 title 'Etot' with lines,\"energy\" using 1:3 title 'Ek' with lines,\"energy\" using 1:4 title 'Ep' with lines" << endl;
    traceEnergy << "pause -1" << endl;
    traceEnergy << "set terminal postscript eps size 3.5,2.62 enhanced color" << endl;
    traceEnergy << "set output \"energy.eps\"" << endl;
    traceEnergy << "replot" << endl;

    system("gnuplot traceenergy.gnu");

    ofstream multiplot("multiplot.gnu");
    multiplot << "set multiplot layout 2,2 rowsfirst" << endl;
    multiplot << "set title \"omega1 as a function of theta1\"" << endl;
    multiplot << "set xlabel \"theta1\"" << endl;
    multiplot << "set xrange [-3.5:3.5]" << endl;
    multiplot << "set ylabel \"omega1\"" << endl;
    multiplot << "set yrange [-14:14]" << endl;
    multiplot << "set style fill transparent solid 0.2 noborder" << endl;
    multiplot << "plot \"phase\" using 2:3 title 'theta1point=f(theta1)' with lines" << endl;

    multiplot << "set title \"omega2 as a function of theta2\"" << endl;
    multiplot << "set xlabel \"theta2\"" << endl;
    multiplot << "set xrange [-3.5:3.5]" << endl;
    multiplot << "set ylabel \"omega2t\"" << endl;
    multiplot << "set yrange [-14:14]" << endl;
    multiplot << "set style fill transparent solid 0.2 noborder" << endl;
    multiplot << "plot \"phase\" using 4:5 title 'omega2=f(theta2)' with lines" << endl;

    multiplot << "set title \"Coordinates as a function of time\"" << endl;
    multiplot << "set xlabel \"time t\"" << endl;
    multiplot << "set xrange [0:" << tMax << "]" << endl;
    multiplot << "set ylabel \"Coordinates\"" << endl;
    multiplot << "set yrange [-1.1:1.1]" << endl;
    multiplot << "set style fill transparent solid 0.2 noborder" << endl;
    multiplot << "plot \"coordinatescolon\" using 1:2 title 'x1' with lines,\"coordinatescolon\" using 1:3 title 'y1' with lines" << endl;

    multiplot << "set title \"Coordinates as a function of time\"" << endl;
    multiplot << "set xlabel \"time t\"" << endl;
    multiplot << "set xrange [0:" << tMax << "]" << endl;
    multiplot << "set ylabel \"Coordinates\"" << endl;
    multiplot << "set yrange [-2.1:2.1]" << endl;
    multiplot << "set style fill transparent solid 0.2 noborder" << endl;
    multiplot << "plot \"coordinatescolon\" using 1:4 title 'x2' with lines,\"coordinatescolon\" using 1:5 title 'y2' with lines" << endl;

    multiplot << "pause -1" << endl;

    system("gnuplot multiplot.gnu");

    ofstream animation("animationpendulum.gnu");
    animation << "set title \"Animation of the pendulum\"" << endl;
    animation << "set terminal gif animate delay " << h * 10000 << endl;
    animation << "set output \"pendulumanimated.gif\"" << endl;
    animation << "set xrange [" << -(length1 + length2) << ":" << (length1 + length2) << "]" << endl;
    animation << "set yrange [" << -(length1 + length2) << ":" << (length1 + length2) << "]" << endl;
    animation << "set style line 11 lc rgb '#808080' lt 1" << endl;
    animation << "set border 3 back ls 11" << endl;
    animation << "set tics nomirror" << endl;
    animation << "set pointintervalbox 3" << endl;
    animation << "do for [i=1:999] {plot \"coordinates\" index i using 2:3 title 'pendulum' w linespoints ls 1,\"coordinatescolon\" every ::i-25::i using 4:5 title 'trajectory' with lines lt rgb \"violet\"}" << endl;

    system("gnuplot animationpendulum.gnu");

    return 0;
}
