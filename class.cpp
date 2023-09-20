#include <stdio.h>
#include"Header.h"
#include <iostream>
#include <cmath>
using namespace  std ;
double g=-9.81 ;

pendulum::pendulum(int NUM,double TIME,double THETA,double OMEGA,double DOMEGA,double LENGTH,double MASS)
{
    num=NUM;
    t=TIME;
    theta=THETA;
    omega=OMEGA;
    domega=DOMEGA;
    length=LENGTH;
    mass=MASS;
}

void pendulum::init(int NUM,double TIME,double THETA,double OMEGA,double DOMEGA,double LENGTH,double MASS)
{
    num=NUM;
    t=TIME;
    theta=THETA;
    omega=OMEGA;
    domega=DOMEGA;
    length=LENGTH;
    mass=MASS;
}

void pendulum::Energy()
{
    double g =-9.81;
    double Ep = -mass*y*g;
    double v2 = pow(dx,2)+pow(dy,2);
    double Ek =0.5*mass*v2;
    energy=Ep+Ek ;
}

void pendulum::ep()
{
    double g=-9.81;
    Ep=-mass*y*g;
}

void pendulum::ek()
{
    double v2 = pow(dx,2)+pow(dy,2);
    Ek=0.5*mass*v2;
}

void pendulum::coordinates()
{
    y=length*cos(theta);
    x=length*sin(theta);
    dy=-omega*length*sin(theta);
    dx=omega*length*cos(theta);
    ddy=-length*(domega*sin(theta)+pow(omega,2)*cos(theta));
    ddx=length*(domega*cos(theta)-pow(omega,2)*sin(theta));
}
