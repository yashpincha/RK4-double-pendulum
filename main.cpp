// Yash Pincha 09/2022

#include<stdlib.h>
#include "Header.h" 
#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

using namespace std ;

void angle(double & theta)
{
    while (theta <-M_PI|| theta>M_PI){
        if(theta<-M_PI){theta=2*M_PI+theta;} 
        if(theta> M_PI){theta=theta-2*M_PI;}
                                     }
}

//omega1, omega2 = dtheta1, dtheta 2. domega1, domega2 = ddottheta1, ddothteta2
double f1(double omega1 , double omega2, double theta1, double theta2, double l1, double l2, double m1, double m2)
{
    double dtheta=theta2-theta1;
    double g=-9.81; // acc.
    double domega1=(pow(omega1,2)*m2*l1*cos(dtheta)*sin(dtheta)+pow(omega2,2)*m2*l2*sin(dtheta)-
    (m1+m2)*g*sin(theta1)+m2*cos(dtheta)*g*sin(theta2))/((m1+m2)*l2-m2*l2*pow(cos(dtheta),2)); // EOM_1
    return domega1;
    
}

// acceleration for num = 2
double f2(double omega1 , double omega2, double theta1, double theta2, double l1, double l2, double m1, double m2)
{
    double dtheta=theta2-theta1;
    double g=-9.81;
    double domega2=(-pow(omega2,2)*m2*l2*cos(dtheta)*sin(dtheta)+(m1+m2)*(g*sin(theta1)*cos(dtheta)-
    l1*pow(omega1,2)*sin(dtheta)-g*sin(theta2)))/((m1+m2)*l2-m2*l2*pow(cos(dtheta),2)); // EOM_2
    return domega2;
}

double g1(double omega)
{
    return omega;
}

double g2(double omega)
{
    return omega ;
}

// updating coords
void coordinate(pendulum & pend_1 , pendulum & pend_2)
{
    pend_1.coordinates();
    pend_2.coordinates();
    pend_2.AddX(pend_1.GetX());
    pend_2.AddY(pend_1.GetY());
    pend_2.AddDx(pend_1.GetDx());
    pend_2.AddDy(pend_1.GetDy());
    pend_2.AddDdx(pend_1.GetDdx());
    pend_2.AddDdy(pend_1.GetDdy());
    
    pend_1.ek();
    pend_1.ep();
    pend_2.ek();
    pend_2.ep();
}

void rk(pendulum & pend_1, pendulum & pend_2, double h)

{
    double t = pend_1.GetT();
    double const omega1=pend_1.GetOmega();
    double const omega2=pend_2.GetOmega();
    double const theta1=pend_1.GetTheta();
    double const theta2=pend_2.GetTheta();
    double const length1 =pend_1.GetLength();
    double const length2=pend_2.GetLength();
    double const mass1=pend_1.GetMass();
    double const mass2=pend_2.GetMass();
    
    double k1=f1(omega1,omega2,theta1,theta2,length1,length2,mass1,mass2);
    double j1=f2(omega1,omega2,theta1,theta2,length1,length2,mass1,mass2);
    double b1=g1(omega1);
    double p1=g2(omega2);
    
    double k2=f1(omega1+h/2.*k1,omega2+h/2.*j1,theta1+h/2*b1,theta2+h/2.*p1,length1,length2,mass1,mass2);
    double j2=f2(omega1+h/2.*k1,omega2+h/2.*j1,theta1+h/2.*b1,theta2+h/2.*p1,length1,length2,mass1,mass2);
    double b2=g1(omega1+h/2.*b1);
    double p2=g2(omega2+h/2.*p1);
    
    double k3=f1(omega1+h/2.*k2,omega2+h/2.*j2,theta1+h/2.*b2,theta2+h/2.*p2,length1,length2,mass1,mass2);
    double j3=f2(omega1+h/2.*k2,omega2+h/2.*j2,theta1+h/2.*b2,theta2+h/2.*p2,length1,length2,mass1,mass2);
    double b3=g1(omega1+h/2.*b2);
    double p3=g2(omega2+h/2.*p2);
    
    double k4=f1(omega1+h*k3,omega2+h*j3,theta1+h*b3,theta2+h*p3,length1,length2,mass1,mass2);
    double j4=f2(omega1+h*k3,omega2+h*j3,theta1+h*b3,theta2+h*p3,length1,length2,mass1,mass2);
    double b4=g1(omega1+h*b3);
    double p4=g2(omega2+h*p3);
    
    //RK4

    double Omega1=omega1+h/6.*(k1+2*k2+2*k3+k4);
    double Omega2=omega2+h/6.*(j1+2*j2+2*j3+j4);
    double Theta1 = theta1+h/6.*(b1+2*b2+2*b3+b4);
    double Theta2 =theta2+h/6.*(p1+2*p2+2*p3+p4);
    
    double Domega1=f1(Omega1,Omega2,Theta1,Theta2,mass1,mass2,length1,length2);
    double Domega2=f2(Omega1,Omega2,Theta1,Theta2,mass1,mass2,length1,length2);
    angle(Theta1);
    angle(Theta2);
    
    pend_1.SetOmega(Omega1);
    pend_2.SetOmega(Omega2);
    pend_1.SetTheta(Theta1);
    pend_2.SetTheta(Theta2);
    pend_1.SetT(t+h);
    pend_2.SetT(t+h);
    pend_1.SetDomega(Domega1);
    pend_2.SetDomega(Domega2);
    coordinate(pend_1,pend_2);
    
}

/*
void smallangle(pendulum & pend_1, pendulum & pend_2, double theta0, double t)
{
    double g=-9.81;
    double  l1=pend_1.GetLength();
    double theta=theta0*cos(sqrt(g/l1)*t);
    double dtheta=-theta0*sqrt(g/l1)*sin(sqrt(g/l1)*t);
    pend_1.SetTheta(theta);
    pend_1.SetOmega(dtheta);
    coordinate(pend_1,pend_2);
}
*/

int main()
{
    pendulum pend_1;
    pendulum pend_2;
    ifstream initialcondition("initialisation.txt");
    ofstream coordinates("coordinates");
    ofstream altcoord("coordinatescolon");
    ofstream valeur_theta("theta");
    ofstream energy("energy");
    ofstream phase("phase");
    int num1 ,num2 ;
    double length1, length2 , mass1 , mass2 , theta1 , theta2 , omega1,omega2, domega1 ,domega2 ;
    double tmin , tmax , h;
    
    initialcondition>>num1>>length1>>mass1>>theta1>>omega1>>domega1;
    initialcondition>>num2>>length2>>mass2>>theta2>>omega2>>domega2;
    initialcondition>>tmin>>tmax>>h;
    
    pend_1.init(num1,tmin,theta1,omega1 ,domega1,length1,mass1);
    pend_2.init(num2,tmin,theta2,omega2, domega2, length2,mass2);
    coordinate(pend_1,pend_2);

    altcoord<<"#"<<"time"<<" "<<"X1"<<" "<<"Y1"<<" "<<"X2"<<" "<<"Y2"<<endl;
    valeur_theta<<"#"<<"time"<<" "<<"theta1"<<" "<<"theta2"<<endl;
    phase<<"#"<<"time"<<" "<<"theta1"<<" "<<"omega1"<<" "<<"theta2"<<" "<<"omega2"<<endl;
    energy<<"#"<<"time"<<" "<<"Etot"<<" "<<"Ek"<<" "<<"Ep"<<endl;
    
    int i=0;
    double t=tmin;
    while(t<tmax)
    {
        if (i%100==0){
            
            altcoord<<pend_1.GetT()<<" "<<pend_1.GetX()<<" "<<pend_1.GetY()<<" "<<pend_2.GetX()<<" "<<pend_2.GetY()<<endl;
            
            coordinates<<pend_1.GetT()<<" "<<"0"<<" "<<"0"<<endl; // origin
            coordinates<<pend_1.GetT()<<" "<<pend_1.GetX()<<" "<<pend_1.GetY()<<endl;
            coordinates<<pend_2.GetT()<<" "<<pend_2.GetX()<<" "<<pend_2.GetY()<<endl;
            coordinates<<" "<<endl;
            coordinates<<" "<<endl;
            
            valeur_theta<<pend_1.GetT()<<" "<<pend_1.GetTheta()<<" "<<pend_2.GetTheta()<<endl ; // theta1 and theta2
            
            phase<<pend_1.GetT()<<" "<<pend_1.GetTheta()<<" "<<pend_1.GetOmega()<<" "<<pend_2.GetTheta()<<" "<<pend_2.GetOmega()<<endl;
            
            // total energy
            double Ek= pend_1.GetEk()+ pend_2.GetEk();
            double Ep= pend_1.GetEp()+pend_2.GetEp();
            double E=  pend_1.GetEnergy()+pend_2.GetEnergy();
            
            energy<<pend_1.GetT()<<" "<<E<<" "<<Ek<<" "<<Ep<<endl;
        }
        
        de(pend_1, pend_2,h);
        t+=h;
        i+=1;
    }
    
    ofstream trac_energy("traceenergy.gnu");
    trac_energy<<"set title \"Energy traced as a function of time\""<<endl;
    trac_energy<<"set xlabel \"time t\""<<endl;
    trac_energy<<"set ylabel \"energy E\""<<endl;
    trac_energy<<"set style fill transparent solid 0.2 noborder"<<endl;
    trac_energy<<"plot \"energy\" using 1:2 title 'Etot' with lines,\"energy\" using 1:3 title 'Ek' with lines,\"energy\" using 1:4 title 'Ep' with lines"<<endl;
    trac_energy<< "pause -1" <<endl;
    trac_energy<<"set terminal postscript eps size 3.5,2.62 enhanced color"<<endl;
    trac_energy<<"set output \"energy.eps\""<<endl;
    trac_energy<<"replot"<<endl;
    
    system("gnuplot traceenergy.gnu");
    ofstream multiplot("multiplot.gnu");
    multiplot<<"set multiplot layout 2,2 rowsfirst"<<endl;
    multiplot<<"set title \"omega1 as a function of theta1\""<<endl;
    multiplot<<"set xlabel \"theta1\""<<endl;
    multiplot<<"set xrange [-3.5:3.5]"<<endl;
    multiplot<<"set ylabel \"omega1\""<<endl;
    multiplot<<"set yrange [-14:14]"<<endl;
    multiplot<<"set style fill transparent solid 0.2 noborder"<<endl;
    multiplot<<"plot \"phase\" using 2:3 title 'theta1point=f(theta1)' with lines"<<endl; // phase
    
    multiplot<<"set title \"omega2 as a function of theta2\""<<endl;
    multiplot<<"set xlabel \"theta2\""<<endl;
    multiplot<<"set xrange [-3.5:3.5]"<<endl;
    multiplot<<"set ylabel \"omega2t\""<<endl;
    multiplot<<"set yrange [-14:14]"<<endl;
    multiplot<<"set style fill transparent solid 0.2 noborder"<<endl;
    multiplot<<"plot \"phase\" using 4:5 title 'omega2=f(theta2)' with lines"<<endl; //phase
    
    multiplot<<"set title \"Coordinates as a function of time\""<<endl;
    multiplot<<"set xlabel \"time t\""<<endl;
    multiplot<<"set xrange [0:"<<tmax<<"]"<<endl;
    multiplot<<"set ylabel \"Coordinates\""<<endl;
    multiplot<<"set yrange [-1.1:1.1]"<<endl;
    multiplot<<"set style fill transparent solid 0.2 noborder"<<endl;
    multiplot<<"plot \"coordinatescolon\" using 1:2 title 'x1' with lines,\"coordinatescolon\" using 1:3 title 'y1' with lines"<<endl; 
    
    multiplot<<"set title \"Coordinates as a function of time\""<<endl;
    multiplot<<"set xlabel \"time t\""<<endl;
    multiplot<<"set xrange [0:"<<tmax<<"]"<<endl;
    multiplot<<"set ylabel \"Coordinates\""<<endl;
    multiplot<<"set yrange [-2.1:2.1]"<<endl;
    multiplot<<"set style fill transparent solid 0.2 noborder"<<endl;
    multiplot<<"plot \"coordinatescolon\" using 1:4 title 'x2' with lines,\"coordinatescolon\" using 1:5 title 'y2' with lines"<<endl; 
    
    multiplot<<"pause -1"<<endl; 
    
    system("gnuplot multiplot.gnu");
    
    ofstream animation("animationpendulum.gnu");
    animation<<"set title \"Animation of the pendulum\""<<endl;
    animation<<"set terminal gif animate delay "<<h*10000<<endl;
    animation<<"set output \"pendulumanimated.gif\"" <<endl;
    animation<<"set xrange ["<<-(length1+length2)<<":"<<(length1+length2)<<"]"<<endl;
    animation<<"set yrange ["<<-(length1+length2)<<":"<<(length1+length2)<<"]"<<endl;
    animation<<"set style line 11 lc rgb '#808080' lt 1"<<endl;
    animation<<"set border 3 back ls 11"<<endl;
    animation<<"set tics nomirror"<<endl;
    animation<<"set pointintervalbox 3"<<endl;
    animation<<"do for [i=1:999] {plot \"coordinates\" index i using 2:3 title 'pendulum' w linespoints ls 1,\"coordinatescolon\" every ::i-25::i using 4:5 title 'trajectoire' with lines lt rgb \"violet\"}"<<endl;

    system("gnuplot animationpendulum.gnu");
    
    return 0;
}
