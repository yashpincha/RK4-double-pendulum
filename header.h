#ifndef Header_h 
#define Header_h

class pendulum{
private:
    double t;
    int num ;
    // position
    double x,dx,ddx;
    double y,dy,ddy;
    // angles
    double theta,omega,domega ;
    // system params (mass and length)
    double length;
    double mass ;

    double energy ;
    double Ek;
    double Ep;
    public :
    
    pendulum(){};
    pendulum(int,double,double,double,double,double,double);
    ~pendulum(){};
      
    // init
    void init(int,double,double,double,double,double,double);
    void coordinates();
    
    //Energy
    void Energy(); 
    void ek();
    void ep();
    
    // Gets
    double GetEnergy() const {return energy;}
    double GetEk() const {return Ek;}
    double GetEp() const {return Ep;}
    int GetNum() const {return num;}
    double GetT() const{return t;}
    double GetX() const{return x;}
    double GetY()const {return y;}
    double GetDx()const {return dx;}
    double GetDdx()const{return ddx;}
    double GetDy()const{return dy;}
    double GetDdy() const{return ddy ;}
    double GetTheta()const{return theta ;}
    double GetOmega()const{return omega ;}
    double GetDomega()const{return domega ;}
    double GetLength() const {return length;}
    double GetMass() const {return mass ;}

    // Sets
    void SetEnergy(double a){energy=a;}
    void SetEk(double a){Ek=a;}
    void SetEp(double a){Ep=a;}
    void SetNum(double a){num=a;}
    void SetT(double a){t=a;}
    void SetX(double a){x=a;}
    void SetY(double a){y=a;}
    void SetDx(double a){dx=a;}
    void SetDy(double a){dy=a;}
    void SetDdx(double a){ddx=a;}
    void SetDdy(double a){ddy=a;}
    void SetTheta(double a){theta=a;}
    void SetOmega(double a){omega=a;}
    void SetDomega(double a){domega=a;}
    void SetLength(double a){length=a;}
    void SetMass(double a){mass=a;}
    
    // Increment += 
    void AddX(double a){x+=a;}
    void AddY(double a){y+=a;}
    void AddDx(double a){dx+=a;}
    void AddDy(double a){dy+=a;}
    void AddDdx(double a){ddx+=a;}
    void AddDdy(double a){ddy+=a;}
};

#endif
