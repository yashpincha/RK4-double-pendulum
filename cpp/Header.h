#ifndef Header_h
#define Header_h

class Pendulum {
private:
    // Time and identification
    double t;
    int num;

    // Position variables
    double x, dx, ddx;
    double y, dy, ddy;

    // Angular variables
    double theta, omega, domega;

    // System parameters (mass and length)
    double length;
    double mass;

    // Energy variables
    double energy;
    double Ek;
    double Ep;

public:
    // Constructors and Destructor
    Pendulum();
    Pendulum(int num, double t, double x, double dx, double ddx, double y, double dy, double ddy,
              double theta, double omega, double domega, double length, double mass);
    ~Pendulum();

    // Initialization
    void init(int num, double t, double x, double dx, double ddx, double y, double dy, double ddy,
              double theta, double omega, double domega, double length, double mass);

    // Coordinate calculations
    void calculateCoordinates();

    // Energy calculations
    void calculateEnergy();
    void calculateEk();
    void calculateEp();

    // Getters
    double getEnergy() const;
    double getEk() const;
    double getEp() const;
    int getNum() const;
    double getT() const;
    double getX() const;
    double getY() const;
    double getDx() const;
    double getDdx() const;
    double getDy() const;
    double getDdy() const;
    double getTheta() const;
    double getOmega() const;
    double getDomega() const;
    double getLength() const;
    double getMass() const;

    // Setters
    void setEnergy(double value);
    void setEk(double value);
    void setEp(double value);
    void setNum(int value);
    void setT(double value);
    void setX(double value);
    void setY(double value);
    void setDx(double value);
    void setDy(double value);
    void setDdx(double value);
    void setDdy(double value);
    void setTheta(double value);
    void setOmega(double value);
    void setDomega(double value);
    void setLength(double value);
    void setMass(double value);

    // Increment operators
    void addToX(double value);
    void addToY(double value);
    void addToDx(double value);
    void addToDy(double value);
    void addToDdx(double value);
    void addToDdy(double value);
};

#endif  // Header_h
