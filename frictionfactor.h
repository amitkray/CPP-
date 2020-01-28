#ifndef FRICTIONFACTOR_H
#define FRICTIONFACTOR_H
#include <iostream>
#include <math.h>

using namespace std;

class frictionfactor {
  double R, e, d;
  public:
    double frict(double R,double e,double d){
        double f;//Friction  Factor
        f = 0.25/(pow((log10(e/3.7*d+5.74/(pow(R,0.9)))),2));
        return f;

    }
};
/*
int main(){
    frictionfactor fric;
    double values = fric.frict(250000,0.000015,0.225);
    cout<<"Values is "<<values;
}
*/

#endif // FRICTIONFACTOR_H
