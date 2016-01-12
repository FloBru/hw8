#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;
void hamilton(double* p, double* q, double& H);
void step(double* p, double* q, double& H, double dt);


int main() {
    double H = -1./2;
    double tend = 20*M_PI;
    double dt = 0.05;
    double t = 0.0;
    double e = 0.6;
    double p[2];
    double q[2];
   
    int steps = tend/dt;

    //Startwerte
    p[0] = 1-e;
    p[1] = 0;
    q[0] = 0;
    q[1] = sqrt((1+e)/(1-e));

    ofstream out("Kepler05.txt");

    for (int i = 0 ; i <= steps; i++) {
        out << t << "\t" << p[0] << "\t" << p[1] << "\t" << q[0] << "\t" << q[1] << "\t" << H << endl;
        
        step(p, q, H, dt);
        
        t += dt;
    }
    out.close();

    return 0;

}
//Hamiltonfunktion
void hamilton(double* p, double* q, double& H) {
    H = 1.0/2.0*(pow(p[0],2) + pow(p[1],2) - 1.0/(sqrt(pow(q[0],2) + pow(q[1],2))));
}

//rechnet die neuen Werte von q,p und H aus 
void step(double* p, double* q, double& H, double dt) {
    p[0] += - dt * (q[0]/pow(q[0]*q[0] + q[1] * q[1], 3.0/2.0));	
    p[1] += - dt * (q[1]/pow(q[0]*q[0] + q[1] * q[1], 3.0/2.0));	
            
    q[0] += dt * p[0];	
    q[1] += dt * p[1];	

    hamilton(q,p,H);
}



