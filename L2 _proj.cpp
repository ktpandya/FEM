#include <iostream>
#include<vector>
#include <iomanip>
#include<math.h>  
#include "Eigen/Dense"
#include "Eigen/Eigen"

using namespace std;

//Code for printing vectors
void print(std::vector <double> a) {
   std::cout << "\n";
   std::cout << std::setprecision(3) << std::fixed;
   for(int i=0; i < a.size(); i++)
       {std::cout << a.at(i) << ' ';
        
        }
    std::cout << "\n \n";
}

//Code for printing matrix
void printm(std::vector<std::vector <double>> vect) {
   std::cout << "\n";
   std::cout << std::setprecision(3) << std::fixed;
   for (int i = 0; i < vect.size(); i++)
        {
            for (int j = 0; j < vect[i].size(); j++)
            {
                cout << vect[i][j] << " ";
            }    
            cout << endl;
        }
  
}

//Defining Functions
double f(double x)
{
    double func= cos(x);
    return func;
}

//Element Force vector
std::vector<double> ComputeElementForceVector(double Xa,double Xb)
{
    vector<double> frc;
    for (int i = 0; i < 2; i++)
    {
        frc.push_back((f(Xa*(1-i) + i * Xb ))*(Xb-Xa)/2);

    }
    return frc;

}

//Element Mass Matrix
std::vector<std::vector<double>> ComputeElementMassMatrix(double Xa,double Xb)
{
    vector<vector<double>> mass;
    mass = {{2*(Xb-Xa)/6,1*(Xb-Xa)/6},{1*(Xb-Xa)/6,2*(Xb-Xa)/6}};
    return mass;

}

int main()
{
    double a = 0;
    double b = 1;
    double n = 10;
    int nel = n - 1;
    vector<double> X_vect1;
    for (int i = 0; i < n+1;i++)
    {
       X_vect1.push_back(a + i/n);
    }
    
    return 0;
}
