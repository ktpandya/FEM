#include <iostream>
#include<vector>
#include <iomanip>
#include<math.h>  
#include "Eigen/Dense"
#include "Eigen/Eigen"
#include<fstream>
#include<iterator>
#include <Python.h>

using namespace std;
//using namespace Eigen;
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
    double func = cos(x);
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

//Global Mass Matrix
Eigen::SparseMatrix<double> GlobalMassMatrix(vector<double> X_vect)
{
   Eigen::SparseMatrix<double> A(X_vect.size(),X_vect.size());
   typedef Eigen::Triplet<double> M;
   std::vector<M> tripletList{};
   tripletList.reserve(3*X_vect.size());
   for(int i = 0; i< X_vect.size()-1;i++)
   {
    for (int k = -1; k<2; k++)
    { 
        if (i==i+k )
        {
            tripletList.push_back(M(i,i+k,ComputeElementMassMatrix(X_vect.at(i),X_vect.at(i+1))[0][0]));
            tripletList.push_back(M(i+1,i+k+1,ComputeElementMassMatrix(X_vect.at(i),X_vect.at(i+1))[1][1]));
        }
        else if (i + k <= X_vect.size() && i != i+k )
        {
            tripletList.push_back(M(i,i+k,ComputeElementMassMatrix(X_vect.at(i),X_vect.at(i+1))[1][0]));
            if (i == X_vect.size()-2 && k == -1)
            {
               tripletList.push_back(M(i+1, i+k+1 ,ComputeElementMassMatrix(X_vect.at(i),X_vect.at(i+1))[1][0])); 
            }
        }
    }
   }
    tripletList.shrink_to_fit();
    A.setFromTriplets(tripletList.begin(),tripletList.end());
    return A;
}

//Global Force matrix
Eigen::VectorXd GlobalForceVector(std::vector<double> X_vect)
{
    Eigen::VectorXd force(X_vect.size());
    for (int i = 0; i <= X_vect.size()-1; i++)
    {
        if( i == 0) 
        {
            force(i) = ComputeElementForceVector(X_vect.at(i),X_vect.at(i+1))[0];
        }
        else if (i == X_vect.size()-1)
        {   
            force(i) = ComputeElementForceVector(X_vect.at(i-1),X_vect.at(i))[1];
        }
        else
        {
            force(i) = ComputeElementForceVector(X_vect.at(i-1),X_vect.at(i))[1] + ComputeElementForceVector(X_vect.at(i),X_vect.at(i+1))[0];
        }
    }
   return force;
}

//Saving vetor to a file
void saveData(string fileName, Eigen::VectorXd  vector)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision ,Eigen::DontAlignCols,", ", "\n");
 
    ofstream file(fileName);
    if (file.is_open())
    {
        file << vector.format(CSVFormat);
        file.close();
    }
    
}

int main()
{
    double a = 0;
    double b = 4 * atan(1);
    double n = 20 ;
    int nel = n - 1;
    int n_int = n;
    vector<double> X_vect1;
    for (int i = 0; i < n;i++)
    {
       X_vect1.push_back(a + i*(b-a)/(n-1));
    }
    Eigen::VectorXd X_vect_Eigen(n_int);
    for (int i = 0 ; i < n ; i++)
    {
        X_vect_Eigen(i) = a + i*(b-a)/(n-1);
    }
 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(GlobalMassMatrix(X_vect1));
    solver.factorize(GlobalMassMatrix(X_vect1));
    Eigen::VectorXd x(X_vect1.size());
    x = solver.solve(GlobalForceVector(X_vect1));
    cout<<"Solution : \n \t"<<endl;
    cout<<x<<endl;
    saveData("Y_Vectors_FEM_HW2.csv",x);
    saveData("X_Vectors_FEM_HW2.csv",X_vect_Eigen);
    
}
