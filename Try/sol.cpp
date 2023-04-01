#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "Eigen/Eigen"

using namespace std;
typedef Eigen::SparseMatrix<double> EigenMatType;
typedef Eigen::VectorXd EigenVecType;

struct Mesh
{
    int                      nNodes;
    std::vector<double> coordinates;
    int                   nElements;
    std::vector<int>   connectivity;
    std::vector<int>        bdNodes;

};

void printmatrix(std::vector<vector<double>> M)
{
    for(int i = 0; i < M.size(); i ++)
    {
        for(int j = 0; j < M.size();j++)
        {
            cout<<M[i][j]<<"\t";
        }
        cout<<"\n";

    }
}

void print_int(std::vector<int> X)
{
    for (int i = 0; i < X.size();i++)
    {
        std::cout<<X[i]<<std::endl;
    }
}

void print_double(std::vector<double> X)
{
    for (int i = 0; i < X.size();i++)
    {
        std::cout<<X[i]<<std::endl;
    }
}

void ReadMesh(std::string coord_filename,std::string conn_filename, Mesh &M)
{
    std::ifstream file_coord;
    std::vector<char> f1;
    std::vector<char> f;

    file_coord.open(coord_filename);
    if(!file_coord)
    {
    std::cerr <<"Problem opening file" << std::endl;
    }
    while(file_coord)
    {
    f1.push_back(file_coord.get());
    }
    for (int i = 50; i < f1.size(); i++)
    {
        f.push_back(f1[i]);
    }

    std::vector<char> x;
    std::vector<double> X;
    std::vector<char> y;
    std::vector<double> Y;
    std::string s_x;
    std::string s_y;
    double d_x;
    double d_y;

    for(int i = 4 ; i < f.size();i++)
    {
        if(f[i] != ' ')
        {   
            x.clear();
            y.clear();

            std::ostringstream out_x;
            std::ostringstream out_y;

            for(int j = i; j < i + 18; j++)
            {
                x.push_back(f[j]);
                y.push_back(f[j+24]);
            }   
            
            for (char c: x) 
            {
                out_x << c;
            }

            for (char c: y) 
            {
                out_y << c;
            }
            
            s_x = out_x.str();
            s_y = out_y.str();
            
            d_x = std::stod(s_x);
            d_y = std::stod(s_y);
            
            
            X.push_back(d_x); //Ans
            Y.push_back(d_y);
            i = i + 70; 
                    
            out_x.clear();
            out_y.clear();
        }
    }
   for(int iter= 0; iter<X.size();iter++)
   {
    M.coordinates.push_back(X[iter]);
    M.coordinates.push_back(Y[iter]); 
   }
    
    

    f1.clear();
    

    std::fstream Inputs;
    Inputs.open("connectivity.dat",std::ios::in);

    
    std::vector<std::vector<int>> connectivityMatrix;

    std::vector<std::string> inp;
    std::vector<std::string> inp2;
    std::string line;
    int current_line = 0;
    while (std::getline(Inputs, line)) 
    {   
     
        inp.push_back(line);
        inp[current_line].erase(0,11);
       
        current_line = current_line+1;  
    }
    for(current_line = 2; current_line < inp.size();current_line++)
    {
        inp2.push_back(inp[current_line]);        
    }

    std::vector<string> n1(inp2.size());
    std::vector<string> n2(inp2.size());
    std::vector<string> n3(inp2.size());
    
    std::vector<int> connectivity(3);

    for(current_line = 0 ; current_line < inp2.size(); current_line++)
    {   
        connectivity.clear();
        n1[current_line] = inp2[current_line];
        n2[current_line] = inp2[current_line];
        n3[current_line] = inp2[current_line];
        n1[current_line].erase(2,10);
   
       

        if(inp2.size()<10)
        {
            n2[current_line].erase(4,10);
            n2[current_line].erase(0,2);
            n3[current_line].erase(0,4);
        }
        else
        {
            n2[current_line].erase(5,10);
            n2[current_line].erase(0,2);
            n3[current_line].erase(0,5);
        }

        size_t pos;
        std::string x = " ", y = "";
        while ((pos = n1[current_line].find(x)) != std::string::npos) 
        {
            n1[current_line].replace(pos, 1, y);
        }
        while ((pos = n2[current_line].find(x)) != std::string::npos) 
        {
            n2[current_line].replace(pos, 1, y);
        }
        while ((pos = n3[current_line].find(x)) != std::string::npos) 
        {
            n3[current_line].replace(pos, 1, y);
        }
        //connectivity.push_back(stoi(n1[current_line]));
        //connectivity.push_back(stoi(n2[current_line]));
        //connectivity.push_back(stoi(n3[current_line]));
        //connectivityMatrix.push_back(connectivity);
        //connectivity.clear();
        M.connectivity.push_back(stoi(n1[current_line]));
        M.connectivity.push_back(stoi(n2[current_line]));
        M.connectivity.push_back(stoi(n3[current_line]));
    }
  
   M.nNodes = M.coordinates.size()/2;
   M.nElements = M.connectivity.size()/3;


   ////###############################################################



   f1.clear();
    

    std::fstream Inputs3;
    Inputs3.open("boundarynodes.dat",std::ios::in);

    
    std::vector<int> bdnodes;

    std::vector<std::string> inp_bd;
    std::vector<std::string> inp2_bd;
   
    current_line = 0;
    while (std::getline(Inputs3, line)) 
    {   
     
        inp_bd.push_back(line);
        inp_bd[current_line].erase(0,22);
        
        
        current_line = current_line+1;  
    }
    for(current_line = 3; current_line < inp_bd.size();current_line++)
    {
        inp2_bd.push_back(inp_bd[current_line]);  
        inp2_bd[current_line-3].erase(3,10);      
    }

    std::vector<string> n(inp2_bd.size());
    
    
   

    for(current_line = 0 ; current_line < inp2_bd.size(); current_line++)
    {   
        
        n[current_line] = inp2_bd[current_line];
        
       
        size_t pos;
        std::string x = " ", y = "";
        while ((pos = n[current_line].find(x)) != std::string::npos) 
        {
            n[current_line].replace(pos, 1, y);
        }
       
        //connectivity.push_back(stoi(n1[current_line]));
        //connectivity.push_back(stoi(n2[current_line]));
        //connectivity.push_back(stoi(n3[current_line]));
        //connectivityMatrix.push_back(connectivity);
        //connectivity.clear();
        M.bdNodes.push_back(stoi(n[current_line]));
        
    }
  
   ////###############################################################
   
}

int Local2GlobalMap(const Mesh& M, int elm_num, int loc_node_num)
{
    int globalindex =0;

    if(elm_num > M.nElements)
    {
        cout<<"Element is not present!!";
        
    }
    else
    {
        globalindex = M.connectivity[3 * (elm_num-1) + loc_node_num-1];
    }
    return globalindex;
}

void nodalCoords(int elm_num,const Mesh &M,vector<double>& nodal_coords)
{
    
    for(int j = 1; j < 4 ; j++)
    {
        nodal_coords.push_back(M.coordinates[ 2 * ( Local2GlobalMap ( M , elm_num , j ) - 1 )  ]);
        nodal_coords.push_back(M.coordinates[ 2 * ( Local2GlobalMap ( M , elm_num , j ) - 1 ) +1 ]);
    }

//    print_double(nodal_coords);
}


void ComputeKe(const std::vector<double>& nodal_coords, double *Ke)
{
    
    vector<double> X = {nodal_coords[0],nodal_coords[2],nodal_coords[4]};
    vector<double> Y = {nodal_coords[1],nodal_coords[3],nodal_coords[5]};
    double Area = 0.5 * ( X[0] * (Y[1] - Y[2]) + X[1] * (Y[2] - Y[0]) + X[2] * (Y[0] - Y[1]));
    vector <double> b = {(Y[1] - Y[2])*0.5/Area,(Y[2] - Y[0])*0.5*Area,(Y[0] - Y[1])*0.5*Area};
    vector <double> c = {(X[2] - X[1])*0.5/Area,(X[0] - X[2])*0.5/Area,(X[1] - X[0])*0.5/Area};
    
    
    for(int i=0;i<3;i++)
   {
       for(int j=0;j<3;j++)
       {
          double *pointer;
          pointer=Ke+(3*i)+j;
          *pointer= b[j] * b[i] + c[j] * c[i];
       }
   }
   
   
  
}

double func(double x, double y)
{
    return -32*x+32*x*x-32*y+32*y*y;
}

void ComputeFe(const std::vector<double>& nodal_coords, double *Fe)
{
    vector<double> X = {nodal_coords[0],nodal_coords[2],nodal_coords[4]};
    vector<double> Y = {nodal_coords[1],nodal_coords[3],nodal_coords[5]};
    double Area = 0.5 * ( X[0] * (Y[1] - Y[2]) + X[1] * (Y[2] - Y[0]) + X[2] * (Y[0] - Y[1]));

    double*pointer;
   
    for(int i= 0; i < 3; i ++)
    {
        pointer = Fe+i;
        *pointer = Area * func(X[i],Y[i])/3;   
    }
}


void SetDirchletBCs(const Mesh& M, EigenMatType& K, EigenVecType& F)
{
    for(int i=0;i<M.bdNodes.size();i++)
    {
        
        for(int j=0;j<M.nNodes;j++)
        {
            if(j==M.bdNodes[i]-1)  //-1 to take care of indexing
            {
                K.coeffRef(M.bdNodes[i]-1,j)=1;
            }
            else
                K.coeffRef(M.bdNodes[i]-1,j)=0;
        }
        F.coeffRef(M.bdNodes[i]-1)=0;
    }

}

void PrintSolution(std::string filename,const Mesh&M,  double*U)
{

    
    std::fstream file;
    file.open(filename,std::ios::out);
    for(int i=0;i<M.nNodes;i++)
    {
        double *pointer;
        pointer=U+i;
        file<<*pointer<<"\n";
    }
    file.close();


}

double Error( const Mesh &M)
{
    
    std::fstream file;
    file.open("output",std::ios::in);

    std::vector<double> inp;
    
    std::string line;
    int current_line = 0;
    while (std::getline(file, line)) 
    {   
     
        inp.push_back(stod(line));
        
       
        current_line = current_line+1;  
    }

    cout<<inp.size()<<endl;
    vector<double> X;
    vector<double> Y;
    double L2 =0;

    for(int i = 0; i < M.nNodes;i++)
    {
        X.push_back(M.coordinates[2*i]);
        Y.push_back(M.coordinates[2*i + 1]);
    }
    vector<double> exact;
    for(int i = 0; i < M.nNodes; i++)
    {     
        //cout<<"Exact = "<<(16 * X[i] * (1-X[i])*Y[i]*(1-Y[i]))<<endl;
       L2 += pow((inp[i] - 16 * X[i] * (1-X[i])*Y[i]*(1-Y[i])),2);

    }
    L2 = sqrt(L2/(M.nNodes+1));
    X.clear();
    Y.clear();
    cout<<"L2 = " <<L2;
    return L2;


}


int main()
{
    Mesh M;
    ReadMesh("coordinates.dat","connectivity.dat",M);
  
    
    vector<double> nodal_coords;
    nodal_coords.clear();
    
    double Ke[3][3];
    double Fe[3];

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList{};
    tripletList.reserve(9*M.nElements);

    Eigen::SparseMatrix<double> A(M.nNodes,M.nNodes);
    Eigen::VectorXd b(M.nNodes);
    
    for(int elm_num = 1; elm_num<=M.nElements;elm_num++)
    {
        nodal_coords.clear();
        
        nodalCoords(elm_num,M,nodal_coords);
        ComputeKe(nodal_coords,&Ke[0][0]);
        ComputeFe(nodal_coords,&Fe[0]);
        
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                tripletList.push_back(T(Local2GlobalMap(M,elm_num,i+1)-1,Local2GlobalMap(M,elm_num,j+1)-1,Ke[i][j]));
            }
        }

        for(int i=1;i<=3;i++)
        {
        
            double f=Fe[i-1];
            b.coeffRef(Local2GlobalMap(M,elm_num,i)-1)+=f;
        }           


    }


    tripletList.shrink_to_fit();
    A.setFromTriplets(tripletList.begin(),tripletList.end());
   
    
cout<<endl<<"Global Stiffness Matrix"<<endl;

for(int i=0;i<M.nNodes;i++)
{

    for(int j=0;j<M.nNodes;j++)
    {
        cout<<A.coeffRef(i,j)<<" ";
    }
cout<<endl;
}

cout<<endl;

    SetDirchletBCs(M,A,b);


    EigenVecType x(M.nNodes);

    Eigen::SparseLU<EigenMatType> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    PrintSolution("output",M,&x(0));

    Error(M);


}