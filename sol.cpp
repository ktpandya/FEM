#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "Eigen/Eigen"

using namespace std;

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
        n2[current_line].erase(5,10);
        n2[current_line].erase(0,2);
        n3[current_line].erase(0,5);
       
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

void nodalCoords(int elm_num,int loc_node_num,const Mesh &M,vector<double>& nodal_coords)
{
    
    int glob_node_num = Local2GlobalMap(M,elm_num,loc_node_num);
    nodal_coords.push_back(M.coordinates[2 * glob_node_num]);
    nodal_coords.push_back(M.coordinates[2 * glob_node_num + 1 ]);

//    print_double(nodal_coords);
}

/*void ComputeKe(const std::vector<double>& nodal_coords, vector<vector<double>>& Ke)
{
    
    vector<double> X = {nodal_coords[0],nodal_coords[2],nodal_coords[4]};
    vector<double> Y = {nodal_coords[1],nodal_coords[3],nodal_coords[5]};
    double Area = 0.5 * ( X[0] * (Y[1] - Y[2]) + X[1] * (Y[2] - Y[0]) + X[2] * (Y[0] - Y[1]));
    vector <double> b = {(Y[1] - Y[2])*0.5/Area,(Y[2] - Y[0])*0.5*Area,(Y[0] - Y[1])*0.5*Area};
    vector <double> c = {(X[2] - X[1])*0.5/Area,(X[0] - X[2])*0.5/Area,(X[1] - X[0])*0.5/Area};
    vector<double> row1 = {b[0]*b[0] + c[0]*c[0], b[0]*b[1] + c[0]*c[1], b[0]*b[2] + c[0]*c[2] };
    vector<double> row2 = {b[1]*b[0] + c[1]*c[0], b[1]*b[1] + c[1]*c[1], b[1]*b[2] + c[1]*c[2] };
    vector<double> row3 = {b[2]*b[0] + c[2]*c[0], b[2]*b[1] + c[2]*c[1], b[2]*b[2] + c[2]*c[2] };
   
   Ke.push_back(row1);
   Ke.push_back(row2);
   Ke.push_back(row3);
   
   //printmatrix(Ke);
}*/

void ComputeKe(const std::vector<double>& nodal_coords,double* Ke)
{
    Ke = 
}


void ComputeFe(const std::vector<double>& nodal_coords)
{
    vector<double> Fe = {0,0,0};

}


int main()
{
    Mesh M;
    ReadMesh("coordinates.dat","connectivity.dat",M);
    //int elm_num = 128;
    //int loc_node_num = 3;
    //int globalindex = Local2GlobalMap(M,elm_num,loc_node_num);

    //Eigen::SparseMatrix<double> GlobalStiffnessMatrix(M.nNodes,M.nNodes);
    vector<double> nodal_coords;
    nodalCoords(126,M,nodal_coords);

    //cout<<nodal_coords.size();
    //ComputeKe(nodal_coords);

    vector<vector<double>> Ke;

    Eigen::SparseMatrix<double> A(M.nNodes,M.nNodes);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList{};
   tripletList.reserve(3*M.nNodes);
   for(int i = 0; i< M.nElements-1;i++)
   {
    for (int k = -1; k<2; k++)
    { 
        if (i==i+k )
        {
            nodalCoords(i,M,nodal_coords);
            ComputeKe(nodal_coords,Ke);
            tripletList.push_back(T(i,i+k,Ke[0][0]));
            tripletList.push_back(M(i+1,i+k+1,Ke[1][1]));
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