#include <iostream>
#include "../Eigen/Eigen/Sparse"
#include <chrono>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

double pi = 4.0*atan(1.0);

struct Mesh
{
    int nNodes;
    vector<double> coordinates;
    int nElements;
    vector<int> connectivity;
    vector<int> bdNodes;
};

void ReadMesh(std::string coord_filename, std::string conn_filename, Mesh& M)
{
    ifstream fin1(coord_filename);
	string data1;
    if(fin1.is_open())		// Checking if input file found and opened successfully
    {
		cout << "Opening coordinates file: Successful."<<endl;
        while(getline(fin1,data1,'\n'))
        {
            M.nNodes=stod(data1);
            cout<<M.nNodes<<endl;
            string rawInput;
            while(getline(fin1, rawInput,' '))
            {
                M.coordinates.push_back(stod(rawInput));
            }
        }
        fin1.close();
        cout << "Closing coordinates file: Successful."<<endl;
    }
    else
    {
        cout << "Opening coordinates file: Failed."<<endl;
		abort();
    }
    // cout<<sizeof(M.coordinates)<<endl;
    // for (int i = 0; i < M.nNodes; i++)
    // {
    //     cout<<M.coordinates[2*i]<<"\t"<<M.coordinates[2*i+1]<<"\n";
    // }


    ifstream fin2(conn_filename);
	string data2;
    if(fin2.is_open())		// Checking if input file found and opened successfully
    {
		cout << "Opening connectivity file: Successful."<<endl;
        while(getline(fin2,data2,'\n'))
        {
            M.nElements=stod(data2);
            cout<<M.nElements<<endl;
            string rawInput;
            while(getline(fin2, rawInput,' '))
            {
                M.connectivity.push_back(stod(rawInput));
            }
        }
        fin1.close();
        cout << "Closiing connectivity file: Successful."<<endl;
    }
    else
    {
        cout << "Opening connectivity file: Failed."<<endl;
		abort();
    }
    // cout<<sizeof(M.connectivity)<<endl;
    // for (int i = 0; i < M.nElements; i++)
    // {
    //     cout<<M.connectivity[3*i]<<"\t"<<M.connectivity[3*i+1]<<"\t"<<M.connectivity[3*i+2]<<"\n";
    // }

}

int Local2GlobalMap(const Mesh& M, int elm_num, int loc_node_num)
{
    return 0;
}

void ComputeKe(const std::vector<double>& nodal_coords, double* Ke)
{

}

void ComputeFe(const std::vector<double>& nodal_coords, double* Fe)
{

}

void SetDirichletBCs(const Mesh& M, double& K, double& F)
{

}

void PrintSolution(std::string filename, const Mesh& M, const double* U)
{
    
}

int main()
{
    struct Mesh mesh;
    ReadMesh("coordinates.dat","connectivity.dat", mesh);
    cout<<"\nCoordinates:\nx\ty\n";
    for (size_t i = 0; i < mesh.nNodes; i++)
    {
        cout<<mesh.coordinates[2*i]<<"\t"<<mesh.coordinates[2*i+1]<<"\n";
    }
    cout<<"\nConnectivity:\nV1\tV2\tV3\n";
    for (int i = 0; i < mesh.nElements; i++)
    {
        cout<<mesh.connectivity[3*i]<<"\t"<<mesh.connectivity[3*i+1]<<"\t"<<mesh.connectivity[3*i+2]<<"\n";
    }
    
    cout<<"Okay\n";
}