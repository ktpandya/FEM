#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>


struct Mesh
{
    int                      nNodes;
    std::vector<double> coordinates;
    int                   nElements;
    std::vector<int>   connectivity;
    std::vector<int>        bdNodes;

};


std::vector<std::string> ReadInputs()
{
    std::fstream Inputs;
    Inputs.open("connectivity.dat",std::ios::in);

    std::vector<std::string> inp;
    std::string line;
    int current_line = 0;
    while (std::getline(Inputs, line)) 
    {   
        
        inp.push_back(line);
        std::cout<<inp[current_line]<<std::endl;
        current_line = current_line+1;  
    }
    
    return inp;
}

void printvector(std::vector<std::vector<char>> X)
{
    for(int i = 0; i < X.size();i++)
    {
        for (int j = 0 ; j < X[i].size();j++)
        {
            std::cout<<X[i][j];
        }
    std::cout<<""<<std::endl;
    }

}

void print(std::vector<char> X)
{
    for (int i = 0; i < X.size();i++)
    {
        std::cout<<X[i]<<std::endl;
    }
}


void ReadMesh2(std::string coord_filename,std::string conn_filename, Mesh &M)
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

    f1.clear();
    
    std::vector<char> f2;
    std::ifstream file_conn;
    file_conn.open(conn_filename);
    if(!file_conn)
    {
    std::cerr <<"Problem opening file" << std::endl;
    }
    while(file_conn)
    {
    f1.push_back(file_conn.get());
    }
    for (int i = 31; i < f2.size(); i++)
    {
        f2.push_back(f1[i]);
    }


 
}


int main()
{
    Mesh M;
    //ReadMesh2("coordinates.dat","connectivity.dat",M);
    std::vector<std::string> inp;
    inp = ReadInputs();

    

}