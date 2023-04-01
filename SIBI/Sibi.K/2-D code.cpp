//Procedure to run the program
// Stores all the files given in the same directory
// First mesh the domain using the matlab file named "firstquestion" by providing Hmax as input
//The run the c++ program named "2-D code"
// The for the visualization of the FEM solution run the matlab file named "eightquestion"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include"Eigen/Eigen"
using namespace std;
const double pi=M_PI;

double gradphik[2][2];
double gradphikinv[2][2];
double N[3][2]={1,0,0,1,-1,-1};
double det;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenMatType;
typedef Eigen::VectorXd EigenVecType;
struct Mesh // Creating Mesh structure to varies data of the Mesh under the object M.
    {
        int nNodes;
        std::vector<double> coordinates;
        int nElements;
        std::vector<double>connectivity;
        int nBoundaryNodes;
        std::vector<double>bdNodes;


    };
    struct Mesh M;

//Reading the mesh information to various struct variables

void ReadMesh(string a,string b,string c,Mesh& M)
{



// To read the coordinates of nodes into M1.coordinates
    std::ifstream ifile(a, std::ios::in);



    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

    double num = 0.0;

    while (ifile >> num)
    {
        M.coordinates.push_back(num);
    }

cout<<endl<<"Displaying the coordinates of nodes in the global node order"<<endl;

    for (int i = 0; i < M.coordinates.size(); ++i)
    {
        std::cout << M.coordinates[i] << std::endl;
    }
    M.nNodes=(M.coordinates.size())/2;
    cout<<endl<<M.nNodes<<endl<<endl<<endl;

//To read the connectivity of every element into M1.connectivity

    std::ifstream jfile(b, std::ios::in);



    if (!jfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

     num = 0.0;

    while (jfile >> num) {
        M.connectivity.push_back(num);
    }

cout<<endl<<"Displaying the connectivity"<<endl;

    for (int i = 0; i < M.connectivity.size(); ++i)
    {
        std::cout << M.connectivity[i] << std::endl;
    }
    M.nElements=(M.connectivity.size())/3;
    cout<<endl<<M.nElements<<endl<<endl<<endl;


//To read the nodes in boundaries

std::ifstream kfile(c, std::ios::in);



    if (!kfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);
    }

     num = 0.0;

    while (kfile >> num) {
        M.bdNodes.push_back(num);
    }

cout<<endl<<"Displaying the nodes in the boundaries"<<endl;

    for (int i = 0; i < M.bdNodes.size(); ++i)
    {
        std::cout << M.bdNodes[i] << std::endl;
    }
    M.nBoundaryNodes=M.bdNodes.size();
    cout<<endl<<M.nBoundaryNodes<<endl<<endl;



}

//Function which names the elements with local nodes and creates a local to global map

double Local2GlobalMap(const Mesh& M,int elm_num,int loc_node_num )
{

int local[3*M.nElements];
double map1[3][M.nElements];

//Numbers the respective global nodes as local nodes using 1,2,3

for(int i=1;i<=3*M.nElements;i+=3)
{
local[i-1]=1;
local[i]=2;
local[i+1]=3;
}
/*for(int i=0;i<3*M.nElements;i++)
{
    cout<<local[i]<<" ";
}*/

//Creating the Local2Global in map1 matrix

for(int i=1;i<=M.nElements;i++)
{
    map1[0][i-1]=M.connectivity[(3*i)-3];
    map1[1][i-1]=M.connectivity[(3*i)-2];
    map1[2][i-1]=M.connectivity[(3*i)-1];
}

/*
for(int i=0;i<3;i++)
{

    for(int j=0;j<M.nElements;j++)
    {
        cout<<map1[i][j]<<" ";
    }
    cout<<endl;
}*/

double c=map1[loc_node_num-1][elm_num-1];

return(c);
}

//Helps to calculate the elements of the element stiffness matrix

double calc(int a,int b)
{
    double l;
    double l1[2];
    double l2[2];
    for(int i=0;i<2;i++)
    {
        l1[i]=((gradphikinv[i][0])*(N[a][0]))+((gradphikinv[i][1])*(N[a][1]));
    }
     for(int i=0;i<2;i++)
    {
        l2[i]=((gradphikinv[i][0])*(N[b][0]))+((gradphikinv[i][1])*(N[b][1]));
    }
    l=l1[0]*l2[0]+l1[1]*l2[1];
    return(l);
}

//Helps to calculate the element stiffness matrix

void ComputeKe(const std::vector<double>&nodal_coords,double*Ke)
{

   for(int i=0;i<2;i++)
   {
       for(int j=0;j<2;j++)
       {
           gradphik[i][j]=nodal_coords[(2*i)+j]-nodal_coords[j+4];
       }
   }
   det=gradphik[0][0]*gradphik[1][1]-gradphik[0][1]*gradphik[1][0];

cout<<endl<<det<<endl;
if(det==0)
{
    cout<<endl<<"Determinant of the jacobian was zero"<<endl;
    exit(0);
}



   gradphikinv[0][0]=(gradphik[1][1]/det);
   gradphikinv[0][1]=-(gradphik[0][1]/det);
   gradphikinv[1][0]=-(gradphik[1][0]/det);
   gradphikinv[1][1]=(gradphik[0][0]/det);


cout<<endl<<"Displaying gradphi inverse"<<endl;

for(int i=0;i<2;i++)
   {
       for(int j=0;j<2;j++)
       {
           cout<<gradphikinv[i][j]<<" ";
       }
       cout<<endl;
   }

cout<<endl<<endl;

for(int i=0;i<3;i++)
   {
       for(int j=0;j<3;j++)
       {
          double *x;
          x=Ke+(3*i)+j;
          *x=0.5*det*calc(i,j);
       }
   }

}

//Helps to calculate the element force vector

void ComputeFe(const std::vector<double>&nodal_coords,double*Fe)
{
double x1=nodal_coords[0]; double y1=nodal_coords[1];
double x2=nodal_coords[2]; double y2=nodal_coords[3];
double x3=nodal_coords[4];double y3=nodal_coords[5];
double*x_;
x_=Fe;
*x_=(- pow(x1,2)/20 - (x1*x2)/30 - (x1*x3)/30 + x1/12 - pow(x2,2)/60 - (x2*x3)/60 + x2/24 - pow(x3,2)/60 + x3/24 - pow(y1,2)/20 - (y1*y2)/30 - (y1*y3)/30 + y1/12 - pow(y2,2)/60 - (y2*y3)/60 + y2/24 - pow(y3,2)/60 + y3/24)*det;
x_=Fe+1;
*x_=(-pow(x1,2)/60 - (x1*x2)/30 - (x1*x3)/60 + x1/24 - pow(x2,2)/20 - (x2*x3)/30 + x2/12 - pow(x3,2)/60 + x3/24 - pow(y1,2)/60 - (y1*y2)/30 - (y1*y3)/60 + y1/24 - pow(y2,2)/20 - (y2*y3)/30 + y2/12 - pow(y3,2)/60 + y3/24)*det;
x_=Fe+2;
*x_=(-pow(x1,2)/60 - (x1*x2)/60 - (x1*x3)/30 + x1/24 - pow(x2,2)/60 - (x2*x3)/30 + x2/24 - pow(x3,2)/20 + x3/12 - pow(y1,2)/60 - (y1*y2)/60 - (y1*y3)/30 + y1/24 - pow(y2,2)/60 - (y2*y3)/30 + y2/24 - pow(y3,2)/20 + y3/12)*det;
}

//Sets the Dirchlet Boundrary conditions to system

void SetDirchletBCs(const Mesh& M, EigenMatType& K, EigenVecType& F)
{
for(int i=0;i<M.nBoundaryNodes;i++)
{
    int h;
    h=M.bdNodes[i]-1;
    for(int j=0;j<M.nNodes;j++)
    {
        if(j==h)
        {
           K.coeffRef(h,j)=1;
        }
        else
            K.coeffRef(h,j)=0;
    }
    F.coeffRef(h)=0;
}

}

//Stores the dof vector in a .dat file named sol

void PrintSolution(std::string filename,const Mesh&M,  double*U)
{

double m;
std::fstream pfile;
pfile.open(filename,std::ios::out);
for(int i=0;i<M.nNodes;i++)
{
double *x__;
x__=U+i;
m=*x__;

pfile<<m<<"\n";

}
pfile.close();


}

//Helps to compete the error(L2 norm error)

void Competeerror(double *dof)
{
    double u[M.nNodes],u_,v_,rmserror,sum=0,k_;
    for (int i=0;i<M.nNodes;i++)
    {
        u_=M.coordinates[2*i];
        v_=M.coordinates[(2*i)+1];
        u[i]=0.5*u_*v_*(1-u_)*(1-v_);
    }
    for(int j=0;j<M.nNodes;j++)
    {    k_=*(dof+j);

        sum=sum+pow((u[j]-k_),2);


        }
    rmserror=sqrt(sum/M.nNodes);

    cout<<endl<<"The rmserror is "<<rmserror<<endl;

}

int main()
{
string file1="coordinates.txt",file2="connectivity.txt",file3="boundarynodes.txt",file4="sol.dat";

ReadMesh(file1,file2,file3,M);

typedef Eigen::Triplet<double>T;
std::vector<T>tripletList{};
tripletList.reserve(9*M.nElements);
EigenVecType ForceVector(M.nNodes);
ForceVector.setZero();

EigenVecType dofvector(M.nNodes);

std::vector<double>nordal;
double ke[3][3];
double fe[3];
double e;
double f;

    for(int k=1;k<=M.nElements;k++)
    {

        for(int i=1;i<=3;i++)
        {

        nordal.push_back(M.coordinates[2*Local2GlobalMap(M,k,i)-2]);
        nordal.push_back(M.coordinates[2*Local2GlobalMap(M,k,i)-1]);

        }
        cout<<endl<<"Displaying the nordal coordinates of the element number "<<k<<endl;

            for(int i=0;i<6;i++)
        {
            cout<<nordal[i]<<" ";

        }
        cout<<endl<<endl;
        ComputeKe(nordal,&ke[0][0]);

        cout<<"Displaying local stiffness matrix of element number "<<k<<endl;

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                cout<<ke[i][j]<<" ";
            }
        cout<<endl;
        }


        for(int i=1;i<=3;i++)
        {
            for(int j=1;j<=3;j++)

            {
                int c,d;
                c=Local2GlobalMap(M,k,i)-1;
                d=Local2GlobalMap(M,k,j)-1;
                tripletList.push_back(T(c,d,ke[i-1][j-1]));
            }
        }


        ComputeFe(nordal,&fe[0]);

        cout<<endl<<"Displaying the forcevector of the element number "<<k<<endl;

        for(int i=0;i<3;i++)
        {
            cout<<fe[i]<<endl;
        }



        for(int i=1;i<=3;i++)
        {
        e=Local2GlobalMap(M,k,i)-1;
        f=fe[i-1];
        ForceVector.coeffRef(e)+=f;
        }              

        nordal.clear();
    }

tripletList.shrink_to_fit();
EigenMatType Stiffness(M.nNodes,M.nNodes);
Stiffness.setFromTriplets(tripletList.begin(),tripletList.end());

cout<<endl<<"Global Stiffness Matrix"<<endl;

for(int i=0;i<M.nNodes;i++)
{

    for(int j=0;j<M.nNodes;j++)
    {
        cout<<Stiffness.coeffRef(i,j)<<" ";
    }
cout<<endl;
}
cout<<endl;
cout<<"Global forcevector"<<endl;

for(int i=0;i<M.nNodes;i++)
{
    cout<<ForceVector.coeffRef(i)<<endl;
}


SetDirchletBCs(M,Stiffness,ForceVector);


cout<<endl<<"Global Stiffness Matrix after BC's"<<endl;

for(int i=0;i<M.nNodes;i++)
{

    for(int j=0;j<M.nNodes;j++)
    {
        cout<<Stiffness.coeffRef(i,j)<<" ";
    }
cout<<endl;
}
cout<<endl;
cout<<"Global forcevector after BC's"<<endl;

for(int i=0;i<M.nNodes;i++)
{
    cout<<ForceVector.coeffRef(i)<<endl;
}


//Soling KU=F

Eigen::SparseLU<EigenMatType> solver;
solver.analyzePattern(Stiffness);
solver.factorize(Stiffness);
dofvector = solver.solve(ForceVector);

cout<<endl<<"dof vector"<<endl;
for(int i=0;i<M.nNodes;i++)
{
    cout<<dofvector.coeffRef(i)<<endl;
}
cout<<endl;
/*
for(int i=0;i<M.nNodes;i++)
{

    double *n;
    n=&dofvector(i);
    cout<<n<<endl;

}
cout<<endl;
// /

PrintSolution(file4,M,&dofvector(0));


Competeerror(&dofvector(0));

return 0;

}


