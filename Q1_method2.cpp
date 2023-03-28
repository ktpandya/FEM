#include <iostream>
#include <chrono>
#include "Eigen/Dense"
#include "Eigen/Eigen"
using Eigen::MatrixXd;

int main()
{
    
    // start clock
    std::chrono::steady_clock::time_point t_begin =std::chrono::steady_clock::now();
    // Your-code-goes-here
    const int N = 5000;
    Eigen::SparseMatrix<double> A(N,N);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList{};
    tripletList.reserve(3*N);
    for(int i=0; i<N; ++i)
        {for(int k=-1; k<=1; ++k)
            {if(i+k>=0 && i+k<N)
                {tripletList.push_back(T(i,i+k,1.0));
                }}}
    tripletList.shrink_to_fit();
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // stop clock
    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    // Track the run time
    double run_time = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin).count();
    std::cout<<"\n Run time is "<<run_time<<"\n";
}