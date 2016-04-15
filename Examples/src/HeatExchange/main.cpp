#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
  
  // Gauss Siedel is initialised with a linear variation
  // of T
  
  for(unsigned int m=0;m <= M;++m)
     theta[m]=(1.-m*h)*(To-Te)/Te;
  
  //Costruzione della matrice tramite 3 vettori
  std::vector<double> a(M), alpha(M),z(M);
  std::vector<double> b(M-1), beta(M-1), y(M);
  std::vector<double> c(M-1), gamma(M-1);
  for(int i=0; i<M-2; ++i){
     a[i]= 2+(hc*hc)*act;
     b[i]= -1;
     c[i]=-1;
  }
  //b[M-3]=-1;
  a[M-1]=1;
  //Costruzione dei vettori alpha, beta e gamma per le matrici L e U
  alpha[0]=a[0];
  for (int k=0; k<M-1; ++k) {
    beta[k]=b[k];
    gamma[k]=c[k]/alpha[k];
    alpha[k+1]=a[k+1]-beta[k]*gamma[k];
  }
  //Risoluzione del sistema
  z[0]=To-Te;
  y[0]=To-Te;
  for(int i=0; i<M-1; ++i){
     z[i+1]=0;
     y[i+1]=z[i+1]-gamma[i]*y[i];
  }
  theta[M]=y[M-1]/alpha[M-1];
  for(int j=M-1; j>0; j--){
     theta[j]=(y[j-1]-beta[j-1]*theta[j+1])/alpha[j-1];
  }
  theta[0]=To-Te;  
  theta[M]=y[M-1]/alpha[M-1];

  //trasformazione per riportarci a T e non a theta
  for (int k=0; k<M+2; ++k){
    theta[k]=theta[k]+Te;
  }


 // Analitic solution
    
    vector<double> thetaa(M+1);
     for(int m=0; m <= M;m++)
        thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));
 
      // writing results with format
      // x_i u_h(x_i) u(x_i) and lauch gnuplot 
 
      Gnuplot gp;

      std::vector<double> coor(M+1);
      std::vector<double> sol(M+1);
      std::vector<double> exact(M+1);
 
      cout<<"Result file: result.dat"<<endl;
      ofstream f("result.dat");
 
      for(int m = 0; m<= M; m++)
        {
 	 // \t writes a tab 
          f<<m*h*L<<"\t"<<theta[m]<<"\t"<<thetaa[m]<<endl;
 	 // An example of use of tie and tuples!
         
 	 std::tie(coor[m],sol[m],exact[m])=
 	   std::make_tuple(m*h*L,theta[m],thetaa[m]);
        }
 
     // Using temporary files (another nice use of tie)
      gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
        "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
        "w l title 'uex'"<<std::endl;
     f.close();
     return status;
}
