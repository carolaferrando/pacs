#include <cmath>
#include <algorithm>
#include "rk45.hpp"
#include <vector>
using namespace std;


class rk45_butcher{
	vector<vector<double>> A;
	vector<vector<double>> B;
	vector<double> C;
public:
	rk45_butcher(){
	    vector <double> a={1./4.};
	    A.push_back(a);
	    C.push_back(1./4);
	    a={3./32, 9./32};
	    A.push_back(a);
	    C.push_back(3./32+9./32);
	    a={1932./2197, -7200./2197, 7296./2197};
	    A.push_back(a);
	    C.push_back(1932./2197 -7200./2197+7296./2197);
	    a={439./216.,- 8., 3680./513.,-845./4104.};
	    A.push_back(a);
	    C.push_back(439./216.- 8.+ 3680./513.-845./4104.);
	    a={-8./27.,2.,-3544./2565., 1859./4104.,-11./40.};
	    A.push_back(a);
	    C.push_back(-8./27.+2.-3544./2565.+ 1859./4104.-11./40.);
		
	    vector <double> b={25./216.,0.,1408./2565.,2197./4104.,-1./5.};
	    B.push_back(b);
	    b={16./135.,0.,6656./12825.,28561./56430.,-9./50.,2./55.};
	    B.push_back(b);
	};
	double geta(int i,int j){
		return A[i-2][j-1];
	};
	double getb(int i,int j){
		return B[i-4][j-1];
	};
	double getc(int i){
		return C[i-2];
	};
		
};

class rk23_butcher{
	vector<vector<double>> A;
	vector<vector<double>> B;
	vector<double> C;
public:
	rk23_butcher(){
	    vector <double> a={1./2.};
	    A.push_back(a);
	    C.push_back(1./2);
	    a={0, 3./4};
	    A.push_back(a);
	    C.push_back(3./4);
	    a={2./9, 1./3, 4./9};
	    A.push_back(a);
	    C.push_back(1);
		
	    vector <double> b={2./9.,1./3,4./9.,0};
	    B.push_back(b);
	    b={7./24.,1./4,1./3.,1./8.};
	    B.push_back(b);
	};
	double geta(int i,int j){
		return A[i-2][j-1];
	};
	double getb(int i,int j){
		return B[i-4][j-1];
	};
	double getc(int i){
		return C[i-2];
	};
		
};


 double rk45_step
  (
   std::function<double (double const &, double const &)> const & f,
   double const & y0, 
   double const & t0, 
   double const & h, 
   double & error,
   rk45_butcher butch
   )
  {
   
    
    double F1 = h * f(t0, y0);
    double F2 = h * f(t0 + butch.getc(2) * h, y0 + butch.geta(2,1) * F1);
    double F3 = h * f(t0 + butch.getc(3) * h, y0 + butch.geta(3,1) * F1 + butch.geta(3,2) * F2);
    double F4 = h * f(t0 + butch.getc(4) * h, y0 + butch.geta(4,1) * F1 + butch.geta(4,2) * F2 + butch.geta(4,3) * F3);
    double F5 = h * f(t0 + butch.getc(5)* h, y0 + butch.geta(5,1) * F1 + butch.geta(5,2) * F2 + butch.geta(5,3) * F3 + butch.geta(5,4) * F4 );
    double F6 = h * f(t0 + butch.getc(6) * h, y0 + butch.geta(6,1) * F1 + butch.geta(6,2) * F2 + butch.geta(6,3) * F3 + butch.geta(6,4) * F4 + butch.geta(6,5) * F5);

    double y4 =   y0 + butch.getb(4,1) * F1 + butch.getb(4,3) * F3 + butch.getb(4,4) * F4 + butch.getb(4,5) * F5;
    double y5 =   y0 + butch.getb(5,1) * F1 + butch.getb(5,3) * F3 + butch.getb(5,4) * F4 + butch.getb(5,5) * F5 + butch.getb(5,6) * F6;
    error = std::abs(y5 - y4);
    return y5;
  }
  
template<typename rkType>

class generic_rk: public rkType{ // vedere se funziona private
	
public:
	std::vector<std::pair<double,double>> solve(
		std::function<double (double const &, double const &)> const & dy,
		double const & t0,
		double const & T,
		double const & y0,
		double const & h_initial, 
		double const & h_max, 
		double const & final_error,
		int & status,
		std::size_t const & maxSteps){
	
			status=0;
			const std::size_t maxReduction=maxSteps;
			// parameters for decreasing/increasing time step
			double const c1=1.0;
			// I need to have a sufficient decrease of the local error
			// to allow time step coarsening
			double const c2=1./64.;

			double length=T-t0;
			//! Make sure that h allows to reach T
			std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
			double h=length/initialNSteps;
			// To avoid underflow we need in any case to limit the time step to a positive number
			// Here I allow h to become 128 time smaller than that giving the maximal number of steps
			double h_min = length/(128*maxSteps);
			// SOme counters
			std::size_t stepsCounter(0);
			// Initial data
			double time(t0);
			double y(y0);
			double errorPerTimeStep=final_error/initialNSteps;
			if (initialNSteps>=maxSteps) throw std::runtime_error("RK45: initial time step h too small!");
			std::vector<std::pair<double,double>> solution;
			solution.emplace_back(std::make_pair(t0,y0));
			double localError;
			double newy;
			rk45_butcher butch;
			while (time<T && stepsCounter <maxSteps)
		 {
		//Do a step
		//adjust h if needed for the last step
		if (time + h > T) h = T-time;
		newy = rk45_step(dy,y,time,h,localError,butch);
		while (h> h_min && localError > c1*errorPerTimeStep)
			{
				// half time step
				h /=2;
				errorPerTimeStep /=2;
				newy = rk45_step(dy,y,time,h,localError,butch);
			}
		if (localError>errorPerTimeStep)status=1;
		//! advance
		y = newy;
		time +=h;
		++stepsCounter;
		solution.emplace_back(std::make_pair(time,y));
		//! check if we reached end
		if(localError<c2*errorPerTimeStep && h<h_max)
			{
				// Double step
				h *=2;
				errorPerTimeStep *=2;
			}
		 }
			//handle exceptions
			if(stepsCounter>=maxSteps && time < T)
		 {
		status=2;
		throw std::runtime_error("RK45: Max number of time steps exceeded");
		 }
			return solution;
		}
};
