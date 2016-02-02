#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cstdlib>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const f0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void step1(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int Nx, const cmplx alpha);
void step2(cmplx* const f,const double dt, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const cmplx alpha = cmplx(0.0, dt / pow(dx,2));
	
	const double eta = 0.2;

	stringstream strm;

	cmplx* f0 = new cmplx[Nx];
	cmplx* f1 = new cmplx[Nx];
	cmplx* h;

	init(f0, eta, dx, dt,Nx);

	writeToFile(f0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {
	
  	        step2(f0,dt*0.5,Nx);
		
		for (int j = 1; j <= Nk-2; j++) {
		  

		  step1(f1,f0,dt,dx,Nx,alpha);
 		  step2(f1,dt,Nx);
		  
		  h = f0;
		  f0 = f1;
		  f1 = h;
		  
		  
		  
		}

                step1(f1,f0,dt,dx,Nx,alpha);
	        step2(f1,dt*0.5,Nx);

      
		h = f0;
		f0 = f1;
		f1 = h;
	    
	    
	    
		strm.str("");
		strm << "psi_" << i;
		writeToFile(f0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const f0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		f0[i] = 2*f/cosh(eta * x);
	}
}

//-----------------------------------------------------------------------------------------------
void step1(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int Nx, const cmplx alpha)
{

  cmplx* d=new cmplx[Nx];
  cmplx* u=new cmplx[Nx];
  cmplx* l=new cmplx[Nx];

  for(int i=0;i<Nx;i++) d[i] = cmplx(1.0,0) - 2.0*alpha;
  for(int i=0;i<Nx;i++) u[i] = alpha;
  for(int i=0;i<Nx;i++) l[i] = alpha;

  
  for(int i=1;i<Nx;i++){
    d[i]  = d[i] - l[i]/d[i-1] * u[i-1];
    f0[i] = f0[i] - l[i]/d[i-1] * f0[i-1];
  }
  
  
  
  f1[Nx-1]=f0[Nx-1]/d[Nx-1];
  for(int i=Nx-2;i>=0;i--){
   f1[i]=(f0[i] - u[i]*f1[i+1])/d[i];
  }
    

  delete[] d;
  delete[] u;
  delete[] l;
}


void step2(cmplx* const f,const double dt, const int Nx){
  
  
  for(int i=0;i<Nx;i++) f[i] = f[i]* cmplx (cos(pow(abs(f[i]),2)*dt), -sin(pow(abs(f[i]),2))*dt); 
  
  
  
  
  
  
}
