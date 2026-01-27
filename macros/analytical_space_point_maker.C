#include "TMath.h"

const int n = 3;
double a[n];
double b[n];
double c[n];

double func(const double* xx){
  double t = xx[0];
  double alpha = xx[1];
  double d2 = 0;
  for (int i=0;i<n;i++){
    double d = a[i]*t*cos(alpha) + b[i]*t*sin(alpha) + c[i];
    d2+=d*d;
  }
  return d2;
}

void analytical_space_point_maker(){
  double thetaTrue = 30*TMath::DegToRad();
  double alphaTrue = 14*TMath::DegToRad();
  double tTrue = tan(thetaTrue);
  double kTrue = tan(alphaTrue);
  double r[] = { 100,  100,  100}; // station minimal radius
  double z[] = {2100, 2110, 2120}; // station pozition
  double phiDeg[] = {10, 15, 20};
  double d[n];
  double phi[n];
  double xTrue[n];
  double yTrue[n];
  double xmeas[n];
  double ymeas[n];
  double sigma = 0.1;
  for (int i=0;i<n;i++){
    d[i] = gRandom->Gaus(0,sigma);
    phi[i] = phiDeg[i]*TMath::DegToRad();
    xTrue[i] = z[i]*tTrue*cos(alphaTrue);
    yTrue[i] = z[i]*tTrue*sin(alphaTrue);
    xmeas[i] = xTrue[i]+d[i]*sin(phi[i]); // TODO shift to min radius
    ymeas[i] = yTrue[i]-d[i]*cos(phi[i]); // TODO shift to min radius
  }
  
  for (int i=0;i<n;i++){
    a[i] = +z[i]*sin(phi[i]);
    b[i] = -z[i]*cos(phi[i]);
    c[i] = -xmeas[i]*sin(phi[i]) + ymeas[i]*cos(phi[i]);
  }  
  
  // minimization algorithm  
  double xx[2];
  xx[0] = tTrue;
  xx[1] = alphaTrue;
  printf("%f\n",func(xx));

  ROOT::Math::Functor f(&func,2);
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(10000);  // for GSL
  minimum->SetTolerance(0.001);
  minimum->SetPrintLevel(1);
  minimum->SetFunction(f);
  minimum->SetVariable(0,"tan(theta)",tTrue, 0.01);
  minimum->SetVariable(1,"alpha",alphaTrue, 0.01);
  minimum->Minimize();
  auto xs = minimum->X();
  double tt = xs[0];
  double aa = xs[1];
  double kk = tan(aa);

  // analytical computation
  double A = 0;
  double B = 0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      double ab = (a[i]*b[j]-b[i]*a[j]);
      A+= c[i]*b[j]*ab;
      B+= c[i]*a[j]*ab;
    }
  }
  printf("A=%f\n",A);
  printf("B=%f\n",B);
  double k = - B/A;
  double num = 0;
  double den = 0;
  for (int i=0;i<n;i++){
    double abk = a[i]+b[i]*k;
    num+=c[i]*abk;
    den+=abk*abk;
  }
  double t = -sqrt(1+k*k)*num/den;
  
  printf("kk=%.10f k=%.10f kTrue=%f\n", kk, k, kTrue);
  printf("tt=%.10f t=%.10f tTrue=%f\n", tt, t, tTrue);
}
