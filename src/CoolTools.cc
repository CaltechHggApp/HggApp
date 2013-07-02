//
// Written by C. Rogan
// Caltech
// 23-06-08
//
// Big ups to Numerical Recipes in C (minim, jacobi transform etc.)
//

#include "CoolTools.hh"

using namespace std;
using namespace stdcomb;

#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j] = g-s*(h+g*tau); a[k][l] = h + s*(g-h*tau)
#define ITMAX 400
#define EPS 1.0e-14
#define GOLD 1.618034
#define CGOLD 0.3819660
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d)
#define MOV3(a,b,c, d,e,f) (a)=(d); (b)=(e); (c)=(f)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define TOL 2.0e-8



CoolTools::CoolTools(){}

CoolTools::~CoolTools() {}

vector<TLorentzVector> CoolTools::BoostVectors(vector<TLorentzVector> input,
                                               TVector3 b){
  vector<TLorentzVector> output;
  for(int i = 0; i < input.size(); i++){
    TLorentzVector v = input[i];
    v.Boost(b);
    output.push_back(v);
  }

  return output;
}

vector<Jet> CoolTools::BoostJets(vector<Jet> input, TVector3 b){
  vector<Jet> output;

  for(int i = 0; i < input.size(); i++){
    TLorentzVector v = input[i].Get4Vector();
    v.Boost(b);
    output.push_back(Jet(v, input[i].EmFrac(), input[i].HadFrac()));
  }

  return output;
}


//itype 0 - all E 1 - ECAL only 2 - HCAL only  
vector<Jet> CoolTools::CaloTowers2Jets(vector<CaloTower> input, int itype){
  vector<Jet> output;

  for(int i = 0; i < input.size(); i++){
    TLorentzVector v = input[i].Get4Vector();
    
    if(itype == 0)
      output.push_back(Jet(v, input[i].EmFrac(), input[i].HadFrac()));
    if(itype == 1){
      if(input[i].EmFrac() > 0.0)
	output.push_back(Jet(v*input[i].EmFrac(), 1.0, 0.0));
    }
    if(itype == 2){
      if(input[i].HadFrac() > 0.0)
	output.push_back(Jet(v*input[i].HadFrac(), 0.0, 1.0));
    }
  }

  return output;
}

vector<TLorentzVector> CoolTools::Get4Vectors(vector<Jet> input){
  vector<TLorentzVector> output;

  for(int i = 0; i < input.size(); i++){
    output.push_back(input[i].Get4Vector());
  }

  return output;
}

void CoolTools::CalcTSphericity(vector<TLorentzVector> input){
  double stot=0, sxx=0, syy=0, sxy=0;

  double sum_x = 0.0;
  double sum_y = 0.0;

  for(int ij = 0; ij < input.size(); ij++){
    stot += input[ij].Pt();
    sxx += input[ij].Px()*input[ij].Px()/input[ij].Pt();
    sxy += input[ij].Px()*input[ij].Py()/input[ij].Pt();
    syy += input[ij].Py()*input[ij].Py()/input[ij].Pt();
    sum_x += input[ij].Px();
    sum_y += input[ij].Py();
  }
  if (stot==0.) stot = 1.;

  double sph_mtrx[3][3];

  sph_mtrx[0][0] = sxx / stot ;
  sph_mtrx[1][0] = sxy / stot ;
  sph_mtrx[2][0] = 0.0;
  sph_mtrx[0][1] = sxy / stot ;
  sph_mtrx[1][1] = syy / stot ;
  sph_mtrx[2][1] = 0.0;
  sph_mtrx[0][2] = 0.0;
  sph_mtrx[1][2] = 0.0;
  sph_mtrx[2][2] = 0.0;

  double e_vec[3];
  double v[3][3];

  jacobi(sph_mtrx, 2, e_vec, v);

  if(e_vec[0] > e_vec[1]){
    ST_e1 = e_vec[0];
    ST_e2 = e_vec[1];
    ST_v1.SetXYZ(v[0][0], v[1][0], 0.0);
    ST_v2.SetXYZ(v[0][1], v[1][1], 0.0);
  } else {
    ST_e1 = e_vec[1];
    ST_e2 = e_vec[0];
    ST_v1.SetXYZ(v[0][1], v[1][1], 0.0);
    ST_v2.SetXYZ(v[0][0], v[1][0], 0.0);
  }
    
  TVector3 PT;
  PT.SetXYZ(sum_x,sum_y,0.0);
  if(PT.Dot(ST_v1) < 0.0){
    ST_v1 = -1.*ST_v1;
  }
  if(PT.Dot(ST_v2) < 0.0){
    ST_v2 = -1.*ST_v2;
  }

}

void CoolTools::CalcSphericity(vector<TLorentzVector> input){
  double stot=0, sxx=0, syy=0, szz=0, sxy=0, sxz=0, syz=0;

  for(int ij = 0; ij < input.size(); ij++){
    stot += input[ij].P()*input[ij].P();
    sxx += input[ij].Px()*input[ij].Px();
    sxy += input[ij].Px()*input[ij].Py();
    sxz += input[ij].Px()*input[ij].Pz();
    syy += input[ij].Py()*input[ij].Py();
    syz += input[ij].Py()*input[ij].Pz();
    szz += input[ij].Pz()*input[ij].Pz();
  }
  if (stot==0.) stot = 1.;

  double sph_mtrx[3][3];

  sph_mtrx[0][0] = sxx / stot ;
  sph_mtrx[1][0] = sxy / stot ;
  sph_mtrx[2][0] = sxz / stot ;
  sph_mtrx[0][1] = sxy / stot ;
  sph_mtrx[1][1] = syy / stot ;
  sph_mtrx[2][1] = syz / stot ;
  sph_mtrx[0][2] = sxz / stot ;
  sph_mtrx[1][2] = syz / stot ;
  sph_mtrx[2][2] = szz / stot ;

  double e_vec[3];
  double v[3][3];

  jacobi(sph_mtrx, 3, e_vec, v);

  int k,j,i;
  double p;

  for(i = 0; i < 2; i++){
    p=e_vec[k=i];
    for(j = i+1; j < 3; j++)
      if(e_vec[j] >= p) p=e_vec[k=j];
    if(k != i) {
      e_vec[k] = e_vec[i];
      e_vec[i] = p;
      for(j = 0; j < 3; j++){
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }
  }
  S_e1 = e_vec[0];
  S_e2 = e_vec[1];
  S_e3 = e_vec[2];

  S_v1.SetX(v[0][0]);
  S_v1.SetY(v[1][0]);
  S_v1.SetZ(v[2][0]);

  S_v2.SetX(v[0][1]);
  S_v2.SetY(v[1][1]);
  S_v2.SetZ(v[2][1]);

  S_v3.SetX(v[0][2]);
  S_v3.SetY(v[1][2]);
  S_v3.SetZ(v[2][2]);
    
}    

void CoolTools::jacobi(double a[3][3], int n, double d[3], double v[3][3]){

  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;

  double b[3];
  double z[3];

  for(ip=0; ip < n; ip++){
    for(iq = 0; iq < n; iq++) v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }
  for(ip = 0; ip < n; ip++){
    b[ip]=d[ip]=a[ip][ip];
    z[ip] = 0.0;
  }
  for(i = 0; i < 50; i++){
    sm = 0.0;
    for(ip = 0; ip < n-1; ip++){
      for(iq = ip + 1; iq < n; iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0){
      return;
    }
    if(i < 3){
      tresh=0.2*sm/(n*n);
    } else {
      tresh = 0.0;
    }
    for(ip = 0; ip < n-1; ip++){
      for(iq = ip+1; iq < n; iq++){
	g = 100.0*fabs(a[ip][iq]);
	if(i > 3 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	   && (double)(fabs(d[iq])+g) == (double)fabs(d[iq])){
	  a[ip][iq] = 0.0;
	} else {
	  if(fabs(a[ip][iq]) > tresh) {
	    h = d[iq]-d[ip];
	    if((double)(fabs(h)+g) == (double)fabs(h))
	      t = (a[ip][iq])/h;
	    else { 
	      theta = 0.5*h/(a[ip][iq]);
	      t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	      if(theta < 0.0) t = -t;
	    }
	    c = 1.0/sqrt(1+t*t);
	    s = t*c;
	    tau = s/(1.0+c);
	    h = t*a[ip][iq];
	    z[ip] -= h;
	    z[iq] += h;
	    d[ip] -= h;
	    d[iq] += h;
	    a[ip][iq] = 0.0;
	    for(j = 0; j < ip-1; j++){
	      ROTATE(a,j,ip,j,iq);
	    }
	    for(j = ip+1; j < iq-1; j++){
	      ROTATE(a,ip,j,j,iq);
	    }
	    for(j = iq+1; j < n; j++){
	      ROTATE(a,ip,j,iq,j);
	    }
	    for(j = 0; j < n; j++){
	      ROTATE(v,j,ip,j,iq);
	    }
	  }
	}
      }
    }
    for(ip = 0; ip < n; ip++){
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  cout << "Too many iterations in routine jacobi!!!" << endl;
}
  

double CoolTools::FoxWolfram_test(vector<Jet> jets, int l){
  double Etot = 0.0;
  for(int i = 0; i < jets.size(); i++){
    Etot += jets[i].e();
  }
  double H = 0.0;
  for(int i = 0; i < jets.size(); i++){
    for(int j = 0; j < jets.size(); j++){
      //if(i == j) continue;
      TVector3 v1 = jets[i].Get3Vector();
      TVector3 v2 = jets[j].Get3Vector();
      H += v1.Mag()*v2.Mag()*Legendre(v1.Dot(v2)/(v1.Mag()*v2.Mag()), l, 0);
    }
  }
  H /= Etot*Etot;
      
  return H;
}

double CoolTools::FoxWolfram(vector<TLorentzVector> input, int l){
  double Etot = 0.0;
  for(int i = 0; i < input.size(); i++){
    Etot += input[i].E();
  }

  double sum = 0.0;

  for(int im = -l; im <= l; im++){
    double s_real = 0.0;
    double s_imag = 0.0;
    for(int i = 0; i < input.size(); i++){
      double real;
      double imag;
      SphericalHarmonic(input[i].Pz()/input[i].P(), input[i].Phi(), 
			l, im, real, imag);
      s_real += real*input[i].P();
      s_imag += imag*input[i].P();
    }
    sum += s_real*s_real+s_imag*s_imag;
  }
  sum *= 4.0*TMath::Pi()/((2.0*double(l)+1.0)*Etot*Etot);
  return sum;
}

double CoolTools::TranFoxWolfram(vector<TLorentzVector> input, int l){
  double Etot = 0.0;
  for(int i = 0; i < input.size(); i++){
    Etot += input[i].Et();
  }

  double sum = 0.0;

  for(int im = -l; im <= l; im++){
    double s_real = 0.0;
    double s_imag = 0.0;
    for(int i = 0; i < input.size(); i++){
      double real;
      double imag;
      SphericalHarmonic(input[i].Pz()/input[i].P(), input[i].Phi(), 
			l, im, real, imag);
      s_real += real*input[i].Et();
      s_imag += imag*input[i].Et();
    }
    sum += s_real*s_real+s_imag*s_imag;
  }
  sum *= 4.0*TMath::Pi()/((2.0*double(l)+1.0)*Etot*Etot);
  return sum;
}

void CoolTools::SphericalHarmonic(double cos_theta, double phi,int l, int m, 
				  double& real, double& imag){
  double coeff = (2.0*double(l)+1)*double(TMath::Factorial(l-abs(m)));
  coeff /= TMath::Pi()*4.0*double(TMath::Factorial(l+abs(m)));
  coeff = sqrt(coeff)*Legendre(cos_theta, l, m);
  real = coeff*cos(double(m)*phi);
  imag = coeff*sin(double(m)*phi);
  double fac = 1.0;
  if(m >= 0){
    if(m%2 == 1) fac = -1.0;
  }
  real *= fac;
  imag *= fac;
  

}
double CoolTools::Legendre(double x, int l, int m){
  m = abs(m);
  
  if(l == 0){
    return 1.0;
  }
  if(l == 1){
    if(m == 0) return x;
    if(m == 1) return (-1.0*sqrt(1-x*x));
  }
  if(l == 2){
    if(m == 0) return (0.5*(3.0*x*x - 1.0));
    if(m == 1) return (-3.0*x*sqrt(1-x*x));
    if(m == 2) return (3.0*(1.0-x*x));
  }
  if(l == 3){
    if(m == 0) return (0.5*x*(5.0*x*x-3.0));
    if(m == 1) return ((3.0/2.0)*(1.0-5.0*x*x)*sqrt(1-x*x));
    if(m == 2) return (15.0*x*(1.0-x*x));
    if(m == 3) return (-15.0*(1.0-x*x)*sqrt(1.0-x*x));
  }
  if(l == 4){
    if(m == 0) return ((1.0/8.0)*(35.0*x*x*x*x-30.0*x*x+3.0));
    if(m == 1) return ((5.0/2.0)*x*(3.0-7.0*x*x)*sqrt(1.0-x*x));
    if(m == 2) return ((15.0/2.0)*(7.0*x*x-1.0)*(1.0-x*x));
    if(m == 3) return (-105.0*x*(1.0-x*x)*sqrt(1.0-x*x));
    if(m == 4) return (105.0*(1.0-x*x)*(1.0-x*x));
  }
  if(l == 5){
    if(m == 0) return ((1.0/8.0)*x*(63.0*x*x*x*x-70.0*x*x +15.0));
    if(m == 1) return ((-1.0/8.0)*(315.0*x*x*x*x-210.0*x*x+15.0)*sqrt(1.0-x*x));
    if(m == 2) return ((1.0/2.0)*(1.0-x*x)*(315.0*x*x*x-105.0*x));
    if(m == 3) return ((-1.0/2.0)*(1.0-x*x)*(945.0*x*x-105.0)*sqrt(1.0-x*x));
    if(m == 4) return (945.0*x*(1.0-x*x)*(1.0-x*x));
    if(m == 5) return (-945.0*(1.0-x*x)*(1.0-x*x)*sqrt(1.0-x*x));
  }
  if(l == 6){
    if(m == 0) return ((1.0/16.0)*(231.0*x*x*x*x*x*x-315.0*x*x*x*x+105.0*x*x-5.0));
    if(m == 1) return ((-1.0/8.0)*(693.0*x*x*x*x*x-630.0*x*x*x+105.0*x)*sqrt(1.0-x*x));
    if(m == 2) return ((1.0/8.0)*(3465.0*x*x*x*x-1890.0*x*x+105.0)*(1.0-x*x));
    if(m == 3) return ((-1.0/2.0)*(3465.0*x*x*x-945.0*x)*(1.0-x*x)*sqrt(1-x*x));
    if(m == 4) return ((1.0/2.0)*(10395.0*x*x-945.0)*(1.0-x*x)*(1.0-x*x));
    if(m == 5) return ((-10395.0)*x*(1.0-x*x)*(1.0-x*x)*sqrt(1.0-x*x));
    if(m == 6) return (10395.0*(1.0-x*x)*(1.0-x*x)*(1.0-x*x));
  }
  int M = l/2;
  double return_sum = 0.0;
  for(int im = 0; im <= M; im++){
    double fac = double(l-2*im);
    if(fac < double(m)) continue;

    double coeff = 1.0;
    if(im%2 == 1) coeff = -1.0;

    coeff *= double(TMath::Factorial(2*l-2*im));
    coeff /= pow(2.0, double(l));
    coeff /= double(TMath::Factorial(im));
    coeff /= double(TMath::Factorial(l-im));
    coeff /= double(TMath::Factorial(l-2*im));
    
    for(int id = 0; id < m; id++){
      coeff *= fac-double(id);
    }
    if(fac > double(m)) coeff *= pow(x, fac-double(m));
    return_sum += coeff;
  }
  if(m%2 == 1) return_sum *= -1.0;
  return_sum *= pow(1.0-x*x, double(m)/2.0);

  return return_sum;
}

//
// Routines for Thrust calculation
//

void CoolTools::CalcThrust(vector<TLorentzVector> input){

  T_list = input;

  //find maximum momentum input to start
  double max_p = -1.0;
  int max_i = -1;
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    if(T_list[ipart].P() > max_p){
      max_p = T_list[ipart].P();
      max_i = ipart;
    }
  }

  double p[2]; 
  p[0] = T_list[max_i].Phi();
  p[1] = T_list[max_i].Theta();

  bool CONV = true;

  int j, its;
  double gg,gam,fp,dgg, g[2], h[2], xi[2], fret;
  double ftol = 1.0e-8;

  fp = T_func(p);
  T_dfunc(p, xi);
  for(j = 0; j < 2; j++){
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for(its=1; its<=ITMAX;its++){
    fret = T_dlinmin(p,xi);
    if(2.0*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS)) {
      CONV = false;
      break;
    }
    fp = fret;
    T_dfunc(p,xi);
    dgg=gg=0.0;
    for(j = 0; j < 2; j++){
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if(gg == 0.0){
      CONV = false;
      break;
    }
    for(j = 0; j < 2; j++){
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  //if(CONV) cout << "Too many iterations in conjug. grad" << endl;

  T_v1.SetXYZ(cos(p[0])*sin(p[1]),
	   sin(p[0])*sin(p[1]),
	   cos(p[1]));
  double sum_num = 0.0;
  double sum_den = 0.0;
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum_num += fabs(T_list[ipart].Vect().Dot(T_v1));
    sum_den += T_list[ipart].P();
  }
  if(sum_den > 0.0)
    T_1 = sum_num/sum_den;
  else 
    T_1 = 0.0;

  temp_rot_1 = T_v1;
  double phi_min = 0.0;
  phi_min = TT_brent(TOL);

  T_v2 = temp_rot_2;
  T_v2.Rotate(phi_min, temp_rot_1);
  T_v3 = temp_rot_2;
  T_v3.Rotate(phi_min+TMath::Pi()/2.0, temp_rot_1);
  

  sum_num = 0.0;
  sum_den = 0.0;
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum_num += fabs(T_list[ipart].Vect().Dot(T_v2));
    sum_den += T_list[ipart].P();
  }
  if(sum_den > 0.0)
    T_2 = sum_num/sum_den;
  else 
    T_2 = 0.0;

  sum_num = 0.0;
  sum_den = 0.0;
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum_num += fabs(T_list[ipart].Vect().Dot(T_v3));
    sum_den += T_list[ipart].P();
  }
  if(sum_den > 0.0)
    T_3 = sum_num/sum_den;
  else 
    T_3 = 0.0;

  T_list.clear();
}

void CoolTools::CalcTranThrust(vector<TLorentzVector> input){
  T_list = input;

  temp_rot_1.SetXYZ(0.0,0.0,1.0);
  double phi_min = TT_brent(TOL);

  TT_v = temp_rot_2;
  TT_v.Rotate(phi_min, temp_rot_1);

  double sum_num = 0.0;
  double sum_den = 0.0;
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum_num += fabs(T_list[ipart].Vect().Dot(TT_v));
    sum_den += T_list[ipart].Pt();
  }
  if(sum_den > 0.0)
    TT = sum_num/sum_den;
  else 
    TT = 0.0;

  sum_num = 0.0;
  sum_den = 0.0;
  temp_rot_2.Rotate(phi_min+TMath::Pi()/2.0, temp_rot_1);
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum_num += fabs(T_list[ipart].Vect().Dot(temp_rot_2));
    sum_den += T_list[ipart].Pt();
  }
  if(sum_den > 0.0)
    TTm = sum_num/sum_den;
  else 
    TTm = 0.0;

  T_list.clear();
}

double CoolTools::TT_func(double px){
  double sum = 0.0;

  TVector3 n = temp_rot_2;
  n.Rotate(px, temp_rot_1);
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum += fabs(T_list[ipart].Vect().Dot(n));
  }
  if(fabs(sum) == 0.0){
    cout << "sum is zero in TT_func !!!!!!!!" << endl;
    return 10000000000000000.0;
  } else {
    return 1/sum;
  }
}
  
double CoolTools::TT_brent(double tol){
  double xmin;
  double theta = temp_rot_1.Theta();
  double theta_p;
  if(theta >= TMath::Pi()/2.0){
    theta_p = theta-TMath::Pi()/2.0;
  } else {
    theta_p = theta+TMath::Pi()/2.0;
  }
  double phi = temp_rot_1.Phi();
  temp_rot_2.SetXYZ(cos(phi)*sin(theta_p),
		    sin(phi)*sin(theta_p),
		    cos(theta_p));
  double ax,bx,cx,fa,fb,fc;
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  ax = 0.0;
  bx = TMath::Pi();
  TT_mnbrak(&ax, &bx, &cx, &fa, &fb, &fc);
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=TT_func(x);
  for(iter = 1; iter <= ITMAX; iter++){
    xm = 0.5*(a+b);
    tol2 = 2.0*(tol1 = tol*fabs(x)+EPS);
    if(fabs(x-xm) <= (tol2-0.5*(b-a))){
      xmin = x;
      return xmin;
    }
    if(fabs(e) > tol1){
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0*(q-r);
      if(q > 0.0) p = -p;
      q = fabs(q);
      etemp=e;
      e=d;
      if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d = p/q;
	u = x+d;
	if(u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1, xm-x);
      }
    } else {
      d = CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu = TT_func(u);
    if(fu <= fx){
      if(u >= x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if(fu <= fw || w == x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if(fu <= fv || v == x || v == w){
	v=u;
	fv=fu;
      }
    }
  }
  cout << "Too Many iterations in TT_brent!!!!" << endl;
  xmin = x;
  
  return xmin;
}




double CoolTools::T_df1dim(double x){

  int j;
  double df1 = 0.0;
  double xt[2], df[2];
  for(int j = 0; j < 2; j++) xt[j] = T_pcom[j]+x*T_xicom[j];
  T_dfunc(xt, df);
  for(int j = 0; j < 2; j++) df1 += df[j]*T_xicom[j];
  return df1;
}

double CoolTools::T_f1dim(double x){
  int j;
  double f, xt[2];

  for(j = 0; j < 2; j++) xt[j] = T_pcom[j]+x*T_xicom[j];
  f = T_func(xt);
  return f;
}

void CoolTools::T_mnbrak(double *ax, double *bx, double *cx, 
			 double *fa, double *fb, double *fc){
  double ulim, u, r, q, fu, dum;
  *fa = T_f1dim(*ax);
  *fb = T_f1dim(*bx);
  if(*fb > *fa){
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc = T_f1dim(*cx);
  while(*fb > *fc) {
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(max(fabs(q-r),TINY),q-r));
    ulim = (*bx)+GLIMIT*(*cx-*bx);

    if((*bx-u)*(u-*cx) > 0.0){
      fu = T_f1dim(u);
      if(fu < *fc){
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      } else if (fu > *fb){
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx)+GOLD*(*cx-*bx);
      fu = T_f1dim(u);
    } else if ((*cx-u)*(u-ulim) > 0.0){
      fu = T_f1dim(u);
      if(fu < *fc){
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	SHFT(*fb,*fc,fu,T_f1dim(u));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0){
      u = ulim;
      fu = T_f1dim(u);
    } else {
      u = (*cx)+GOLD*(*cx-*bx);
      fu = T_f1dim(u);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}

double CoolTools::T_func(double p[2]){
  double sum = 0.0;

  TVector3 n(cos(p[0])*sin(p[1]),
	     sin(p[0])*sin(p[1]),
	     cos(p[1]));
  for(int ipart = 0; ipart < T_list.size(); ipart++){
    sum += fabs(T_list[ipart].Vect().Dot(n));
  }
  if(fabs(sum) == 0.0){
    cout << "sum is zero in T_func !!!!!!!!" << endl;
    return 10000000000000000.0;
  } else {
    return 1/sum;
  }
}

void CoolTools::T_dfunc(double x[2], double dfx[2]){
  dfx[0] = 0.0;
  dfx[1] = 0.0;

  TVector3 n(cos(x[0])*sin(x[1]),
	     sin(x[0])*sin(x[1]),
	     cos(x[1]));

  double phi = x[0];
  double the = x[1];

  for(int ipart = 0; ipart < T_list.size(); ipart++){
    double px = T_list[ipart].Px();
    double py = T_list[ipart].Py();
    double pz = T_list[ipart].Pz();

    double den = fabs(T_list[ipart].Vect().Dot(n));
    double num0 = sin(the)*sin(the)*sin(phi)*cos(phi)*(py*py-px*px);
    double num1 = sin(the)*cos(the)*(cos(phi)*cos(phi)*px*px+
				    sin(phi)*sin(phi)*py*py-pz*pz);
    if(den > 0.0){
      dfx[0] += num0/den;
      dfx[1] += num1/den;
    }
  }
  double t = T_func(x);
  dfx[0] *= -t*t;
  dfx[1] *= -t*t;

  return;
}

double CoolTools::T_dlinmin(double p[2], double xi[2]){
  int j;
  double xx,xmin,fx,fb,fa,bx,ax,fret;

  for(j = 0; j < 2; j++){
    T_pcom[j] = p[j];
    T_xicom[j] = xi[j];
  }

  ax = 0.0;
  xx = 1.0;
  
  T_mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
  fret = T_dbrent(ax, xx, bx, TOL, &xmin);

  for(j = 0; j < 2; j++){
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  return fret;
}

double CoolTools::T_dbrent(double ax, double bx, double cx, 
			   double tol, double *xmin){
  
  int iter, ok1, ok2;
  double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=T_f1dim(x);
  dw=dv=dx=T_df1dim(x);
  for(iter=1; iter <= ITMAX; iter++){
    xm = 0.5*(a+b);
    tol1=tol*fabs(x)+EPS;
    tol2=2.0*tol1;
    if(fabs(x-xm) <= (tol2-0.5*(b-a))){
      *xmin = x;
      return fx;
    }
    if(fabs(e) > tol1){
      d1 = 2.0*(b-a);
      d2=d1;
      if(dw != dx) d1=(w-x)*dx/(dx-dw);
      if(dv != dx) d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde = e;
      e=d;
      if(ok1 || ok2){
	if(ok1 && ok2)
	  d = (fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d = d1;
	else 
	  d = d2;
	if(fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if(u-a < tol2 || b-u < tol2)
	    d = SIGN(tol1, xm-x);
	} else {
	  d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if(fabs(d) >= tol1){
      u=x+d;
      fu = T_f1dim(u);
    } else {
      u = x+SIGN(tol1,d);
      fu = T_f1dim(u);
      if(fu > fx){
	*xmin = x;
	return fx;
      }
    }

    du = T_df1dim(u);
    if(fu <= fx){
      if(u >= x) a = x; else b = x;
      MOV3(v,fv,dv, w,fw,dw);
      MOV3(w,fw,dw, x,fx,dx);
      MOV3(x,fx,dx, u,fu,du);
    } else {
      if(u < x) a = u; else b = u;
      if(fu <= fw || w == x){
	MOV3(v,fv,dv, w,fw,dw);
	MOV3(w,fw,dw, u,fu,du);
      } else  if(fu < fv || v == x || v == w){
	MOV3(v,fv,dv, u,fu,du);
      }
    }
  }
  cout << "Too many iterations in routint dbrent!!!" << endl;
  return 0.0;
}

void CoolTools::TT_mnbrak(double *ax, double *bx, double *cx, 
			  double *fa, double *fb, double *fc){
  double ulim, u, r, q, fu, dum;
  *fa = TT_func(*ax);
  *fb = TT_func(*bx);
  if(*fb > *fa){
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc = TT_func(*cx);
  while(*fb > *fc) {
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(max(fabs(q-r),TINY),q-r));
    ulim = (*bx)+GLIMIT*(*cx-*bx);

    if((*bx-u)*(u-*cx) > 0.0){
      fu = TT_func(u);
      if(fu < *fc){
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      } else if (fu > *fb){
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx)+GOLD*(*cx-*bx);
      fu = TT_func(u);
    } else if ((*cx-u)*(u-ulim) > 0.0){
      fu = TT_func(u);
      if(fu < *fc){
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	SHFT(*fb,*fc,fu,TT_func(u));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0){
      u = ulim;
      fu = TT_func(u);
    } else {
      u = (*cx)+GOLD*(*cx-*bx);
      fu = TT_func(u);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}

//////////////// CLEO CONES //////////////////////

/// 3D CLEO cones w.r.t the Sphericity Axis
void CoolTools::MakeCLEOCones_SphAxis(int ncons, vector<TLorentzVector> input) {
  MakeCLEOCones(ncons, input, S_v1);
}

/// 3D CLEO cones w.r.t the Thrust Axis
void CoolTools::MakeCLEOCones_ThrAxis(int ncons, vector<TLorentzVector> input) {
  MakeCLEOCones(ncons, input, T_v1);
}

/// Transverse CLEO cones w.r.t the transverse Sphericity Axis
void CoolTools::MakeTCLEOCones_SphAxis(int ncons, vector<TLorentzVector> input) {
  MakeCLEOCones(ncons, input, ST_v1);
}

/// Transverse CLEO cones w.r.t the transverse Thrust Axis
void CoolTools::MakeTCLEOCones_ThrAxis(int ncons, vector<TLorentzVector> input) {
  MakeCLEOCones(ncons, input, TT_v);
}

/// Transverse CLEO cones w.r.t the pT component of a generic axis
void CoolTools::MakeTCLEOCones(int ncons, vector<TLorentzVector> input, TVector3 axis) {
  
  vector<TLorentzVector> input_T;
  for(int i=0; i<int(input.size()); i++) 
    input_T.push_back(TLorentzVector(input[i].X(), input[i].Y(), 0., input[i].T()));
  
  TVector3 axis_T(axis.X(), axis.Y(), 0.);
  
  MakeCLEOCones(ncons, input_T, axis_T);
}
 
/// 3D CLEO cones w.r.t a generic axis
void CoolTools::MakeCLEOCones(int ncons, vector<TLorentzVector> input, TVector3 axis) {

  const double piover2 = asin(1);
  const double pi = 2.* piover2;
  const double angwidth = piover2/ncons;

  // delete previously calculated Energy flows
  _eflowarray.clear();
  _momsum = 0.;

  for(int i=0; i<ncons; i++) 
    _eflowarray.push_back(0.);
  
  // loop over the all candidates
  for(int i = 0; i < input.size(); i++){

    // this candidate's 3-momentum   
    TVector3 p = input[i].Vect();
    double pmag = p.Mag();
    
    // the angle with the event axis
    double ang =  axis.Angle(p);
    if( ang > piover2 ) ang = pi - ang;
    int l = (int)(ang/angwidth);
    if( l < 0 ) {
      cout << "Error in MakeCLEOCones: negative cone number. Candidate ignored " << endl;
    } else if( l >= ncons ) {
      cout << "Error in MakeCLEOCones: too big cone number. Candidate ignored  " << endl;
    } else {
      // contribution of this track
      _eflowarray[l] += pmag;
      // total momentum
      _momsum += pmag;
    }
  }
}

/// Get The momentum in CLEO cone # ncon
double CoolTools::GetMomCLEOCones(int ncon) {
  double eflow;
  if(ncon >= int(_eflowarray.size())) {
    cout << "Error in GetMomCLEOCones: too big cone number. " << endl;
    eflow = -999.;
  } else if(ncon < 0) {
    cout << "Error in GetMomCLEOCones: negative cone number. " << endl;
    eflow = -999.;
  } else {
    eflow = _eflowarray[ncon];
  }
  return eflow;
}

/// Get the total momentum flow
double CoolTools::GetMomFlowCLEOCones() {
  return _momsum;
}
