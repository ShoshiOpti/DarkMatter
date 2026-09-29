"""Construct the supplied constrained positive mean envelope, including an exterior tail."""
from pathlib import Path
import os,sys,json
os.environ['OPENBLAS_NUM_THREADS']='1'
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'output/scientific_runtime'))
import numpy as np
from scipy.optimize import root
from scipy.sparse import diags,csc_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg import solve
from scipy.interpolate import PchipInterpolator

CORE=HERE/'inputs/core_profile.csv'
def build(n=599,R=48.):
 x=np.arange(1,n+1)*R/(n+1);dr=x[0];g=.3;gp=1.825;gm=-.61875
 raw=np.genfromtxt(CORE,delimiter=',',names=True);N=raw['enclosed_mass_integral'][-1]
 u=np.interp(x,raw['r_dimensionless'],raw['u'],right=0);U=x*u
 T=diags([np.full(n-1,-.5/dr**2),np.full(n,1/dr**2),np.full(n-1,-.5/dr**2)],[-1,0,1]).toarray()
 P=-dr/np.maximum(x[:,None],x[None,:])
 def fun(y):
  z=y[:-1];mu=y[-1];phi=P@(z*z)
  return np.r_[T@z+(phi+g*(z/x)**2-mu)*z,dr*np.sum(z*z)-N]
 def jac(y):
  z=y[:-1];mu=y[-1];phi=P@(z*z)
  j=np.zeros((n+1,n+1));j[:-1,:-1]=T+np.diag(phi+3*g*(z/x)**2-mu)+2*z[:,None]*P*z[None,:]
  j[:-1,-1]=-z;j[-1,:-1]=2*dr*z
  return j
 fit=root(fun,np.r_[U,-.81833],jac=jac,tol=1e-10);res=float(np.max(abs(fun(fit.x))))
 assert res<1e-8
 U=fit.x[:-1];mu=fit.x[-1];u=U/x;phi=P@U**2
 L0=T+np.diag(phi+g*u*u-mu)
 Ls=L0+np.diag((gp-g+gm)*u*u);Lp=L0+np.diag((gp-g-gm)*u*u)
 vals,vec=eigs(csc_matrix(Lp@Ls),k=1,sigma=0,tol=2e-12)
 omega=np.sqrt(vals[0].real);A=vec[:,0].real;B=Ls@A/omega
 fac=np.sqrt(N/(dr*np.sum((A*A+B*B)/2)));A*=fac;B*=fac
 nmode=(A*A+B*B)/2;amode=(A*A-B*B)/2
 mat=np.zeros((n+1,n+1));mat[:-1,:-1]=L0+np.diag(2*g*u*u)+2*U[:,None]*P*U[None,:]
 mat[:-1,-1]=-U;mat[-1,:-1]=2*dr*U
 rhs=np.r_[-U/x**2*(gp*nmode+gm*amode)-U*(P@nmode),-dr*np.sum(nmode)]
 sol=solve(mat,rhs);v=sol[:-1];dw=sol[-1];eps=.1
 positive=(U+eps**2*v)**2+eps**2*nmode
 corr=2*U*v+nmode
 record={'grid':n,'box':R,'dr':dr,'mu':float(mu),'omega':float(omega),'delta_omega':float(dw),
  'N_radial':float(N),'core_residual':res,'mean_linear_residual':float(np.max(abs(mat@sol-rhs))),
  'number_correction_relative':float(abs(dr*np.sum(corr))/N),'epsilon':eps,'r_star':.8,'g0':g,
  'core_central_estimate':float(u[0]),'positive_normalization_before':float(dr*sum(positive)/N),
  'mass_fraction_beyond_r24_before_tail':float(dr*sum(positive[x>24])/ (dr*sum(positive))),
  'scope':'Isolated dimensionless template; per-galaxy baryonic backreaction and shared dimensional microphysics are not included.'}
 return x,U**2,positive,float(np.sqrt(2*(abs(mu)-omega))),record

def main():
 x,core,resp,k,rec=build();x2,c2,r2,k2,rec2=build(1199)
 np.savetxt(HERE/'field_profiles.csv',np.column_stack((x2,c2,r2)),delimiter=',',header='x,core_radial_density,response_radial_density',comments='')
 # Compare normalized cumulative profiles; both include the entire numerical domain.
 def cdf(x,f,z):
  p=PchipInterpolator(np.r_[0,x,48],np.r_[0,f,0]).antiderivative();return p(z)/p(48)
 z=np.geomspace(.2,30,500)
 rec['refinement_max_core_CDF_error']=float(max(abs(cdf(x,core,z)-cdf(x2,c2,z))))
 rec['refinement_max_response_CDF_error']=float(max(abs(cdf(x,resp,z)-cdf(x2,r2,z))))
 (HERE/'profile_diagnostics.json').write_text(json.dumps({'coarse':rec,'fine':rec2,'response_tail_decay_rate':k2,'status':'pass'},indent=2)+'\n')
 print(json.dumps({'coarse':rec,'fine':rec2},indent=2))
if __name__=='__main__':main()
