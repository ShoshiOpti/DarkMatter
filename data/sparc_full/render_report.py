"""Render the complete SPARC atlas and a standalone LaTeX results report."""
from pathlib import Path
import os,sys,json,csv,math,shutil,subprocess
import test_sparc as t
ROOT=t.ROOT;HERE=t.HERE;OUT=(ROOT/'output/pdf/sparc_full_profile_2026_09_29') if (ROOT/'Cubic_Scalar_GR_Intersection.tex').exists() else HERE.parent/'report';OUT.mkdir(parents=True,exist_ok=True)
FIG=OUT/'figures';FIG.mkdir(exist_ok=True);TMP=ROOT/'tmp/pdfs/sparc_full_profile_2026_09_29';TMP.mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(TMP/'mpl_cache')
sys.path.insert(0,str(ROOT/'research/core_contact_bridge_2026_09_28/plot_runtime'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
data=t.load();rows=list(csv.DictReader((HERE/'predictions.csv').open()));idx={}
for r in rows:idx.setdefault((r['galaxy'],r['phase'],r['model']),[]).append(float(r['predicted']))
colors={'baryons':'#737a82','scalar':'#326b91','new_analytic':'#af3370','mean_response':'#bd782a','all_new':'#288069'}
labels={'baryons':'Baryons only','scalar':'Scalar halo (inner fit)','new_analytic':'New analytic envelope (inner fit)','mean_response':'Mean response (inner fit)','all_new':'New analytic envelope (all rows fitted)'}
styles={'baryons':'--','scalar':':','new_analytic':'-','mean_response':'-.','all_new':'--'}

def sheet(names,title,atlas=False):
 fig,axes=plt.subplots(3,2,figsize=(10.8,11.5));fig.subplots_adjust(left=.08,right=.98,bottom=.12,top=.925,hspace=.36,wspace=.23)
 mods=['baryons','new_analytic','mean_response','all_new'] if atlas else ['baryons','scalar','new_analytic','mean_response','all_new']
 for ax,n in zip(axes.flat,names):
  d=data[n];r=d['r'];edge=r[~d['outer']][-1];ax.axvspan(edge,r[-1]*1.025,color='#eef0f2',zorder=0)
  ymax=max(d['y']+d['err']);ymin=min(0.,min(d['y']-d['err']))
  for mod in mods:
   phase='all' if mod=='all_new' else 'inner';m='new_analytic' if mod=='all_new' else mod;v=idx[n,phase,m]
   ax.plot(r,v,color=colors[mod],ls=styles[mod],lw=1.6 if mod=='new_analytic' else 1.1);ymax=max(ymax,max(v))
  for outer in [False,True]:
   mask=d['outer']==outer;ax.errorbar(r[mask],d['y'][mask],d['err'][mask],fmt='o',ms=3,capsize=1.3,color='#1b2730',mfc='white' if outer else '#1b2730',elinewidth=.6,zorder=8)
  flag='' if d['quality'] else ' [outside original quality cut]'
  ax.set(title=n+flag,xlabel='Radius [kpc]',ylabel='Speed [km/s]',xlim=(0,r[-1]*1.025),ylim=(ymin*1.1,ymax*1.08));ax.title.set_fontsize(9);ax.grid(axis='y',alpha=.17)
 for ax in list(axes.flat)[len(names):]:ax.set_visible(False)
 fig.suptitle(title,x=.08,ha='left',fontsize=15,fontweight='bold')
 handles=[Line2D([0],[0],color=colors[m],ls=styles[m],lw=1.7,label=labels[m]) for m in mods]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.047),ncol=2,frameon=False,fontsize=8)
 fig.text(.5,.029,'Open points / shading: outer evaluation block. Solid profiles retain their full exterior mass.',ha='center',fontsize=8)
 fig.text(.5,.014,'Catalogue D, inclination and stellar M/L fixed. Green full-data fit is descriptive, not held-out evidence.',ha='center',fontsize=8)
 return fig

names=['DDO154','UGC11557','F568-3','NGC4559','UGC06614','NGC2955']
f=sheet(names,'Figure 3 revisited: complete-mass profiles');f.set_size_inches(8.2,10.5);f.subplots_adjust(bottom=.17);f.savefig(FIG/'figure3.pdf');f.savefig(HERE/'figure3_updated.png',dpi=180);plt.close(f)
with PdfPages(OUT/'SPARC_all_175_atlas.pdf') as pdf:
 for i in range(0,len(data),6):
  f=sheet(list(data)[i:i+6],f'SPARC full-profile atlas | {i+1}-{min(i+6,175)} of 175',True);pdf.savefig(f);plt.close(f)
res=json.loads((HERE/'RESULTS.json').read_text());nu=json.loads((HERE/'nuisance_sensitivity.json').read_text());ver=json.loads((HERE/'verification.json').read_text());guard=json.loads((HERE/'guard_sensitivity.json').read_text())
scores=list(csv.DictReader((HERE/'galaxy_scores.csv').open()));look={(x['galaxy'],x['phase'],x['model']):float(x['loss']) for x in scores}
states=list(csv.DictReader((HERE/'fit_states.csv').open()))
fig,axs=plt.subplots(1,2,figsize=(10.8,4.5),constrained_layout=True)
xx=np.array([look[n,'inner','baryons'] for n in data]);yy=np.array([look[n,'inner','new_analytic'] for n in data]);good=np.array([data[n]['quality'] for n in data])
axs[0].scatter(xx[good],yy[good],s=19,c='#af3370',alpha=.75,label='Original 131');axs[0].scatter(xx[~good],yy[~good],s=19,facecolors='none',edgecolors='#555',label='Additional 44')
axs[0].plot([1e-5,1e5],[1e-5,1e5],color='gray',ls=':');axs[0].set(xscale='log',yscale='log',xlabel='Baryon outer mean standardized loss',ylabel='New-envelope outer loss',title='Each point is one galaxy');axs[0].legend(fontsize=8)
outside=[float(s['mass_outside_last_observation_fraction']) for s in states if s['phase']=='inner' and s['model']=='new_analytic']
axs[1].hist(outside,bins=np.linspace(0,1,21),color='#288069',alpha=.85);axs[1].set(xlabel=r'$1-M(<R_{\max})/M_{\rm total}$',ylabel='Number of galaxies',title='Exterior mass retained, not discarded');axs[1].axvline(np.median(outside),color='black',ls=':',label=f'Median {np.median(outside):.3f}');axs[1].legend(fontsize=8)
fig.savefig(FIG/'summary.pdf');fig.savefig(HERE/'full_profile_summary.png',dpi=180);plt.close(fig)

def fmt(x,d=4):return f'{x:.{d}g}'
def table_summary(phase):
 out=[]
 for sample,lab in [('all175','175 galaxies'),('quality131','131 quality galaxies')]:
  for x in res['summaries'][sample]:
   if x['phase']==phase and x['model'] in ['baryons','scalar','new_analytic','mean_response']:
    name={'baryons':'Baryons only','scalar':'Scalar halo','new_analytic':'New analytic envelope','mean_response':'Positive mean response'}[x['model']]
    out.append(f"{lab} & {name} & {x['mean_loss']:.4f} & {x['mae']:.3f} \\\\")
 return '\n'.join(out)
def statrows(sample):
 out=[]
 for x in res['inference'][sample]:
  name=('Analytic' if x['model']=='new_analytic' else 'Mean response')+' / '+('baryons' if x['baseline']=='baryons' else 'scalar')
  out.append(f"{name} & {x['mean_difference']:.4f} & ${x['Holm_p']:.5g}$ & {x['Holm_Z']:.3f} \\\\")
 return '\n'.join(out)
tex=r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=23mm]{geometry}
\usepackage[T1]{fontenc}\usepackage{lmodern,microtype,amsmath,amssymb,booktabs,longtable,graphicx,hyperref}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\allowdisplaybreaks\emergencystretch=2em
\title{Full SPARC tests of complete-mass\\two-field halo-envelope approximations}
\author{Reproducible project calculation}\date{29 September 2026}
\begin{document}\maketitle
\begin{abstract}
We test the supplied positive mean-response construction and second-order analytic envelope on all 175 SPARC galaxies and 3,391 rotation measurements. Halo mass profiles extend to infinity; signed baryonic force templates are retained. We report complete-curve fits and a separate inner-fit/outer-evaluation procedure containing 761 outer points, including the historical quality subset of 131 galaxies and 659 outer points. With catalogue calibration fixed, the new analytic envelope reduces equal-galaxy outer loss from 461.1992 to 44.9633. Its grouped bootstrap comparison reaches the simulation resolution floor: the reported four-comparison Holm-adjusted Monte Carlo value is $p=4\times10^{-5}$, equivalent to two-sided $Z=4.107$. When distance, inclination and stellar mass-to-light ratio are profiled with common priors, losses become 81.9405 and 25.5581, with exploratory $p=0.01083$, $Z=2.548$. The mean-response improvement over a fitted scalar halo is much smaller than the improvement over baryons. These are retrospective conditional template tests, not a fit with universal microscopic parameters or a detection of uniquely two-field physics.
\end{abstract}
\tableofcontents

\section{What is being tested}
The user-supplied \texttt{two\_field\_friedmann\_extension.tex} provides the controlling equations. We independently reconstruct its constrained positive mode envelope and implement its analytic core-plus-envelope prescription. The homogeneous Friedmann solution does not provide a galaxy's radial density or occupation. No term proportional to a homogeneous expansion rate is added to the local rotation force.

The new analytic profile uses fixed, previously theory-compressed mode-0 shape constants. Each galaxy fits only total halo mass and length. The scalar and mean-response comparators fit the same two coordinates. The mean response has a supplied excitation amplitude $\varepsilon=0.1$, not an occupation estimated from an outer rotation curve. No galaxy is dropped because its fit is poor.

These are isolated spherical halo \emph{templates} added to SPARC baryonic forces. They do not solve a separate baryon-forced two-field equilibrium in each galaxy, nor impose common $m,M_\phi,\beta_*$. The data can test these declared observable profiles without supplying those missing physical constraints. All halo models here still use GR in its weak-field limit; the null is \emph{baryon-only} GR, not GR with a dark halo.

\section{Equations and the entire mass profile}
\subsection{Catalogue forward prediction}
Let $F_X=V_X|V_X|$, preserving outward gas-force contributions. With $d=D/D_0$, $t_i=\sin i/\sin i_0$, $R=dR_{\rm cat}$, the prediction is
\begin{equation}
 V_{\rm cat}^2=d t_i^2(F_{\rm gas}+\Upsilon_dF_{\rm disk}+\Upsilon_bF_{\rm bulge})
       +t_i^2\frac{GM_h}{R}\,F_h(R/L).
\end{equation}
The primary catalogue comparison fixes $d=t_i=1$, $\Upsilon_d=0.5$, $\Upsilon_b=0.7$. SPARC's gas template already includes helium. Disk and gas force templates are not replaced by spherical enclosed baryon masses; their signed geometry is retained \cite{sparc}.

For the analytic envelope, use
\begin{align}
 F_h(x)&=(1-f_e)F_c(x/s_c)+f_e F_e^{(2)}(x/s_e),\\
 F_e^{(2)}(z)&=\frac{2}{\pi}\left(\arctan z-\frac{z}{1+z^2}\right)
       +\frac{8q}{3\pi}\frac{z^3}{(1+z^2)^3},\\
 (s_c,f_e,s_e,q)&\simeq(1.105453894,0.050356819,1.488431634,0.00059489).
\end{align}
The last constants are taken at full precision from the archived theory-compression record; no new SPARC outcome selects them. $F_c$ is the recalculated unit-mass scalar profile. The corresponding envelope density is
\begin{equation}
 \rho_e(R)=\frac{M_e}{\pi^2s^3}
 \left[\frac{1}{(1+z^2)^2}+\frac{2q(1-z^2)}{(1+z^2)^4}\right],\qquad z=R/s.
\end{equation}
It is positive for $0\le q<1$, integrates to $M_e$, and has an $R^{-4}$ tail. The second-order term changes the distribution, not the total mass. The exact old dPIE is also refitted under precisely the same present protocol, using a stable divided-difference integral.

\subsection{Positive mode response with compensating core backreaction}
On the supplied $r_*=0.8$, $g_0=0.3$ core, write $L_0u=0$ and
\begin{align}
 L_s&=L_0+(g_+-g_0+g_-)u^2,&L_p&=L_0+(g_+-g_0-g_-)u^2,\\
 L_pB&=\Omega A,&L_sA&=\Omega B,
\end{align}
where $g_+=1.825$, $g_-=-0.61875$. Normalize
$\int(A^2+B^2)/2=N$ and set $\mathcal N=(A^2+B^2)/2$, $\mathcal A=(A^2-B^2)/2$.
We solve the supplied constrained mean problem
\begin{align}
 (L_0+2g_0u^2)v+mu\overline\Phi_2-\delta\omega\,u
   &=-u(g_+\mathcal N+g_-\mathcal A),\\
 \nabla^2\overline\Phi_2&=4\pi Gm(2uv+\mathcal N),\\
 \int(2uv+\mathcal N)\,d^3x&=0.
\end{align}
The fitted template is obtained from
\begin{equation}
 n_+(R)=\mathcal Z[(u+\varepsilon^2v)^2+\varepsilon^2\mathcal N],\qquad
 \int n_+\,d^3x=N,\quad \varepsilon=0.1.
\end{equation}
Both the core response and its Poisson feedback are retained. This is a second-order mean-envelope approximation, not an exact eternally stationary excited halo. Changing its overall fitted mass between galaxies is a separate state fit; the internal correction itself conserves number.

\subsection{Exterior mass is retained, but shell theorem still applies}
For every spherical profile,
\begin{equation}
 \Phi_h(R)=-G\left[\frac{M(<R)}R+4\pi\int_R^\infty\rho_h(s)s\,ds\right],
 \qquad R\Phi_h'(R)=\frac{GM(<R)}R.
\end{equation}
Exterior shells contribute to potential depth but cancel from the interior radial force. Replacing $M(<R)$ by total $M_h$ would incorrectly count exterior spherical mass as inward acceleration. No profile is cut off at the outer-data boundary or the final observed radius. The weak-field guard uses the \emph{full} central potential, including exterior mass.

For the recalculated wave profiles, a positive exponential density continuation is matched at dimensionless radius 30: $\rho(x)=\rho(30)(30/x)^2e^{-2\kappa(x-30)}$. It is integrated analytically to infinity. The scalar uses $\kappa=\sqrt{2|\omega_c|}$; the mode response uses $\sqrt{2(|\omega_c|-\Omega)}$. This is a declared asymptotic continuation, not a new measured exterior density. Its added mass fractions are $2.20\times10^{-27}$ and $2.07\times10^{-22}$. The analytic envelope's much larger algebraic tail is retained exactly. Unknown exterior baryons are not invented beyond the supplied SPARC mass model.

\section{Observational and inference protocol}
The complete local SPARC tables contain 175 galaxies and 3,391 rows. The historical cut $Q<3$, $i\ge30^\circ$, $N_g\ge8$ retains 131 galaxies and 3,034 rows. The other 44 remain in the full-sample analysis and atlas, explicitly flagged. They are not promoted to high-quality circular-motion measurements.

Two fits are supplied. The all-row fit uses the entire observed curve and is descriptive. The radial-prediction fit uses the inner $N_g-\max(2,\lceil0.2N_g\rceil)$ rows and evaluates every outer row. This gives 761 outer rows for all 175 galaxies and 659 for the quality subset. Holding observations out of optimization does not remove any outer mass from the force model. Some short curves have only two inner points for two fitted halo parameters and are weakly identifying; this motivates the separate quality result.

The length scan covers $R_{\min}/100\le L\le100R_{\max}$. At each length the nonnegative squared-amplitude fit is convex and solved by its monotone derivative, with positivity at every supplied radius. The guard is $|\Phi_h(0)|/c^2\le10^{-3}$. A 65-point logarithmic scan and local refinement search the length coordinate. A doubled scan on the six historical examples and bound-active cases changes predicted speeds by at most $2.69\times10^{-5}\,\mathrm{km\,s^{-1}}$. Bound hits are retained and reported.

\paragraph{A baseline failure that must remain visible.}
UGC01281 has two central catalogue rows with negative net baryonic force at the adopted stellar masses. A circular speed does not exist for that baryon-only prediction. The all-row descriptive score uses a flagged $V_b=\sqrt{\max(V_b^2,0)}$ convention there. Both rows lie in the inner block, so no outer fixed-calibration loss uses this substitution. The galaxy and all its points remain present. Separate admissible-subset sensitivity records are supplied; no failed point is silently removed. Calibration-profile sensitivities using the zero-speed convention still do not repair a globally inadmissible baryonic equilibrium.

For each galaxy, the outer loss is
\begin{equation}
 \ell_g=\frac1{N_{g,\rm out}}\sum_{i\in\rm out}
       \left(\frac{V_{{\rm pred},i}-V_{{\rm obs},i}}{\sigma_i}\right)^2,
 \qquad\Delta=\frac1{N_g}\sum_g(\ell_{g,\rm new}-\ell_{g,\rm base}).
\end{equation}
These are equal-galaxy means, not reduced chi-squares. Negative $\Delta$ favors the new profile. We resample calibration groups, keeping all distance-method-4 Ursa Major galaxies together; all other galaxies form individual groups. There are 148 groups in the full sample and 111 in the quality sample.

From $B=99,999$ group-pairs draws, we report the approximate null-centered two-sided bootstrap value
\begin{equation}
 p_{\rm MC}=\frac{1+\#\{|\Delta^*-\Delta|\ge|\Delta|\}}{B+1},
 \qquad Z_{\rm eq}=\Phi^{-1}(1-p/2).
\end{equation}
$Z_{\rm eq}$ is an unsigned Gaussian-equivalent conversion; it is not $\sqrt{\Delta\chi^2}$. Four primary comparisons (two new profiles against baryons and against the scalar) receive Holm adjustment within each sample. Exact two-sided group-sign diagnostics are supplementary: they test frequency of improvement, not mean magnitude. They additionally assume independent groups with fair null signs.

The public velocity errors omit inclination systematics, and a full radial covariance matrix is unavailable. Group resampling does not manufacture that covariance. The bootstrap and sign values are conditional approximations for this non-random, historically inspected sample and model family; they do not account for the project's historical search or establish a universal discovery significance.

\section{Results on all radial observations and outer predictions}
\subsection{Complete-curve fits: descriptive accuracy}
\begin{center}\small
\begin{tabular}{llrr}\toprule
Sample & Model & Mean all-row loss & MAE (km/s)\\\midrule
@@ALLTABLE@@
\bottomrule\end{tabular}\end{center}
These optimized all-row errors are not assigned predictive p-values. Adding a halo is expected to improve fitted residuals; the radial evaluation below is the more demanding comparison.

\subsection{Inner-fit, outer-evaluation results}
\begin{center}\small
\begin{tabular}{llrr}\toprule
Sample & Model & Mean outer loss & MAE (km/s)\\\midrule
@@OUTTABLE@@
\bottomrule\end{tabular}\end{center}
For all 175 galaxies, the analytic profile improves outer loss in 167 galaxies, worsens it in five, and has three numerical ties. The positive mean response improves 168, worsens four, and has three ties. Small average differences between halo models should not be confused with the much larger baryon-to-halo change.

\begin{figure}[p]\centering\includegraphics[width=\linewidth]{figures/figure3.pdf}
\caption{The same six historical Figure 3 examples under the \emph{new matched catalogue calibration and fitting protocol}. Open points are outer evaluation rows. The green all-row fit is supplied to distinguish curve fitting from radial prediction. All 175 galaxies appear in the companion atlas. No line is represented as a cosmological state prediction.}\end{figure}

\subsection{Grouped significance and equivalent Gaussian significance}
\begin{center}\small
\begin{tabular}{lrrr}\toprule
All 175: comparison & Mean difference & Holm $p_{\rm MC}$ & $Z_{\rm eq}$\\\midrule
@@STATS175@@
\bottomrule\end{tabular}\end{center}
\begin{center}\small
\begin{tabular}{lrrr}\toprule
Quality 131: comparison & Mean difference & Holm $p_{\rm MC}$ & $Z_{\rm eq}$\\\midrule
@@STATS131@@
\bottomrule\end{tabular}\end{center}
For the analytic-versus-baryon all-sample test, zero bootstrap draws exceed the observed centered statistic. The reported $p=4\times10^{-5}$ is therefore a multiplicity-adjusted Monte Carlo floor, not a measured tail probability of exactly that size and not permission to extrapolate to arbitrarily large $Z$. The quality-sample comparison has one exceedance before adjustment.

The analytic-versus-baryon mean contrast has approximate 95\% group-resampling interval $[-599.363,-277.933]$ for all 175 and $[-721.530,-333.020]$ for the quality subset. Removing the baseline-inadmissible galaxy only as a separately labelled sensitivity leaves the direction and interval conclusion intact; adjusted Monte Carlo values are $8\times10^{-5}$ and $2\times10^{-4}$, respectively.

The supplementary adjusted group-sign values for analytic versus baryons are $p=3.60\times10^{-35}$, $Z=12.374$ (141 improving versus five worsening non-tied groups) and $p=2.82\times10^{-29}$, $Z=11.233$ (108 versus two) for the two samples. These answer a different question: improvement is widespread. They must not replace the mean-loss p-values or be advertised as a 12-sigma detection of the action, its cubic, or dark-field occupation.

\subsection{Matched nuisance profiling reduces the evidence strength}
As an additional exploratory check, both the baryon and analytic-envelope models fit the same distance, inclination and common stellar multiplier on inner rows. Gaussian priors use catalogue distance and inclination errors; the common stellar multiplier has a 0.1-dex log prior around one. Disk and bulge reference masses scale together. Bounds are $0.05\le D/D_0\le5$, $5^\circ\le i\le90^\circ$, and $0.1\le\mathrm{stellar\ multiplier}\le4$. Both models use the same values and bounds. This is profiling, not Bayesian evidence or marginalization.
\begin{center}\small
\begin{tabular}{lrrrr}\toprule
Sample & Baryon loss & Envelope loss & $p_{\rm MC}$ & $Z_{\rm eq}$\\\midrule
All 175 & 81.9405 & 25.5581 & 0.01083 & 2.548\\
Quality 131 & 101.9967 & 31.8958 & 0.01453 & 2.444\\\bottomrule
\end{tabular}\end{center}
These exploratory sensitivity p-values are unadjusted for the extra checks. All 350 nuisance fits report optimizer success and satisfy the declared halo guard. The all-sample mean difference is $-56.3825$, with interval $[-102.056,-27.349]$. Profiling only the stellar multiplier gives intermediate losses 223.3679 versus 26.3115 and $p=0.00038$. Therefore a single calibration-independent significance number is not justified.

\section{Outer mass, extrapolation and actual improvement}
\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{figures/summary.pdf}
\caption{Left: held-out loss for every galaxy. Right: fitted halo mass outside its final measured radius, retained in the profile. Unobserved mass is a model extrapolation, not an independently measured inventory.}\end{figure}
The median analytic-envelope mass fraction outside $R_{\max}$ is 0.02693; its largest fitted value is 0.999999905. The median fitted total mass is $1.87\times10^{10}M_\odot$, but one extrapolation reaches $2.81\times10^{16}M_\odot$. The latter is not an acceptable measured galaxy mass: the inner curve can constrain a central force combination while leaving the exterior extent very weakly identified. The primary analytic fit has six length-bound hits and ten potential-guard hits. Neither deleting outer mass nor using its total as enclosed mass repairs that degeneracy. Independent lensing, satellites, environmental tides or explicit abundance/size priors are needed.

Tightening the full central-potential guard from $10^{-3}$ to $10^{-4}$ and $10^{-5}$ gives full-sample outer losses 44.8524 and 43.8148, compared with 44.9633 initially. Thirteen fits then hit the tighter guard. These are sensitivity results; the primary model is not reselected using the favorable outer scores.

The new analytic and old exact dPIE prescriptions give outer losses 44.963256 and 44.963251 under the same current fits. Their maximum relative composite mass-profile difference is $7.27\times10^{-8}$. The improved formula is a stable analytic approximation, not the cause of a new observational gain. The mean-response profile improves on the scalar by only 0.20468 mean-loss units (about 0.447\%) in the all-sample test. The analytic mixture improves by 0.87570 (about 1.91\%). The larger advantage against baryons is principally evidence that an additional fitted gravitating profile describes these data better than baryons alone.

Historical loss 18.8200 used inherited kinematic/nuisance calibrations and a different length-fitting procedure. It must not be subtracted from the present 44.9633 or 58.0592 to attribute a change to new dynamics. The directly refitted old-exact row is the matched formula comparison.

Independent per-galaxy mass and length fits also do not impose $Lv_s=\hbar/m$ with the same particle mass across the sample. A universal-microphysics test must recompute baryon-forced states, enforce shared physical scales, supply occupations and boundary conditions, and then repeat the prediction and statistical pipeline. That is a stronger test than the completed full-SPARC envelope test here.

\section{Numerical checks and reproducibility}
The mode calculation solves a discrete stationary core at fixed total number, the lowest internal bound mode, and the constrained mean backreaction on $0<x<48$. Grids with 599 and 1,199 interior points have spacings 0.08 and 0.04. The fine solution has $\omega_c=-0.8184853845$, $\Omega=0.2699598986$, and $\delta\omega=0.5353323128$. The core and mean linear residuals are below $1.9\times10^{-13}$ and $1.1\times10^{-12}$; the second-order number correction is below $4.8\times10^{-15}$ of the total. Grid refinement changes the normalized scalar and response CDFs by at most $5.56\times10^{-4}$ and $5.48\times10^{-4}$ on the tested radial interval. Small algebraic residuals are not continuum-error certificates.

Refitting all 175 galaxies at the coarser field resolution changes the mean response-minus-scalar outer contrast from $-0.204684$ to $-0.205032$, a difference of $0.000348$ (about 0.17\% of the effect). For the quality subset the change is $0.000452$. This checks the numerical stability of the small profile correction; it does not address physical occupation uncertainty.

Independent quadratures verify the analytic density normalization and its cumulative mass for $q=0,0.00059489359,0.1,0.8$ to below $2.3\times10^{-16}$. Recomputing every halo prediction as $GMF(R/L)/R$ reproduces the saved velocities to $2.1\times10^{-13}\,\mathrm{km\,s^{-1}}$. All profiles are monotone cumulative masses. Input SHA-256 hashes, all fitted parameters, all observed and predicted rows, sample flags, nuisance fits, bound checks, p/Z records and scripts are included.

From the reproduction directory, run
\begin{verbatim}
python -B build_profiles.py
python -B test_sparc.py
python -B robustness.py
python -B audit.py
python -B check_mesh.py
python -B render_report.py
\end{verbatim}
The scripts require NumPy, SciPy and Matplotlib; compilation uses pdfLaTeX. The delivered reproduction package includes the raw SPARC tables and theory-compression input. A complete 30-page atlas contains all 175 galaxies, with no poor fit removed.

\begingroup\small\setlength{\parskip}{0pt}
\begin{thebibliography}{9}
\bibitem{sparc} F. Lelli, S. S. McGaugh and J. M. Schombert, \emph{SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves}, AJ 152, 157 (2016), \href{https://arxiv.org/abs/1606.09251}{arXiv:1606.09251}. Public tables: \url{https://astroweb.case.edu/SPARC/}.
\bibitem{extension} User-supplied \emph{Coupled two-field--Friedmann solutions and rotation-curve envelope approximations}, 29 September 2026. Source and diagnostics are included in the reproduction inputs.
\bibitem{project} Project manuscripts \emph{Cubic Scalar--GR Intersection} and Supplement; consolidated study version 4; archived theory compression (28 September 2026), included as \texttt{inputs/compression.json}.
\end{thebibliography}
\endgroup
\end{document}
'''
tex=tex.replace('@@ALLTABLE@@',table_summary('all')).replace('@@OUTTABLE@@',table_summary('inner')).replace('@@STATS175@@',statrows('all175')).replace('@@STATS131@@',statrows('quality131'))
(OUT/'SPARC_full_profile_report.tex').write_text(tex,encoding='utf-8')
for j in range(3):
 p=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory',str(TMP),'SPARC_full_profile_report.tex'],cwd=OUT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (TMP/f'pass{j+1}.txt').write_text(p.stdout,encoding='utf-8')
 if p.returncode:print(p.stdout[-6000:]);raise SystemExit(p.returncode)
shutil.copy2(TMP/'SPARC_full_profile_report.pdf',OUT/'SPARC_full_profile_report.pdf')
print('Produced report and atlas:',OUT)
print('\n'.join(l for l in (TMP/'SPARC_full_profile_report.log').read_text(errors='replace').splitlines() if any(x in l for x in ['Warning','Overfull','Underfull','Output written'])))
