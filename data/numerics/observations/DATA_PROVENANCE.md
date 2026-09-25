# Observational data and comparison protocol

Retrieved 2026-09-23 from the SPARC team site. This note records provenance and conventions for the delivered observational comparison.

## Primary sources

- [SPARC team data portal](https://astroweb.case.edu/SPARC/).
- [SPARC metadata, SPARC_Lelli2016c.mrt](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt).
- [SPARC radial mass models, MassModels_Lelli2016c.mrt](https://astroweb.case.edu/SPARC/MassModels_Lelli2016c.mrt).
- Lelli, McGaugh and Schombert, *Astronomical Journal* **152**, 157 (2016), [doi:10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157); [author PDF](https://astroweb.case.edu/ssm/papers/AJv152n157.pdf). Sections 3.2.2 and 3.3 give the selection, velocity error model and signed component convention.
- [CDS journal catalogue documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/AJ/152/157?format=html&tex=true), including units and flags.
- McGaugh, Lelli and Schombert, *Physical Review Letters* **117**, 201101 (2016), [doi:10.1103/PhysRevLett.117.201101](https://doi.org/10.1103/PhysRevLett.117.201101); [author PDF](https://arxiv.org/pdf/1609.05917). The "Stellar Mass-to-Light Ratios" discussion explicitly fixes disk and bulge values at 0.5 and 0.7 solar units, respectively.
- Lelli, McGaugh, Schombert and Pawlowski, *Astrophysical Journal* **836**, 152 (2017), [author PDF](https://arxiv.org/pdf/1610.08981), Table 1 and Section 4.1: the same 0.5/0.7 convention and its dependence on population-synthesis assumptions.
- Navarro, Frenk and White, *Astrophysical Journal* **490**, 493-508 (1997), [doi:10.1086/304888](https://doi.org/10.1086/304888); [author PDF](https://arxiv.org/pdf/astro-ph/9611107). Equations (1) and (3) specify the density and circular-speed profiles.

## Reproducible data facts

Raw files are in `data/` next to this note. SHA256:

- metadata: `5AA0501F6B0D881FA579030E315E7B5B6EF561A5BD3A07472F9929C7E5728243`
- radial table: `9108994B12CC401B94A1768BECA61C53EC354779385C9C9CC571049F3043244C`

A direct parse gives 175 galaxies and 3391 radial rows. The published galaxy-level selection Q<3 and inclination>=30 degrees retains 153 galaxies and **3168 rows in these downloaded files**. This count must be reported from the actual download, not replaced by the 3149 pre-cut or 2693 precision-selected rows quoted in an earlier RAR-paper version. Do not add a velocity-error cut unless separately declared. The smallest selected galaxy has four rows. Every radial velocity error is positive (minimum 0.2 km/s).

Whitespace parsing is appropriate. The metadata file has 19 fields per data row; its printed byte ranges do not match its current row widths. `results/selected_metadata.json` records the selected galaxy metadata. The original tables retain distance methods, uncertainty columns and kinematic reference codes; these raw files are authoritative.

## Column and sign conventions

Radial-table fields: galaxy name; adopted distance [Mpc]; radius [kpc]; observed speed and its error [km/s]; gas, stellar disk and bulge velocity components [km/s]; disk and bulge surface brightness [solar luminosity/pc^2]. Stellar components are tabulated at M/L=1. Gas already includes the 1.33 helium correction. Do not apply that correction again.

The baryonic radial-force contribution is

`B(R) = Vgas*abs(Vgas) + Upsilon_d*Vdisk*abs(Vdisk) + Upsilon_b*Vbul*abs(Vbul)`.

Negative component values encode an outward radial force from a centrally depressed mass distribution. The raw release contains 361 negative gas-component rows. With Upsilon_d=0.5 and Upsilon_b=0.7, the selected sample contains two negative *total* B values, both UGC01281: R=0.08 kpc gives -33.0975 (km/s)^2; R=0.23 kpc gives -27.92125. A baryons-only circular speed is undefined there. Never conceal this by clipping B to zero or taking an absolute value. Mark these baseline predictions undefined; if using circular-speed chi-square, report its 3166-row domain and use a common valid-row domain when comparing aggregate baseline differences. Augmented models can retain these rows if their total squared speed is positive.

## Uncertainties and sample choice

The catalogue speed errors combine fit uncertainty and the approaching/receding asymmetry estimate: e_V^2=e_fit^2+[(Vapp-Vrec)/4]^2. They omit inclination uncertainty. Distance, inclination and M/L produce shared systematic shifts, and the release provides no radial covariance matrix. Catalogue-conditional fits with those nuisance quantities fixed are reproducible descriptive comparisons, not a complete likelihood or model evidence calculation. Q=3 flags substantial asymmetries/noncircular motions; Q=2 allows smaller irregularities. Very slow gas systems require additional attention to pressure support; the release does not supply a uniform pressure-correction likelihood.

The declared descriptive protocol fits every galaxy passing Q<3 and i>=30, and displays six evenly spaced ranks of per-galaxy median observed speed, ties by name. The display selection is fixed before considering fit performance. Full fit tables, boundary flags and residual diagnostics are supplied; all optimizers succeeded and unfavorable fits are retained. The included M/L sensitivity changes disk and bulge ratios together. Distance/inclination marginalization and a complete uncertainty model remain necessary for inferential claims; formal fitting errors would not substitute for them.

## NFW comparator and interpretation

Use rho_h(r)=rho_s/[x(1+x)^2], x=r/r_s, with rho_s,r_s>0. Direct integration gives

`Vh^2(r)=4*pi*G*rho_s*r_s^2*[log(1+x)-x/(1+x)]/x`.

With v_h^2=4*pi*G*rho_s*r_s^2, the two fitted positive scales can be (r_s,v_h). Evaluate the bracket stably at small x; its quotient is x/2-2*x^2/3+3*x^3/4+O(x^4). The unrestricted two-scale template imposes no cosmological concentration relation and should be called a GR+NFW halo comparator, not the unique prediction of GR or of LambdaCDM. Report large-r_s boundary cases: in that limit Vh^2 approaches a constant times r and the separate scales are weakly identified.

Baryons-only, baryons+NFW and baryons+the manuscript's scalar-density template all use the weak-field gravitational source relation. GR alone does not determine a galaxy's matter profile. An additive standalone scalar-core template ignores baryonic backreaction and does not establish a universal microscopic parameter fit. Internal material stress must not be inserted as an independent force on gas without a derived coupling or as a second gravitational source already counted by the scalar density.

## Which existing figures can become observational?

The circular-speed comparison can use SPARC observations and its baryonic decomposition directly. A separate controlled calculation can compare the model's internal stress with a precisely defined baseline, but that calculation remains theoretical. The signed-angle KL-response curve cannot be replaced by SPARC: a single observed rotation curve supplies neither controlled changes of the internal angle nor repeated distributions of the specified tracer readout. The entropy/ergotropy figures illustrate exact formulas and have no direct corresponding measurement in this dataset.
