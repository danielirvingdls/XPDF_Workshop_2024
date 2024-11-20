'pdf.inc, by Phil Chater
'Copy this file to your TOPAS folder and add the following line to your local.inc 
'#include "C:\TOPAS-6\pdf.inc" 

'===INSTRUMENTAL/PROCESSING FUNCTIONS===
macro dQ_damping(dQ,dQv)
{
	'Models dampening of the PDF intensity as a function of r due to a constant width reciprocal space peak shape
	'Ref: Acta Cryst. (1992). A48, 336-346
	#m_argu dQ 'Instrumental FWHM of S(Q) data
	If_Prm_Eqn_Rpt(dQ, dQv, min 0.001 max 0.2, del 0.0001)
	scale_phase_X = Exp(-0.5  X^2 (CeV(dQ,dQv)/2.35482)^2);
}

macro dQ_lor_damping(dQ,dQv,lor,lorv)
{
	'A modification of dQ_dampening to correct for a Lorentzian component to the reciprocal space peak shape
	#m_argu dQ 'Instrumental FWHM of S(Q) data
	#m_argu lor 'Lorentzian contribute to peak shape (between 0 and 1)
	If_Prm_Eqn_Rpt(dQ, dQv, min 0.001 max 0.2, del 0.0001)
	If_Prm_Eqn_Rpt(lor, lorv, min 0 max 1)
	scale_phase_X = (1-CeV(lor,lorv)) Exp(-0.5  X^2 (CeV(dQ,dQv)/2.35482)^2) + CeV(lor,lorv) Exp(-0.5 CeV(dQ,dQv) X);
}

macro convolute_Qmax_Sinc(Qmax,Qmaxv)
{
	'Convolution to include Sinc function due to termeination of S(Q) at finite Qmax during Fourier transform
	'Use in combination with start_X 1 to avoid very low-r data 
	'Ref: J. S. Chung and M. F. Thorpe, Phys. Rev. B, 55(3), 1997, 1545
	#m_argu Qmax
	If_Prm_Eqn_Rpt(Qmax, Qmaxv, min 1 max 100, del 0.01)
	local #m_unique conv_max = (5 CeV(Qmax,Qmaxv) - Mod(5 CeV(Qmax,Qmaxv),1))/5 2 Pi / CeV(Qmax,Qmaxv);
	pdf_convolute = If(Abs(X) > Yobs_dx_at(Xo),(Sin(CeV(Qmax,Qmaxv) X)/(X)),CeV(Qmax,Qmaxv));
		min_X = Min(-Xo,-conv_max) ;
		max_X = Max( Xo, conv_max) ;
}
macro convolute_Qmax_Sinc(Qmax) { convolute_Qmax_Sinc( ,Qmax) }

macro convolute_alpha(alpha,alphav)
{
	'Convolution to account for a Q dependent broadening term in the S(Q) data
   #m_argu alpha
   If_Prm_Eqn_Rpt(alpha, alphav, min 0.00001 max 0.1, del 0.001)
   local #m_unique fwhm = Max(CeV(alpha,alphav) Xo,0.00001);
   local #m_unique conv_max = Ln(0.001) / fwhm;
   pdf_convolute = Gauss(0, fwhm);
   	min_X = Min(-Xo,-conv_max) ;
   	max_X = Max( Xo, conv_max) ;
   	convolute_X_recal = If(Xo,1,1) 1;
}
macro convolute_alpha(alpha) { convolute_alpha( ,alpha) }

macro convolute_Lorch(Qmax,Qmaxv)
{
	'Convolution to account for the use of a Lorch function prior to Fourier transform of the S(Q)
   #m_argu Qmax
   If_Prm_Eqn_Rpt(Qmax, Qmaxv, min 5 max 60, del 0.1)
   local #m_unique fwhm = (Pi / CeV(Qmax,Qmaxv))/Ln(2);
   pdf_convolute = Gauss(0, fwhm);
   	min_X = Ln(0.001) / fwhm;
   	max_X =-Ln(0.001) / fwhm;
   	convolute_X_recal = If(Xo,1,1) 1;
}
macro convolute_Lorch(Qmax) { convolute_Lorch( ,Qmax) }

macro convolute_SoperLorch(d_zero,d_zerov,d_one,d_onev)
{
	'Covolution to account for the use of the Soper-Lorch function prior to Fourier transform of the S(Q)
	'd_zero is the "Width of broadening in r space", as used in GudrunN
	'Ref: A. K. Soper and E. R. Barney, J. Appl. Cryst. (2012). 45, 1314–1317
   #m_argu d_zero
   #m_argu d_one
   If_Prm_Eqn_Rpt(d_zero, d_zerov, min 0.001 max 1.0, del 0.01)
   If_Prm_Eqn_Rpt(d_one, d_onev, min 0 max 0.1, del 0.01)
   local #m_unique fwhm = 2.3548 (CeV(d_zero,d_zerov)/Ln(2) Ln(2 Cosh( (CeV(d_one,d_onev)/CeV(d_zero,d_zerov)) Ln(2) Xo)));
   pdf_convolute = Gauss(0, fwhm);
   	min_X = Ln(0.001) / fwhm;
   	max_X =-Ln(0.001) / fwhm;
   	convolute_X_recal = If(Xo,1,1) 1;
}
macro convolute_SoperLorch(d_zero,d_zerov) {	convolute_SoperLorch(d_zero,d_zerov, ,0) }
macro convolute_SoperLorch(d_zero) {	convolute_SoperLorch( ,d_zero, ,0) }

'===BEQ PEAK SHAPE FUNCTIONS===
macro beq_rcut(rcut,beqlo,beqlov,beqhi,beqhiv)
{
	#m_argu beqlo
	#m_argu beqhi
	#m_argu rcut
	If_Prm_Eqn_Rpt(beqlo, beqlov, min 0.001 max 10, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqhi, beqhiv, min 0.001 max 10, del = 0.05 Val;)
	beq = IF X < rcut THEN
				CeV(beqlo, beqlov)
			ELSE
				CeV(beqhi, beqhiv)
			ENDIF;
}

macro beq_rcut_rlo_spherical(rcut,rcutv,beqcut,beqcutv,rlo,rlov,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
{
	#m_argu rcut
	#m_argu beqcut
	#m_argu rlo
	#m_argu beqlo
	#m_argu beqhi
	#m_argu radius
	If_Prm_Eqn_Rpt(rcut, rcutv, min 0.1 max 10, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqcut, beqcutv, min 0.001 max 10, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(rlo, rlov, min 0.1 max 100, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqlo, beqlov, min 0.001 max 10, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqhi, beqhiv, min 0.001 max 100, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(radius, radiusv, min 1 max 1000, del = 0.05 Val;)
	beq = 
		IF X > 0.1 THEN
			IF X > (2 CeV(radius, radiusv)) + CeV(rlo, rlov) THEN
				CeV(beqhi, beqhiv)
			ELSE
				IF X > CeV(rlo, rlov) THEN
					CeV(beqlo, beqlov) + (CeV(beqhi, beqhiv)-CeV(beqlo, beqlov)) 
						(1-(Pi (X-CeV(rlo, rlov))^2 ((0.25 ((X-CeV(rlo, rlov))/CeV(radius, radiusv))^3)-(3 (X-CeV(rlo, rlov))/CeV(radius, radiusv))+4)/(4 Pi (X-CeV(rlo, rlov))^2)))
				ELSE
					IF X > CeV(rcut, rcutv) THEN
						CeV(beqlo, beqlov)
					ELSE
						CeV(beqcut, beqcutv)
					ENDIF
				ENDIF
			ENDIF
		ELSE
			0.01
		ENDIF;
}
macro beq_rcut_spherical(rcut,rcutv,beqcut,beqcutv,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
{
	beq_rcut_rlo_spherical(rcut,rcutv,beqcut,beqcutv,!rlo,0,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
}
macro beq_spherical(beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
{
	beq_rcut_rlo_spherical(!rcut,0,!beqcut,0,!rlo,0,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
}
macro beq_rlo_spherical(rlo,rlov,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
{
	beq_rcut_rlo_spherical(!rcut,0,!beqcut,0,rlo,rlov,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)	
}

macro pairs_spherical(rlo,rlov,beqlo,beqlov,beqhi,beqhiv,radius,radiusv)
{
	#m_argu rlo
	#m_argu beqlo
	#m_argu beqhi
	#m_argu radius
	If_Prm_Eqn_Rpt(rlo, rlov, min 0.1 max 100, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqlo, beqlov, min 0.001 max 10, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(beqhi, beqhiv, min 0.001 max 100, del = 0.05 Val;)
	If_Prm_Eqn_Rpt(radius, radiusv, min 1 max 1000, del = 0.05 Val;) 
	pdf_gauss_fwhm = 
		IF X > (2 CeV(radius, radiusv)) + CeV(rlo, rlov) THEN
			CeV(beqhi, beqhiv)
		ELSE
			IF X > CeV(rlo, rlov) THEN
				CeV(beqlo, beqlov) + (CeV(beqhi, beqhiv)-CeV(beqlo, beqlov)) 
					(1-(Pi (X-CeV(rlo, rlov))^2 ((0.25 ((X-CeV(rlo, rlov))/CeV(radius, radiusv))^3)-(3 (X-CeV(rlo, rlov))/CeV(radius, radiusv))+4)/(4 Pi (X-CeV(rlo, rlov))^2)))
			ELSE
				CeV(beqlo, beqlov)
			ENDIF
		ENDIF;
}
macro pairs_spherical(rlo,rlov,beqlo,beqhi,radius,radiusv)
{
	pairs_spherical(rlo,rlov, ,beqlo, ,beqhi,radius,radiusv)
}

macro beq_r_r2(beqc,beqv,r,rv,r2,r2v)
{
	#m_argu beqc
	#m_argu r
	#m_argu r2
	If_Prm_Eqn_Rpt(beqc, beqv, min 0.001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(r, rv, min 0.001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(r2, r2v, min 0.00001 max 10, del = 0.01 Val;)
	beq = CeV(beqc, beqv) + CeV(r, rv) X + CeV(r2, r2v) X^2;
}

macro beq_empirical(v1c,v1v,v2c,v2v,v3c,v3v,v4c,v4v)
{
	#m_argu v1c
	#m_argu v2c
	#m_argu v3c
	#m_argu v4c
	If_Prm_Eqn_Rpt(v1c, v1v, min 0.000001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(v2c, v2v, min 0.000001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(v3c, v3v, min -100.0   max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(v4c, v4v, min -100.0   max 10, del = 0.01 Val;)
	beq = CeV(v1c, v1v) + CeV(v2c, v2v) X + CeV(v3c, v3v) / X + CeV(v4c, v4v) / X^2;
}

macro beq_PDFfit2(uiso,uisov,rcut,rcutv,sratio,sratiov,delta1,delta1v,delta2,delta2v,qbroad,qbroadv)
{
	#m_argu uiso
	#m_argu rcut
	#m_argu sratio
	#m_argu delta1
	#m_argu delta2
	#m_argu qbroad
	If_Prm_Eqn_Rpt(uiso, uisov,  min 0.000001 max 1, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(rcut,  rcutv,   min 0.000001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(sratio,sratiov, min 0.000001 max 10, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(delta1,delta1v, min 0.000001 max 1, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(delta2,delta2v, min 0.000001 max 1, del = 0.01 Val;)
	If_Prm_Eqn_Rpt(qbroad,qbroadv, min 0.000001 max 1, del = 0.01 Val;)
	beq = If(X < CeV(rcut, rcutv), 0.5 CeV(sratio, sratiov), 1) CeV(uiso, uisov) 8 Pi^2 (1 - Min(((CeV(delta1, delta1v)/X) + (CeV(delta2, delta2v)/X^2) - (CeV(qbroad, qbroadv)^2 X^2)),1));	
}

'===SAMPLE FUNCTIONS===
macro spherical_damping(r,rv)
{
	#m_argu r
	If_Prm_Eqn_Rpt(r, rv, min 1 max 1000, del = 0.01 Val;)
	scale_phase_X = 
	IF X > 2 CeV(r,rv) THEN
		0
	ELSE
		If(X>0.01,(Pi X^2 ((0.25 (X/CeV(r,rv))^3)-(3 X/CeV(r,rv))+4)) / (4 Pi X^2 1),1)
	ENDIF;
}

macro local_sphere(r,rv)
{
	spherical_damping(r,rv)
}
macro local_sphere(rv) { local_sphere( ,rv) }

macro longrange_sphere(r,rv)
{
	#m_argu r
	If_Prm_Eqn_Rpt(r, rv, min 1 max 1000, del = 0.01 Val;)
	scale_phase_X = 
	IF X > 2 CeV(r,rv) THEN
		1
	ELSE
		If(X>0.01,1-((Pi X^2 ((0.25 (X/CeV(r,rv))^3)-(3 X/CeV(r,rv))+4)) / (4 Pi X^2 1)),0)
	ENDIF;
}
macro longrange_sphere(rv) { longrange_sphere( ,rv) }

'===PDF SCALING FUNCTIONS===
macro pdf_D_dash_of_r_bkg(num_density,num_densityv)
{
	#m_argu num_density
	If_Prm_Eqn_Rpt(num_density, num_densityv, min 0 max 1, del 0.001)
	fit_obj = -4 Pi X CeV(num_density,num_densityv);
}

macro G_dash_of_r_scaling(num_density,num_densityv)
'Can be used either in xdd with the total sample number density or in a phase with the phase number density
{
	'Rescaling to G'(r) according to Keen, J. Appl. Cryst. (2001). 34, 172-177
	'num_density is the number density in atoms per cubic Angstrom for the phase (within str) or sample (within xdd) 
	#m_argu num_density
	If_Prm_Eqn_Rpt(num_density, num_densityv, min 0 max 1, del 0.001)
	scale_phase_X = If(X>0.01,1/(4 Pi CeV(num_density,num_densityv) X),0);
   fit_obj = If(X>0.01,1,0);
}

macro G_of_r_scaling(num_density,num_densityv,weighted_av,weighted_avv)
{
	'Rescaling to G(r) according to Keen, J. Appl. Cryst. (2001). 34, 172-177
	'num_density is the number density in atoms per cubic Angstrom for the phase (within str) or sample (within xdd)
	'weighted_av is the weighted average structure factor, <F>^2
	#m_argu num_density
	#m_argu weighted_av
	If_Prm_Eqn_Rpt(num_density, num_densityv, min 0 max 1, del 0.001)
	If_Prm_Eqn_Rpt(weighted_av, weighted_avv, min 0 max 2, del 0.001)
	scale_phase_X = If(X>0.01,CeV(weighted_av, weighted_avv) 1/(4 Pi CeV(num_density,num_densityv) X),0);
   fit_obj = If(X>0.01,0,-CeV(weighted_av, weighted_avv));
}

'===BRAGG PEAKSHAPE FUNCTIONS===
macro pkshape_dQ_alpha(dQ,dQv,alpha,alphav,lor,lorv)
{	
	'Bragg peak shape of the form (dQ^2 + alpha^2 Q^2)^0.5 for 2theta data
	'	Q = (2 Pi/D_spacing)
	'	deltaQ/Q = (dQ^2 + alpha^2 Q^2)^0.5 / Q
	'	delta2Th = 2 deltaQ/Q Tan(Th)
	#m_argu dQ 		'Fixed peak width, in Q
	#m_argu alpha 	'Q dependent peak width term
	#m_argu lor		'Lorch fraction
	If_Prm_Eqn_Rpt(dQ, dQv, min 0.001 max 0.1, del 0.0001)
	If_Prm_Eqn_Rpt(alpha, alphav, min 0.001 max 0.05, del 0.001)
	If_Prm_Eqn_Rpt(lor, lorv, min 0 max 1)
	peak_type pv
		pv_lor = CeV(lor,lorv);		
		pv_fwhm = 2 (CeV(dQ,dQv)^2 + CeV(alpha,alphav)^2 (2 Pi/D_spacing)^2)^0.5 / (2 Pi/D_spacing)  Tan(Th) Rad;
}

macro pkshape_dQ(dQ,dQv,lor,lorv)
{
	pkshape_dQ_alpha(dQ,dQv,!alpha,0,lor,lorv)
}

'===MISCELLANEOUS FUNCTIONS===
macro Out_PDF(file)
{
   xdd_out file load out_record out_fmt out_eqn
   {
      " %11.6f  " = X; 
   	" %11.6f  " = Yobs;
      " %11.6f  " = Ycalc;
      " %11.6f\n" = Yobs-Ycalc;
   }
}

macro pdf_peak(xo, xov, a, av, fwhm, fwhmv)
'For use when peak fitting to PDF without a phase present
{
   #m_argu xo
   #m_argu a
   #m_argu fwhm
   If_Prm_Eqn_Rpt(xo, xov, )
   If_Prm_Eqn_Rpt(a, av, )
   If_Prm_Eqn_Rpt(fwhm, fwhmv,  min 0.001 del 0.001)   
   fit_obj = CeV(a,av) Gauss(CeV(xo,xov), CeV(fwhm,fwhmv));
   fit_obj_phase 99
}

macro pdf_CN(xo, xov, CN, CNv, beq, beqv)
'Calculates coordination numbers for monatomic phases
'Peak width is the equivalent of two atoms of the same beq value
{
   #m_argu CN
   #m_argu xo
   #m_argu beq
   If_Prm_Eqn_Rpt(CN, CNv,  min 0.00001)
   If_Prm_Eqn_Rpt(xo, xov, )
   If_Prm_Eqn_Rpt(beq, beqv,  min 0.00001 del 0.001)   
   fit_obj = (CeV(CN,CNv) / CeV(xo,xov)) Gauss(CeV(xo,xov), 2.3548 (2 CeV(beq,beqv))^0.5); ''2 (2 Ln(2))^0.5 = 2.3548
}

macro pdf_gauss_fwhm_beqs(beq1, beq2)
{
   pdf_gauss_fwhm = 2.3548 ((beq1 + beq2) / (8 Pi^2))^0.5;
}

macro fwhm_beqs(beq1, beq2)
{
   2.3548*((beq1/(8 Pi^2))^2 + (beq2(8 Pi^2))^2)^0.5
}

macro erf_fwhm(r,fwhm,fwhmv)
'Defines a crossover (error) function from 0 to 1 centered at r=0, with an equivalent Gaussian width fhwm
'NOTE: Erf_Approx is the error function for a sigma of ln(2). 1.632237 is ln(2)*2(2ln(2))^0.5 to convert this to fwhm.
{
	#m_argu r
	#m_argu fwhm
	If_Prm_Eqn_Rpt(fwhm, fwhmv,  min 0.01 max 100 del 0.001)
	(Erf_Approx(r 1.632237/CeV(fwhm,fwhmv))+1)/2
}
macro erf_fwhm(r,fwhm) { erf_fwhm(r, ,fwhm) }

'=======CHANGE LOG=======
'===Updated 02/06/2017===
'	beq_spherical_rmin_rcut added
'	pkshape_dQ_alpha updated
'	dQ_lor_damping added
'===Updated 08/06/2017===
'  beq_r_min updated to beq_rcut for consistency
'  beq_spherical functions updated to refer to same master function
'===Updated 13/06/2017===
'  convolution functions fixed using convolute_X_recal = If(Xo,1,1) 1;
'===Updated 14/06/2017===
'  beq_rcut_rlo_spherical corrected for rlo != 0 and made more robust
'===Updated 17/07/2019===
'  single parameter macros for local_sphere and longrange_sphere added
'===Updated 18/07/2019===
'  beq_rcut_rlo_spherical updated to fix problem with beq_rcut_spherical
'  min value for alpha changed to 0.0001 as smaller values causes TOPAS to run very slowly
'===Updated 09/06/2022===
'  convolute_alpha updated to use Max(CeV(alpha,alphav) Xo,0.00001) for use with TOPAS v7
'  beq_rcut_rlo_spherical updated to use IF X > 0.1 THEN for use with TOPAS v7