%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SISTER: Starshade Imaging Simulation Toolkit for Exoplanet Reconnaissance
%
% Copyright (c) <2019>, California Institute of Technology ("Caltech"). U.S. Government sponsorship acknowledged.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% • Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% • Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
% • Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ fluxgridst, fluxgrid ] = glint_distrib4imaging( opt, sz_scn_cnv )
% compute distribution of solar glint on a starshade
% Uses the lab data processed and calibrated from scatterometerScaleSingle.m

% fluxgrid(st) are square grids with the same pixel scale than the convolved scene.
% fluxgrid(st) is dimensionless, it is a flux ratio. Later on, in convolve_with_one_wavelength.m it is multipilied by the
% Sun's apparent flux as well as convolved with the telescope optical response (at the true distance of the starshade to 
% the telescope, so for large telescopes the PSF is 'blurred' wrt the far distance response).  
% Starshade will be centered
% Author: Stuart B. Shaklan, Stuart.B.Shaklan@jpl.nasa.gov
% Modified: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com
% 101519: Stuart B. Shaklan 
% Modifications from original glint_distrib4imaging.
%  We now read in ZS, ZP, ZSst, and ZPst (stealth) rather than generating
%  them from the lab data.  That function is done outside of SISTER.
% ZS, ZP are 1801 x 701,  First dimension is -90 to 90 deg in 0.1 deg
% steps, Second dimension is 15-85 deg in 0.1 deg steps
% Units of Z files are fractional flux per m of edge, per m^2 of aperture, m^2 of
% distance (divide by distance squared to get fractional flux at telescope)
% So, to get Jy per meter squared at the telescope, compute
%    solar Jy/m2 incident on the edge * ZS (or ZP) / (range to starshade)^2 * lambda/635 nm (lab wavelength) * length of edge.
% Length of edge is computed in fluxgrid_f.
% 013020: Stuart B. Shaklan
% Added the case of displaced starshades

%  NOTE:
%  Diffraction term goes as lambda, not 1/lambda. Sommerfeld is E~ 1
%  sqrt(k*z) = sqrt(lambda/Z) so power goes as lambda/Z.
%  Reflection term is geometrical, so independent of lambda.
%  HERE WE ASSUME WE ARE DOMINATED BY DIFFRACTION. SO ASSUME LAMBDA
%  DEPENDENCE.

load( [ opt.scn_dr '/solar_glint/' opt.slr_glnt.zs_fl] )  % EACH ARRAY REPRESENTS THE FRACTION OF S OR P LIGHT THAT IS SCATTERED.
ZS = Z_out/2; %SBS 102719:  ZS is fraction of S light, which is half of the sunlight
load( [ opt.scn_dr '/solar_glint/' opt.slr_glnt.zp_fl] )
ZP = Z_out/2;  %SBS 102719:  ZS is fraction of S light, which is half of the sunlight

% New petal configuration. For instance, with shredded edges (stealth data)
load( [ opt.scn_dr '/solar_glint/' opt.slr_glnt.zs2_fl] )
ZSst = Z_out/2 ; %SBS 102719:  ZS is fraction of S light, which is half of the sunlight
load( [ opt.scn_dr '/solar_glint/' opt.slr_glnt.zp2_fl] )
ZPst = Z_out/2; %SBS 102719:  ZS is fraction of S light, which is half of the sunlight

% Rescale for telescopeRange (increases as range) and wavelength
% increases with R
% (goes as lam because diffraction is 1-D from the cylinder)
% Correcting for wavelength and distance between the Starshade and the telescope
ZS = ZS * ( 1 / opt.dst_strshd_tlscp_m )^2 * ( opt.lmbd_tmp_nm / 635 ) ;
ZP = ZP * ( 1 / opt.dst_strshd_tlscp_m )^2 * ( opt.lmbd_tmp_nm / 635 ) ;
ZSst = ZSst * ( 1 / opt.dst_strshd_tlscp_m )^2 * ( opt.lmbd_tmp_nm / 635 ) ;
ZPst = ZPst * ( 1 / opt.dst_strshd_tlscp_m )^2 * ( opt.lmbd_tmp_nm / 635 ) ;
% Rotate starshade
% Loading the petals' position. If there are some perturbations, load the perturbed locus (for WFIRST non-ideal starshade, the solar glint flux differences are very small <1e-3 relative between non-ideal and ideal petals)
  if ( opt.lcs.do )
  load( [ opt.path_occulter opt.locus.perturbed_filename ] )
  else
  load( [ opt.path_occulter opt.starshade.nominal_filename ] )
  end

% Petal locus
  if ~exist( 'xVals', 'var' )
  % This is very fast ~0.01 sec all
  load( sprintf( '%s/locus/in/%s.mat', opt.scene_dir, opt.starshade.nominal_filename ) )
  vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, opt.n_ptl, {}) ;
  tmpxVals = [];
  tmpyVals = [];
  tmpzVals = [];
    for j = 1 : opt.n_ptl
    tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)];
    tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)];
    tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)];
    end
  xVals = [tmpxVals tmpxVals(1)];
  yVals = [tmpyVals tmpyVals(1)];
  zVals = [tmpzVals tmpzVals(1)];
  % If they did not exist, store them
  fl_nm_xyzvals = sprintf( '%s/locus/in/%s_xyzVals.mat', opt.scene_dir, opt.starshade.nominal_filename ) ;
    if exist( fl_nm_xyzvals, 'file' ) ~= 2
      if sum( zVals )
      save( fl_nm_xyzvals, 'xVals', 'yVals', 'zVals' ) ;
      else
      save( fl_nm_xyzvals, 'xVals', 'yVals' ) ;
      end
    end
  end

% telback below uses ones(1,length(diffx)
  if size( xVals, 1 ) ~= 1
  xVals = xVals' ;
  end
% This should always be consistent with xVals
  if size( yVals, 1 ) ~= 1
  yVals = yVals' ;
  end
% If the starshade is non-spinning, we only consider 1 position (alpha may still rotate the petals of the starshade by some angle). For the case of a spinning starshade, we consider all necessary rotations and average out the result. PS: about 0.1 s per loop element in a MacBook Pro 2017, 3.1 GHz Intel Core i7. Below, one petal is fully rotated if the pupil is ideal. 
  n_rot = 0 ;
    if strcmp( opt.starshade.mode, 'spinning' )
    dlt_ang_dg = 0.5 ;  % This step size keeps the overall precision better than 1% (since the solar glint may be subtracted in a post-processing step, the residual error on the solar glint fluctuations -Poisson- is even smaller)
    % Bulding the list of all angles of a petal, excluding 0, which is already computed, and equivalent rotations that would fully rotate a petal, i.e. equivalent to 0.
    ssrot = dlt_ang_dg : dlt_ang_dg : 360 / opt.n_ptl / 2  ; % half of the angles of a petal, excluding 0.
    % avoiding repeating the same intial/final position if it happens:
      if ssrot( end ) == 360 / opt.n_ptl / 2
      ssrot = [ -ssrot ssrot( 1 : end - 1 ) ] ; % +/- angles for a petal
      else
      ssrot = [ -ssrot ssrot ] ; % +/- angles for a petal
      end
    n_rot = numel( ssrot ) ;
    end

% For the zero angle first
[ fluxgrid_0 fluxgridst_0 ] = fluxgrid_f( xVals, yVals, opt, ZS, ZP, ZSst, ZPst, sz_scn_cnv, 0 ) ;
% Half of the non-zero angles
  for ii = 1 : n_rot
  [ fluxgrid_tmp fluxgridst_tmp ] = fluxgrid_f( xVals, yVals, opt, ZS, ZP, ZSst, ZPst, sz_scn_cnv, ssrot( ii ) ) ;
    if ii == 1
    fluxgrid = fluxgrid_tmp * 0 ;
    fluxgridst = fluxgrid ;
    end
  fluxgrid = fluxgrid + fluxgrid_tmp ;
  fluxgridst = fluxgridst + fluxgridst_tmp ;
  end % ii
% Averaging the accumulated result: 0 +/- rest of the non-zero angles
  if ( n_rot )
  fluxgrid = ( fluxgrid_0 + fluxgrid ) / ( 1 + n_rot ) ; 
  fluxgridst = ( fluxgridst_0 + fluxgridst ) / ( 1 + n_rot ) ; 
  else
  fluxgrid = fluxgrid_0 ;
  fluxgridst = fluxgridst_0 ;
  end

% Sample plot (run scene_3 and stop here):
 if ( 0 )
 clf
 opt_plt.mnmx = [-1e-24,5e-24];opt_plt.zm=4;plt_mnmx( fluxgrid,opt_plt ) ;
 xlabel( 'PIX (3 mas)', 'FontSize', 14 )
 ylabel( 'PIX (3 mas)', 'FontSize', 14 )
 hold on
 title( sprintf( 'New glint%sdistrib4imaging (%i nm)', '\_', opt.lmbd_tmp_nm ), 'FontSize', 14 )
 end

%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%
% Subfunction that derives the flux from the solar glint for a given orientation of the petals
function [ fluxgrid fluxgridst ] = fluxgrid_f( xVals, yVals, opt, ZS, ZP, ZSst, ZPst, sz_scn_cnv, ssrot )
stealth_exclusion_angle = 2 ; % degrees. +/- this theta angle of sun will be replaced with stealth edges.
xVals_rtd =   xVals*cosd( ssrot ) + yVals*sind( ssrot ) ;
yVals_rtd = - xVals*sind( ssrot ) + yVals*cosd( ssrot ) ;
% star-starshade-sun angle. (Stuart Shaklan 11/17/18:  Phi = 90 has the sun off to the side of the starshade.  Phi = 0 has the sun behind the starshade.  Missions will generally restrict phi angles between 45 and 85 degrees).
phi = round( opt.slr_glnt.phi_dg * 10, 1 ) / 10 ;  % round to nearest 0.1 deg, anticipating that thetaS and phiS arrays are in 0.1 deg steps.
diffx = diff(xVals_rtd) ;
diffy = diff(yVals_rtd) ;
diffnorm = sqrt(diffx.^2 + diffy.^2) ;
diffx = diffx./diffnorm ;
diffy = diffy./diffnorm ;
% vector from telescope to angle phi backward from starshade.
% once we have that, we can see if within 0.25 deg of sun.
% telback = [-diffy*sind(phi); -diffx*sind(phi); cosd(phi)*ones(1,length(diffx))];  % this is in the plane of the starshade. Not used.
% Now go backwards to the phi of the star. relative x and
% y components stay the same, but need to renormalize
% because of non-zero z component.  Assumes telescope is
% in negative z direction.  Sun is in positive z plane.
alpha = opt.slr_glnt.alph_dg ; % rotation of sun around the plane.  zero is along x axis.
%k0 = [sind(phi)*cosd(alpha); sind(phi)*sind(alpha); cosd(phi)] * ones(1,length(diffx));  % sun vector.  Not used.
theta = alpha+90+180/pi*atan2(diffx,diffy);
theta = round(theta,1);  % round to nearest 0.1 deg. Anticipating that thetaS and phiS are arrays in 0.1 deg increments.
theta(theta>270)=theta(theta>270)-360;  % always want theta values between -90 and 270 for leading/trailing discrimination
theta(theta<-90)=theta(theta<-90)+360;
% SBSSBSSBS   013020
% for starshades that are shifted (e.g. formation flying errors), our check
% on minimum radius for glint (edges are blocked at radii smaller than r_m)
% needs to account for the position of the starshade.
xmean = mean(xVals_rtd);
ymean = mean(yVals_rtd);
r=sqrt((xVals_rtd-xmean).^2 + (yVals_rtd-ymean).^2);  % shift the definition of r to the center of the starshade wherever that may be.
r=r(1:end-1);  % because we did a diff to get the edge slopes.
% SBSSBSSBS

% Each piece of edge is diffnorm long.  Should be 2 mm for the 192,000
% point starshade (24 petals, each 8 m long, 8000 points per petal.

% scale by the length of the differential segments
%ZSsc = ZS * diffnorm(1);
%ZPsc = ZP * diffnorm(1);  % this is the power ratio out/in for each tiny length (2 mm) of starshade edge.
ZSsc = ZS ;%  use each individual diffnorm when computing fluxedge.  SBS 102319
ZPsc = ZP ;
ZTsc = ZSsc + ZPsc;  % total power ratio out/in for both polarizations, into telesocpe which has 1 m^2 aperture (that's set in scatterometerScaleSingle.m)
% same for stealth
%ZSscst = ZSst * diffnorm(1);
%ZPscst = ZPst * diffnorm(1);  % this is the power ratio for each tiny length (2 mm) of starshade edge.
ZSscst = ZSst ;%  use each individual diffnorm when computing fluxedge.   SBS 102319
ZPscst = ZPst ;
ZTscst = ZSscst + ZPscst;  % total power ratio for both polarizations, into telescope which has 1 m^2 aperture (that's set in scatterometerScaleSingle.m)
a = find(theta>-90 & theta < 90 & r>opt.slr_glnt.r_m); % angles outside this range are 'trailing edge' illumination
npts = length(theta);
fluxedge = zeros(1,npts);
fluxedgest=fluxedge;
% Given theta and phi for each point on the edge (theta varies around the edge, phi is fixed), index points us to the right spot in the 2D array.
% Recall from the comments at the top of the code: "ZS, ZP are 1801 x 701. First dimension is -90 to 90 deg in 0.1 deg
% steps. Second dimension is 15-85 deg in 0.1 deg steps."
indx = round(((theta(a))+90)*10 + 1) + 1801*round((phi-15)*10);
%fluxedge(a) = ZTsc(indx); % It should be fluxedge(a)=ZTsc(indx) * 2e-3;
fluxedge(a) = ZTsc(indx).*diffnorm(a);  % Given theta and phi for each point on the edge (theta varies around the edge, phi is fixed), index points us to the right spot in the 2D array.
fluxedgest = fluxedge ;
% stealth edges
a=find(abs(theta)<stealth_exclusion_angle & r>opt.slr_glnt.r_m);
indx = round(((theta(a))+90)*10 + 1) + 1801*round((phi-15)*10);
%fluxedgest(a)=ZTscst(indx); % It should be fluxedgest(a)=ZTscst(indx) * 2e-3;
fluxedgest(a)=ZTscst(indx) .* diffnorm(a);
% stealth_length = length(a)*mean(diffnorm(a)); % Not used.

% Let's bin the locus into the pixel pitch of the image (the PSF basis has already been rebinned from opt.px_psf_mas to opt.px_scn_mas in convolve_with_one_wavelength).
xmas = opt.px_scn_mas * ( round( xVals_rtd( 1 : npts ) / opt.dst_strshd_tlscp_m * 180 / pi * 3600 * 1000 / opt.px_scn_mas ) ) ;  % This sets the size of the starshade in the gridded arrays.
ymas = opt.px_scn_mas * ( round( yVals_rtd( 1 : npts ) / opt.dst_strshd_tlscp_m * 180 / pi * 3600 * 1000 / opt.px_scn_mas ) ) ;
% The PSF has to be bigger than the range of xmas/ymas
hlf_psf_1 = ( sz_scn_cnv( 1 ) - 1 ) / 2 ;
hlf_psf_2 = ( sz_scn_cnv( 2 ) - 1 ) / 2 ;
fluxgrid = zeros( sz_scn_cnv ) ;
fluxgridst = fluxgrid;
  for jj=1:npts
  xindx = round( ( xmas( jj ) + hlf_psf_2 * opt.px_scn_mas ) / opt.px_scn_mas + 1 ) ;
  yindx = round( ( ymas( jj ) + hlf_psf_1 * opt.px_scn_mas ) / opt.px_scn_mas + 1 ) ;
  fluxgrid(yindx,xindx)=fluxgrid(yindx,xindx)+fluxedge(jj);   % ratio of W/m^2 out/in
  fluxgridst(yindx,xindx)=fluxgridst(yindx,xindx)+fluxedgest(jj);
  end
