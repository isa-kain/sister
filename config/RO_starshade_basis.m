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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = RO_starshade_basis( opt )
% Generates a PSF basis for WFIRST and a spinning starshade. See below to find some more specific options.  

  if ~exist( 'opt', 'var' )
  opt = [ ] ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%
% PSF BASIS AND IMAGING %
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% Starshade %
%%%%%%%%%%%%%
% Spinning/Non-spinning starshade mode.
opt.starshade.mode = 'spinning' ; % 'spinning', 'non-spinning'
opt.starshade.nominal_filename = 'UH17' ;
opt.starshade.number_of_petals = 48 ;
opt = new_occulter_from_matlab_file( opt ) ;

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% Pupil. Ideal if it is spinning, including the secondary and struts. Struts are considered in a simple way in this case: the blocked area of the pupil is a secondary with equivalent area as the secondary+struts. See get_default_options.m when setting opt.secondary for details.
% If it is non-spinning, it uses input_scenes/pupil_D1Kpix_256.fits which is for WFIRST including secondary and struts in the optical model and when computing the collecting area.
% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
% opt.Nx_pupil_pix = 100 ; % original run
opt.Nx_pupil_pix = 16 ;
% Wavelength band
opt.lambda_band_nm_min = 400 ;
opt.lambda_band_nm_max = 700 ;
opt.delta_lambda_psf_nm = 10 ;
% Pixel scale of the PSF basis
opt.px_psf_mas = 3 ;
% Extension of the PSF on the image plane
opt.n_lambda_over_d = 7 ;

opt.diameter_telescope_m = 39; % ELT diameter
opt.distance_starshade_telescope_m = 195000000; % minimum altitude during observations taken from BAAS paper
opt.r_stationary_mas = 128 ; % taken from scene_15? not sure what analysis backs up this number.

% The IWA has to be set. For WFIRST (NI2) is 72 mas. For HabEx (TV3) is 70 mas.
opt.geo_iwa_mas = 58.34 ;


