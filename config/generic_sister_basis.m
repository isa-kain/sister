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
function opt = generic_sister_basis( opt )
% 
% 1) Parameters that are necessary to create the PSF basis 
% opt.starshade.mode (also for imaging). Default: 'spinning'.
% opt.starshade.nominal_filename (also for imaging). Default: 'NI2'.
% opt.starshade.number_of_petals (also for imaging). Default: 24.
% opt.distance_starshade_telescope_m. Default values only when starshade.nominal_filename has the charcacters 'NI2' or 'NW2'.
% opt.diameter_telescope_m. Default values only when starshade.nominal_filename has the characters 'NI2' or 'NW2'.
% opt.pupil_filename. Default: ideal (means circular). If it is NI2, it defaults to pupil_D1Kpix that corresponds to WFIRST, including struts)
% PS: opt.pupil_filename is any Matlab or FITS file with a 2D array defining the pupil. The file is to be stored in 'input_scenes/locus/in'.
% opt.secondary_size. Default: 0, no secondary if the pupil is ideal. It is the ratio between the diameter of the secondary and primary mirrors.
% opt.Nx_pupil_pix. Default: 64.
% opt.lambda_band_nm_min. Minimum wavelength of the imaging bands of the instrument. Default values only when starshade.nominal_filename has the charcacters 'NI2' or 'NW2'.
% opt.lambda_band_nm_max. Maximum wavelength of the imaging bands of the instrument. Default values only when starshade.nominal_filename has the charcacters 'NI2' or 'NW2'.
% opt.delta_lambda_psf_nm. The wavelength spacing used to generate the PSF basis. Default value is 10 nm.
% opt.psf_spacing_mas. Minimum spacing between two parallel rays (each ray gives rise to a PSF. See sister_basis.m). Default: 1.
% opt.px_psf_mas. The pixel scale of the PSF basis. Default: 1.
% opt.n_lambda_over_d. The size of the two-dimensional array of the PSF on the image plane that is stored is defined as twice the number of times of lambda/D in both directions. Default: 7. So that the PSF size is +/- 7 lambda/D in each direction. 
% opt.r_stationary_mas. Default values only when starshade.nominal_filename has the characters 'NI2' or 'NW2'.
% opt.geo_iwa_mas. Default. Default values only when starshade.nominal_filename has the characters 'NI2' or 'NW2'.
%
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
% File containing the x,y values of the Starshade petals in meters. If z values are not provided, they are set to 0. This file is to be found in /input_scenes/locus/in/.
% Load NI2_24petal.mat to see an example of xVals, yVals.
opt.starshade.nominal_filename = 'some_xy_occulter' ;
% Default number of petals: 24
opt.starshade.number_of_petals = 24 ;
% Distance Starshade-telescope (see set_starshade_distance.m in get_default_options.m for some pre-defined cases)
opt.distance_starshade_telescope_m = 3.724225668350351e+07 ; % meters

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% Diameter of the (primary mirror) of the telescope in meters
opt.diameter_telescope_m = 1 ; 
% Any Matlab or FITS file with a 2D array defining the pupil. 
opt.pupil_filename = 'ideal' ; % Ideal pupil (no struts. If NI2, there's a secondary of 0.32x primary diameter, consistent with WFIRST. If TV3, the value is 0.1125, consistent with HabEx)
% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 16 ;
%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%
% When introducing a new occulter, one has to provide the different imaging bands in order to (i) create the necessary PSF basis, and (ii) check that the imaging wavelength, defined thereafter, is in-band. If more than one imaging band, use a vector array to define them [ a, b, c ]. The cases of WFIRST and HabEx have been added in starhade_imaging.m (get_imaging_options.m) as default cases when the occulter has NI2 and NW2, respectively, in their names.
opt.lambda_band_nm_min = 400 ; 
opt.lambda_band_nm_max = 470 ; 
% The wavelength spacing used to generate the PSF basis
opt.delta_lambda_psf_nm = 10 ; 
% Pixel scale of the PSF basis (usually 3 mas is good enough. It depends on the scale of the camera's pixels)
opt.px_psf_mas = 3 ; 
% Minimum spacing between two parallel rays (each ray gives rise to a PSF. See sister_basis.m)
opt.psf_spacing_mas = 1 ;
% Extension of the PSF on the image plane
opt.n_lambda_over_d = 7 ;

% The distance from which the PSF becomes stationary. For WFIRST (NI2) and HabEx (NW2) there's a default value for each band coming from some analysis. For a new occulter, the user has to provide the value based on some own analysis.
opt.r_stationary_mas = 140 ; 

% The IWA has to be set. For WFIRST (NI2) and HabEx (NW2) there's a default value for each band.
opt.geo_iwa_mas = 30 ;

