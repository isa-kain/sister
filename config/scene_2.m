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
function opt = scene_2( opt )
% Starlight only from a sun-like star at 5 pc with a non-ideal (perturbed) starshade. 
%
% Parameters shown below are the subset of all possibilities used in this run.
% Hint: exit debug mode (K> in a matlab session) typing 'dbquit all'
%
% Run it as:
% clear opt ; opt.run = 'scene_2' ; sister( opt ) ; % Blue band
% clear opt ; opt.lambda_imaging_1_nm = 615 ; opt.lambda_imaging_2_nm = 800 ; opt.run = 'scene_2' ; sister( opt ) ; % Green band
%

%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
% opt.starshade.mode='non-spinning' ;  % Default is 'spinning'

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% For NI2 (WFIRST), SISTER takes into account automatically the collecting area blocked by the secondary (0.32x diameter of the primary mirror)
% and the struts (7.17% loss). Total loss of 17.41%. For TV3 (HabEx), SISTER takes into account the secondary mirror (450 mm. The primary is 4000 mm).
% There's no specs about the struts in the HabEx public report of 09/2019, although the effect should really be of second order since the secondary
% is blocking 1.3% only. These settings are consistent with  WFIRST and HabEx mission specs.

% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 16 ;

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%

% First wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_1_nm' )
  opt.lambda_imaging_1_nm = 425 ; % Default value depends on the telescope. 
  end
% Last wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_2_nm' )
  opt.lambda_imaging_2_nm = 552 ; % Default value depends on the telescope.
  end
% Wavelength steps to be simulated
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  opt.delta_lambda_imaging_nm = 10 ; % Set it here to a constant value, or provide an array of wavelengths (see sister_imaging.m/get_lambda_scene.m)
  end

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera (final pixel scale of the output images)
opt.pix_camera_mas = 3 ; % Default is 0.4 lambda_imaging_1_nm/D. Enough sampling for the imaging band.

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%
opt.noise.do = 0 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 24 * 3600 ; % Default 3600 sec

%%%%%%%%%%%%%%%%%%%%%%%
% Non-ideal starshade %
%%%%%%%%%%%%%%%%%%%%%%%
% Setting on/off perturbed shapes
opt.locus.do = 1 ; % Default is 0
% File containing the new locus of points. For td5 I assume it would have the same name as the nominal filename plus some label. For example: 'td5_30petal_locus_p1_5mm.mat'
opt.locus.perturbed_filename = 'NI2_test_case_1em10' ; % Default is none.
% One can set the relative precision of the non-ideal spinning PSF (see SISTER Handbook for details)
% opt.locus.perturbed_psf_precision = 1e-11 ;

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0. It would read a pre-existing scene
% Field of View of the astrophysical scene in mas
opt.scene.fov_diam_mas = 550 ; % Default is 5000 mas

%%%%%%%%
% Star %
%%%%%%%%
opt.star.type = 'Sun' ; % Default is 'Sun'
opt.star.distance_to_earth_pc = 5 ; % Default is 10 pc
% Besides the star type, it's necessary to set its brightness
opt.star.app_mag_v = 4.81 - 5*log10( 10 / opt.star.distance_to_earth_pc ) ; % Default 4.81 (Sun at 10 pc)

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
opt.plot.min = 0 ;
opt.plot.max = 300 ;
opt.save_plot = 1 ; % If the image is to be stored (by default, they are. If '0', the plot is not stored)

% Over-plotting the shape of the petals
opt.plot.starshade_petals.do = 1 ; % Default is 0, not overplotted
opt.plot.starshade_petals.n_values = 0.0025 ; % the portion of the elements of the locus array that will be overplotted. See SISTER handbook for more details.
opt.plot.starshade_petals.color = 'w' ; % Default, w(hite). Any Matlab color.


%%%%%%%%%%%%%%%%%%%%
% Storing products %
%%%%%%%%%%%%%%%%%%%%

% Saving the simulation
  if ~isfield( opt, 'save_output' )
  opt.save_output = 0 ;
  % Input/output with FITS files instead of matlab files (used in starshade_band_imaging.m)
  opt.fitsio = 0 ; % Default is 0, and matlab files are used for I/O.
  end
