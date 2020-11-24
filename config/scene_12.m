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
function opt = scene_5( opt )
% In this example, we show how to include a non-ideal starshade and compare it with an ideal one (scene_5.m). We image an external astrophysical scene, with local zodiacal light in some direction, and with extragalactic objects (from https://asd.gsfc.nasa.gov/projects/haystacks/downloads.html)
%
% Parameters shown below are the subset of all possibilities used in this run.
% Hint: exit debug mode (K> in a matlab session) typing 'dbquit all'
%
% Run it as:
% clear opt;opt.run={'scene_12','scene_5'};opt.noise.do=1;opt.plot.min=700;opt.plot.max=1300;sister(opt); % Green band
% 

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = 'haystacks_with_background' ;
  end

% Important remark: if you are changing the configuration file *and* storing the matlab data (opt.save_output=1, found at the end), 
% set redo equal to 1 to make sure the imaging is processed again. Otherwise it will use the results obtained before.
% opt.redo=1 ; % Default 0

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
  opt.lambda_imaging_1_nm = 615 ; % Default value depends on the telescope.
  end
% Last wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_2_nm' )
  opt.lambda_imaging_2_nm = 800 ; % Default value depends on the telescope.
  end
% Wavelength steps to be simulated
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  opt.delta_lambda_imaging_nm = 10 ; % Set it here to a constant value, or provide an array of wavelengths (see sister_imaging.m/get_lambda_scene.m)
  end

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera. Default:
% opt.pix_camera_mas = 0.4 * opt.lambda_imaging_1_nm * 1e-9 / opt.diameter_telescope_m * 180 / pi * 3600e3 ;

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%
  if ~isfield( opt, 'noise' )
  opt.noise.do = 0 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
  end
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 24 * 3600 ; % Default 3600 sec
% Time per frame in seconds
opt.noise.exp_time_frame_sec = 400 ; % Default 600 sec

%%%%%%%%%%%%%%%%%%%%%%%
% Non-ideal starshade %
%%%%%%%%%%%%%%%%%%%%%%%
% Setting on/off perturbed shapes
opt.locus.do = 1 ; % Default is 0
% File containing the new locus of points. For td5 I assume it would have the same name as the nominal filename plus some label. For example: 'td5_30petal_locus_p1_5mm.mat'
opt.locus.perturbed_filename = 'NI2_test_case_1em10' ; % Default is none.
opt.locus.perturbed_psf_precision = 1e-11 ;

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%

% Reading an external scene. This example comes from the Haystacks Project (https://asd.gsfc.nasa.gov/projects/haystacks/haystacks.html)
opt.scene.name = 'modern_cube_zodi1inc60dist10' ;

%%%%%%%%%%%%%%%%%%%%%%%%
% Local zodiacal light %
%%%%%%%%%%%%%%%%%%%%%%%%
opt.local_zodi.do = 1 ; % Sets whether local zodiacal light will be added to the scene. Default is 0, not added.
% Deriving its surface brightness for some line of sight
opt.local_zodi.helio_eclipitic_longitude_deg = 134 ; % The helio-ecliptic longitude in degrees. Data are derived from STScI (see SISTER handbook). The value of the surface brightness can be read off opt.local_zodi.mag_v_arcsec2 after the run
opt.local_zodi.eclipitic_latitude_deg = 62 ; % Notice it does not have helio in the name. The ecliptic latitude in degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Astrophysical background, star position and kinematics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.extragalactic_background.add = 1 ; % Default: 0, not added
opt.extragalactic_background.pos_arc_ra_mas = -1750 ;  % default 0, same center as external extragalactic background data
opt.extragalactic_background.pos_arc_dec_mas = 2600 ; % default 0, same center as external extragalactic background data
% For proper motion and parallax (only used if there's an extragalactic background field set on)
% Position of the exo-planetary system
opt.star.ra_deg  = 80 ; % Default is 0 deg
opt.star.dec_deg = 60 ; % Default is 0 deg
% Proper motion (parallax is defined from the distance to the star)
opt.star.pm_ra_mas_yr  = 150 ; % Default is 0, no proper motion
opt.star.pm_dec_mas_yr = -100 ; % Default is 0, no proper motion

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
  if ~isfield( opt, 'plot' )
  opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
  end
% Signal+noise (see starshade_band_imaging.m, % Combinations of signal and noise, for other options. Multiple options are fine: e.g., [ 1, 3, 4 ]. 3 figures will be created )
opt.plot.combination_list = 2 ;
% A maximum limit for the linear color scale (noise is the RMS of the noise contribution over all the FOV). Alternatively, use opt.plot.min and opt.plot.max to set the range of the plot to any specified value.
opt.plot.n_sigma_noise = 10 ; % Default is 10
% Zoom on the FOV to be plotted
opt.plot.zoom_list = [ 1, 3 ] ; % Default is 1 (no zooom in)
% Storing the (individual) images
opt.save_plot = 1 ; % If the image is to be stored (by default, they are. Set it to 0 not to store them.)

%%%%%%%%%%%%%%%%%%%%
% Storing products %
%%%%%%%%%%%%%%%%%%%%

% Saving the simulation
  if ~isfield( opt, 'save_output' )
  opt.save_output = 0 ;
  % Input/output with FITS files instead of matlab files (used in starshade_band_imaging.m)
  opt.fitsio = 0 ; % Default is 0, and matlab files are used for I/O.
  end
