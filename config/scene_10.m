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
function opt = scene_10( opt )
% In this example, we will show how to select the positions of the planets directly with the Keplerian parameters.
%
% Parameters shown below are the subset of all possibilities used in this run.
% Hint: exit debug mode (K> in a matlab session) typing 'dbquit all'
%
% Run it as:
% clear opt;opt.run='scene_10';sister(opt); % Blue band
% clear opt;opt.run='scene_10';opt.kepler.inclination_rad=75/180*pi;sister(opt);
% 

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = '' ;
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
% Jitter of the telescope alignment
opt.jitter.do = 1 ;
opt.jitter.rms_mas = 15 ;

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
% Pixel scale of the camera. Default:
% opt.pix_camera_mas = 0.4 * opt.lambda_imaging_1_nm * 1e-9 / opt.diameter_telescope_m * 180 / pi * 3600e3 ;

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%

opt.noise.do = 0 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 1.5 * 3600 ; % Default 3600 sec

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0
% Field of View of the astrophysical scene in mas
opt.scene.fov_diam_mas = 3000 ; % Default is 5000 mas

%%%%%%%%
% Star %
%%%%%%%%

% Host star from ExoCat. The ExoCat excel file is found in the 'input_scenes' folder. Any HIP identifier can be used. Some of them have some 
% particular names ('COMMON' in the excel sheet): 'tau ceti' (3.65 pc, 0.52 L_Sun), 'alf Arietis' (20.18 pc, 85.89 L_Sun), '82 Eridani' (6.04 pc, 0.69 L_Sun), 'epsilon Eridani' (3.21 pc, 0.35 L_Sun), alf Cen A (1.24 pc, 1.61 L_Sun), ... The code is not case-sensitive for the star name.
opt.star.name = 'HIP 8102' ; % equivalent to 'tau ceti'

%%%%%%%%%%%
% Planets %
%%%%%%%%%%%

% Do Kepler orbits (PS: if set on, overwrites the planets.add options set before)
opt.kepler.do = 1 ; 
% Example using Solar System planets
opt.kepler.system = 'solar' ;
% Optional variations of the solar system (NaN will leave the value from the Solar System)
opt.kepler.phase_function = { 'Lambertian' } ;
% Orbital inclination with respect the invariant plane of the extra-solar system. Make sure it is consistent with the inclination of the exozodiacal light, not necessarily the same.
  if ~isfield( opt.kepler, 'inclination_rad' )
  opt.kepler.inclination_rad = 0 ;
  end
% Moving the planets along their orbits
opt.kepler.argument_periapsis_rad = [ 0,  pi / 2,  pi,  pi / 4, 3 * pi / 4, 0, 0 ] ;
opt.kepler.longitude_ascending_node_rad = pi / 2 ; % Default is pi/2, horizontal orbit.
% Starting date for the orbital motion
opt.kepler.time_year.initial = 2021 ; % Default 0h00m00s 01/01/2025
% End date for the orbital motion
opt.kepler.time_year.final = 2021 ; % Default 0h00m00s 01/01/2030. Hint: set it to the same value as time_year.initial to get a single image to check that the display is the one intended for the whole movie (see below 'plot options').
% Time ellpased between two simulations
opt.kepler.time_year.interval = 1 / 87.65432 ; % Default 1 year (careful with time intervals that are a multiple of the period of some planet since it would repeat the location of that planet in the following image). 

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
  if ~isfield( opt, 'plot' )
  opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
  end
% Signal+noise (see starshade_band_imaging.m, % Combinations of signal and noise, for other options. Multiple options are fine: e.g., [ 1, 3, 4 ]. 3 figures will be created )
opt.plot.combination_list = 1 ;
% Whether to add a label for the planets
opt.plot.add_planet_label = 1 ;
  if ( opt.plot.add_planet_label )
  % Optional labels for planets (as many as planets)
  opt.plot.planet_label = { 'V', 'E', 'M', 'J', 'S', 'U', 'N' } ; % Default B, C, ...
  % Positional shift for the label of each planet wrt the position of the planet (times FWHM)
  opt.plot.planet_label_pos_fwhm = [ 1.5, 1.5, 0, 2.5, -2, 0.2, 0.2 ] ; % Default all 1
  end
opt.plot.min = -10 ;
opt.plot.max = 20 ;
% Zoom on the FOV to be plotted
opt.plot.zoom_list = 1 ; % Default is 1 (no zooom in)
% Video (automatically saves the individual images)
opt.video.do = 0 ; % If the same initial and final epochs are used, there's no video.
% Storing the (individual) images (if a video is stored, one may likely choose not to store the individual images)
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
