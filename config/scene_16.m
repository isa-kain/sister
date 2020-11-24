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
function opt = scene_16( opt )
% Same scene as scene_1, but showing explicitly all necessary parameters to run a general starshade-telescope configuration.
% Reminder: secene_1 is the starlight only from a sun-like star at 5 pc. 
%
% Run it as:
% clear opt ; opt.run = 'scene_16' ; sister( opt ) ; % Blue band

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = 'test' ;
  end

%% All the parameters defined in the imaging simulation must be consistent with the ones used to generate the PSF basis %%
%% That is, SISTER will only check whether the geometric IWA, the sampling of pupil and the pass band correspond to the PSF basis, but it will *not* check whether the diameter of the telescope, the distance starshade-telescope, are consistent. It is the responsibility of the user to make sure PSF basis and imaging simulations correspond to the same instrumental specs. %%


%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
% opt.starshade.mode='non-spinning' ;  % Default is 'spinning'
opt.starshade.nominal_filename = 'NI2' ; 

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% For NI2 (WFIRST), SISTER takes into account automatically the collecting area blocked by the secondary (0.32x diameter of the primary mirror)
% and the struts (7.17% loss). Total loss of 17.41%. For TV3 (HabEx), SISTER takes into account the secondary mirror (450 mm. The primary is 4000 mm).
% There's no specs about the struts in the HabEx public report of 09/2019, although the effect should really be of second order since the secondary
% is blocking 1.3% only. These settings are consistent with  WFIRST and HabEx mission specs.

opt.pupil_filename = 'ideal' ; % 'ideal' means circular with an optional secondary. It can also be any Matlab or FITS file with a 2D array defining the pupil.
% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 16 ;
% Diameter of the telescope in meters
opt.diameter_telescope_m = 2.4 ;
% Ratio of the diameter of the secondary mirror to the main mirror (0 means there's no obscuration due tot he secondary)
opt.secondary_size = 0.417257964788079 ; % WFIRST value is 0.32 but accounting for the struts, which cover almost an 7.17% of the collecting area, this is the effective value of an equivalent secondary without struts. % https://wfirst.ipac.caltech.edu/sims/Param_db.html#telescope
% The distance between the Starshade and the telescope in meters
opt.distance_starshade_telescope_m = 3.7242257e+07 ;
% The geometric IWA in milli arcsec (diamater of the starshade / 2 / distance starshade-telescope)
opt.geo_iwa_mas = 72 ;
% Optical throughput (it may be a file with wavelength and throughput, or a constant value).
opt.optical_throughput = 0.56 ;
% QE of the detector (it may be a file with wavelength and throughput, or a constant value). 
opt.qe = sprintf( '%sinput_scenes/EMCCD201_20_QE.txt', sister_installation_path() ) ; ; % PS: 0.7851 is an effective value across the blue band 
% Radius for the stationary region in milli arcsec.
opt.r_stationary_mas = 150 ;

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%
% Instrument pass band:
opt.lambda_band_nm_min = 425 ;
opt.lambda_band_nm_max = 552 ;
% The wvelength spacing used to generate accurate PSF responses
opt.delta_lambda_psf_nm = 10 ;

% Imaging wavelengths (anything between the instrumental pass band)
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
opt.star.app_mag_v = 4.81 - 5*log10( 10 / opt.star.distance_to_earth_pc ) ; % Default 4.81 (Sun at 10 pc, http://mips.as.arizona.edu/~cnaw/sun.html)

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
% Overplotting a circle around the starshade tips
opt.plot.starshade_circle.do = 1 ; % Default is 0, not overplotted
opt.plot.starshade_circle.color = 'w' ; % Default, w(hite). Any Matlab color.
opt.plot.starshade_circle.line_width = 1.5 ; % Default, 1.5. Any value. 1.5 may be good for dots, and 1 for lines.
opt.plot.starshade_circle.line_style = ':' ; % Default, dots. Any Matlab option: '-',  '--',  ':',  '-.',  or 'none'.

% Storing the plot
opt.save_plot = 1 ; % Default value 1: they are stored. If 0, the plot is not stored.

%%%%%%%%%%%%%%%%%%%%
% Storing products %
%%%%%%%%%%%%%%%%%%%%

% Saving the simulation
  if ~isfield( opt, 'save_output' )
  opt.save_output = 0 ;
  % Input/output with FITS files instead of matlab files (used in starshade_band_imaging.m)
  opt.fitsio = 0 ; % Default is 0, and matlab files are used for I/O.
  end

