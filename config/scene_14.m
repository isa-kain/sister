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
function opt = scene_14( opt )
% Example of extracting a data cube for each planet. We use a previous example where we added an exo-kuiper belt (including noise, specific exoplanets, star, exozodi, and some other detector options.
% Run it as:
% clear opt ; opt.run='scene_14'; sister( opt ) ;

% Preliminary
  if ~nargin
  opt = [] ;
  end

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = '' ;
  end

%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
opt.starshade.mode='spinning' ;  % Default is 'spinning'
% HabEx (TV3)
opt.starshade.nominal_filename = 'TV3' ;

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% Getting some of the parameters from the locus file (when they are generated at JPL, they have some more info)
opt = new_occulter_from_matlab_file( opt ) ;
% For NI2 (WFIRST), SISTER takes into account automatically the collecting area blocked by the secondary (0.32x diameter of the primary mirror)
% and the struts (7.17% loss). Total loss of 17.41%. For TV3 (HabEx), SISTER takes into account the secondary mirror (450 mm. The primary is 4000 mm).
% There's no specs about the struts in the HabEx public report of 09/2019, although the effect should really be of second order since the secondary
% is blocking 1.3% only. These settings are consistent with  WFIRST and HabEx mission specs.

% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 32 ;
% optical throughput
opt.optical_throughput = sprintf( '%sinput_scenes/HabEx_transmission_VIS.txt', sister_installation_path() ) ; % Data provided by Stefan Martin, JPL, 12/19. They may change in the future.

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%
opt.lambda_band_nm_min = 300 ;
opt.lambda_band_nm_max = 1000 ;
% For HabEx, this PSF basis was precomputed in steps of 50 nm
opt.delta_lambda_nm = 50 ;
% Pixel scale of the PSF objects. It should be the same or less than the pixel scale of the scene. 
  if ~isfield( opt, 'px_psf_mas' )
  opt.px_psf_mas = 3 ; % Default is 1 mas
  end
% Extension of the PSF on the image plane
  if ~isfield( opt, 'n_lambda_over_d' )
  opt.n_lambda_over_d = 7 ; % default is 7
  end

% The distance from which the PSF becomes stationary. For WFIRST (NI2) and HabEx (NW2) there's a default value for each band coming from some analysis. For a new occulter, the user has to provide the value based on some own analysis.
opt.r_stationary_mas = 150 ;

% First wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_1_nm' )
  opt.lambda_imaging_1_nm = 450 ; % Default value depends on the telescope. 
  end
% Last wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_2_nm' )
  opt.lambda_imaging_2_nm = 550 ; % Default value depends on the telescope. 
  end
% Wavelength steps to be simulated
% For this example, The exo-Kuiper was derived at 500 nm, planets have constant contrast, the star is G0V. A flat response would incur in only 0.33% error, compared with simulating many steps (ln(550/450)/(100/500)=0.33%). One can check if that's the error given that the only thing that is caled with wavelength is the exo-zodiacal light.
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  opt.delta_lambda_imaging_nm = 10 ; % Set it here to a constant value, or provide an array of wavelengths (see sister_imaging.m/get_lambda_scene.m)
  end

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera. Default:
opt.pix_camera_mas = 13 ;
% Detector QE
opt.qe = sprintf( '%sinput_scenes/HabEx_QE_VIS.txt', sister_installation_path() ) ; % % Data provided by Stefan Martin, JPL, 12/19. They may change in the future.

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%
opt.noise.do = 1 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 10 * 3600 ; % Default 3600 sec
% Time per frame in seconds
opt.noise.exp_time_frame_sec = 100 ; % Default 600 sec
% Detector's gain (conversion of photons to photo-electrons)
opt.noise.gain = 1 ; % Default is 1
% Read noise in photo-electron units per detector pixel

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0. It would read a pre-existing scene
% Field of View of the astrophysical scene in mas
opt.scene.fov_diam_mas = 12000 ; % Default is 5000 mas

%%%%%%%%%%%%%%%%%%%%%%%%
% Local zodiacal light %
%%%%%%%%%%%%%%%%%%%%%%%%
% Recall that, by default, local zodi is not added. It should be present in any astrophysical simulation. To add it, set opt.local_zodi.do=1; in the call or configuration file.
opt.local_zodi.do = 1 ; % sets whether local zodiacal light is added
% This option will not be used unless opt.local_zodi.do=1;
opt.local_zodi.mag_v_arcsec2 = 22.5 ; % The surface brightness of the local zodiacal light, expressed in V mag per arcsec^2. Default value is 23.

%%%%%%%%
% Star %
%%%%%%%%
opt.star.name = 'beta_cvn' ; %  (https://en.wikipedia.org/wiki/Beta_Canum_Venaticorum: 8.44 pc, apparent V 4.26, G0V type)
opt.star.type = 'G0V' ; % Default is 'Sun'. Non case-sensitive.
opt.star.distance_to_earth_pc = 8.44 ; % Default is 10 pc
% Besides the star type, it's necessary to set its brightness
opt.star.app_mag_v = 4.26 ; % Beta VCn has V 4.26. Default 4.81 (Sun at 10 pc, http://mips.as.arizona.edu/~cnaw/sun.html)

%%%%%%%%%%%
% Planets %
%%%%%%%%%%%

% Example of adding 'static' planets. That is, without an orbital motion
opt.planets.add.do = 1 ; % Default is 0
opt.planets.add.pos_arc_ra_mas = [ -136.9372 88.1161 -193.8981 837.8042 628.3532 ] ;
opt.planets.add.pos_arc_dec_mas = [ -107.4348 -184.1232 423.8389 837.8042 -628.3532 ] ;
opt.planets.add.flux_ratio = [ 1.49e-10 1.911e-9 3.0e-10 3.0e-10 1.60e-10 ] ;

%%%%%%%%%%%%
% Exo-zodi %
%%%%%%%%%%%%

% Times solar zodiacal light
opt.exozodi.factor = 3 ; 
% Inclination as viewed from the telescope (it does not need to be the same value as the orbital inclination)
opt.exozodi.inclination_deg = 60 ; 
% Position angle of the exoplanetary system
opt.exozodi.position_angle_deg = 135 ;

%%%%%%%%%%%%%%%%%%%
% Exo-Kuiper belt %
%%%%%%%%%%%%%%%%%%%

opt.exokuiper.do = 1 ; % 0. By default is not added.
% Kuiper exterior radius in AU
opt.exokuiper.radius_int_au = 19 ;
% Kuiper exterior radius in AU
opt.exokuiper.radius_ext_au = 31 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Astrophysical background, star position and kinematics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.extragalactic_background.add = 1 ; % Default: 0, not added
opt.extragalactic_background.pos_arc_ra_mas = -1750 ;  % default 0, same center as external extragalactic background data
opt.extragalactic_background.pos_arc_dec_mas = 2600 ; % default 0, same center as external extragalactic background data

%%%%%%%%%%%%%%%%%%%%
% SISTER Data Cube %
%%%%%%%%%%%%%%%%%%%%
opt.cube.do = 1 ; % Default is 0. SISTER does not output a FITS file with th spectral results of the simulation.
opt.cube.fits = 0 ; % Whether the output is a FITS or Matlab file. By default, a Matlab file.
opt.cube.background = 0 ; % Sets whether the background is included or subtracted. If 1, it is included. if 0, it is subtracted.
opt.cube.planets.do = 1 ; % Extract cubes for each planet
opt.cube.planets.n_fwhm = 3 ; % Area to be stored with the planet at the center in terms of FWHM.

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
% Signal+noise (see starshade_band_imaging.m, % Combinations of signal and noise, for other options. Multiple options are fine: e.g., [ 1, 3, 4 ]. 3 figures will be created )
opt.plot.combination_list = 1 ;
% Optional labels for planets (as many as planets)
opt.plot.planet_label = { 'a', 'b', 'c', 'd', 'e' } ; % Default B, C, ...
% Positional shift for the label of each planet wrt the position of the planet (times FWHM)
opt.plot.planet_label_pos_fwhm = [ 7, 7, 5, -4, 6 ] ; % Default all 1
% A maximum limit for the linear color scale (noise is the RMS of the noise contribution over all the FOV). If the image is noiseless, the plotting script
% will choose a default min/max. Alternatively, use opt.plot.min and opt.plot.max to set the range of the plot to any specified value.
opt.plot.n_sigma_noise = 10 ; % Default is 10
opt.plot.min = 0 ;
opt.plot.max = 1000 * opt.noise.exp_time_total_sec / ( 86400 ) ; % 4000 is for 24 hour integration time
% Zoom on the FOV to be plotted
opt.plot.zoom_list = 1 ;
% No title
opt.plot.title = '' ;
% No super title
opt.plot.suptitle = '' ;
% No colorbar
opt.plot.colorbar = 0 ;
% N (upwards) E (left)
opt.plot.noflip = 1 ;
% Storing the (individual) imags
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

