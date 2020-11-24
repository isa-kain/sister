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
function opt = new_occulter_from_matlab_file( opt )
% Function that fills out some basic parameters from a new occulter file in matlab form
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

  if ~isfield( opt, 'starshade' )
  disp( 'Provide the option opt.starshade.nominal_filename. Returning.' )
  return
  else
    if ~isfield( opt.starshade, 'nominal_filename' )
    disp( 'Provide the option opt.starshade.nominal_filename. Returning.' )
    return
    end
  end

% Check there is indeed a matlab file where it should be
pth_occltr = sister_installation_path() ;
fl_occltr = [ pth_occltr 'input_scenes/locus/in/' opt.starshade.nominal_filename '.mat' ] ;
  if exist( fl_occltr, 'file' ) ~= 2
  disp( sprintf( 'The matlab occulter file %s is not found.', fl_occltr ) )
  return
  else
  disp( sprintf( 'Considering the occulter %s.', fl_occltr ) )
  end

% 
% Parameters that are necessary to create the PSF basis (default values are meant to simplify the imaging interface)
% Loading the matlab occulter file
load( fl_occltr )

% The attenuation variabe must be present
  if ~exist( 'r', 'var' )
  disp( sprintf( 'The matlab occulter file %s does not contain the attenuation profile r.', fl_occltr ) )
  return
  end

% Number of petals. Default number of petals: 24
  if exist( 'numPetals', 'var' )
  opt.starshade.number_of_petals = numPetals ; 
  disp( sprintf( 'The number of petals is %i', numPetals ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable numPetals.', fl_occltr ) )
  return
  end

% Distance Starshade-telescope (see set_starshade_distance.m in get_default_options.m for some pre-defined cases)
  if exist( 'Z', 'var' )
  opt.distance_starshade_telescope_m = Z ; % meters
  % It is usually 10s of thousand km
    if ( Z < 1e7 )
    disp( sprintf( 'The distance Starshade-telescope seems to be too small (<1e4 km): %f m. Make sure it is what is expected. Returning', Z ) )
    return
    end
  disp( sprintf( 'The distance Starshade-telescope is %f m', Z ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable Z.', fl_occltr ) )
  return
  end

% Telescope diameter (main mirror). Default values only when starshade.nominal_filename is a Matlab file with these parameters defined in it (e.g., for 'NI2', 'NW2', 'TV3', or 'UH17).
  if exist( 'telescopeDiameter', 'var' )
  opt.diameter_telescope_m = telescopeDiameter ;
  disp( sprintf( 'The diameter of the telescope is %f m', telescopeDiameter ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable telescopeDiameter.', fl_occltr ) )
  return
  end

% Geometric IWA. Default values only when starshade.nominal_filename is a Matlab file with these parameters defined in it (e.g., for 'NI2', 'NW2', 'TV3', or 'UH17).
  if exist( 'occulterGeoIWA', 'var' )
  opt.geo_iwa_mas = occulterGeoIWA * 3600e3 * 180 / pi ;
  disp( sprintf( 'The geometric IWA is %f mas', opt.geo_iwa_mas ) )
  disp( sprintf( 'The diameter of the Starshade is %3.0f m', 2 * occulterGeoIWA * Z ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable occulterGeoIWA.', fl_occltr ) )
  return
  end

% Minimum wavelength of the imaging bands of the instrument. Default values only when starshade.nominal_filename is a Matlab file with these parameters defined in it (e.g., for 'NI2', 'NW2', 'TV3', or 'UH17).
  if exist( 'lowerScienceWavelength', 'var' )
  opt.lambda_band_nm_min = round( lowerScienceWavelength * 1e9 ) ;
  disp( sprintf( 'The minimum instrumental wavelength is %04.2f nm', opt.lambda_band_nm_min ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable lowerScienceWavelength.', fl_occltr ) )
  return
  end

% Maximum wavelength of the imaging bands of the instrument. Default values only when starshade.nominal_filename is a Matlab file with these parameters defined in it (e.g., for 'NI2', 'NW2', 'TV3', or 'UH17).
  if exist( 'upperScienceWavelength', 'var' )
  opt.lambda_band_nm_max = round( upperScienceWavelength * 1e9 ) ;
  disp( sprintf( 'The maximum instrumental wavelength is %04.2f nm', opt.lambda_band_nm_max ) )
  else
  disp( sprintf( 'The matlab occulter file %s does not contain the variable upperScienceWavelength.', fl_occltr ) )
  return
  end

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = opt.starshade.nominal_filename ;
  end

