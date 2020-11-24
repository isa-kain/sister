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

function [ scn_dt_cnv_ph_s, scn_dt_cnv_no_plnt_ph_s, opt_img, fl_dt_out ] = convolve_with_one_wavelength( i_bnd, i_lmbd, opt_img )
% Convolution of a scene with the Starshade PSF
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com
% Don't convolve the scene if it exists, or delete the output files if required. Only using tg_nm_1 to recycle the individual long convolutions if the diffuse part does not change.
  if ~isdir( sprintf( '%ssingle_wavelength/', opt_img.out_dr ) )
  mkdir( sprintf( '%ssingle_wavelength/', opt_img.out_dr ) ) ;
  end
% Whether this is from a config file, or not:
  if isfield( opt_img, 'case' )
    if isnumeric( opt_img.case )
    fl_out = sprintf( '%ssingle_wavelength/sister_%i_%snm', opt_img.out_dr, opt_img.case, strrep( sprintf( '%04.2f', opt_img.lmbd_tmp_nm ), '.', 'p' ) ) ;
    end
    if ischar( opt_img.case )
    fl_out = sprintf( '%ssingle_wavelength/sister_%s_%snm', opt_img.out_dr, opt_img.case, strrep( sprintf( '%04.2f', opt_img.lmbd_tmp_nm ), '.', 'p' ) ) ;
    end
  else
    if isfield( opt_img, 'add_label' ) && length( opt_img.add_label )
    fl_out = sprintf( '%ssingle_wavelength/starshade_%s_%s_%s_%s_%snm', opt_img.out_dr, opt_img.starshade.mode, opt_img.scn.nm, opt_img.tg_nm_1, opt_img.add_label, strrep( sprintf( '%04.2f', opt_img.lmbd_tmp_nm ), '.', 'p' ) ) ;
    else
    fl_out = sprintf( '%ssingle_wavelength/starshade_%s_%s_%s_%snm', opt_img.out_dr, opt_img.starshade.mode, opt_img.scn.nm, opt_img.tg_nm_1, strrep( sprintf( '%04.2f', opt_img.lmbd_tmp_nm ), '.', 'p' ) ) ;
    end
  end
fl_dt_out = [ fl_out '.mat' ] ;

% Get the stationary region radius and geometric IWA
opt_img.r_st_mas = set_r_stationary_mas( opt_img ) ;
% Allowing the user to define a different stationary radius
  if isfield( opt_img, 'r_stat_mas_new' )
  % Check (in pinciple, the stationary radius set in set_r_stationary_mas is the one used to build the imaging basis. If the new one set by the user is greater than that, the basis will not be constructed.
    if ( opt_img.r_stat_mas_new > opt_img.r_st_mas )
    disp( sprintf( 'The new value for the stationary radius %3.2f mas is greater than the one used to build the basis %3.2f mas. If this new value is correct, one has to create the missing elements of the imaging basis by modifying set_r_stationary_mas and running sister_basis. Stopped', opt_img.r_stat_mas_new, opt_img.r_stat_mas ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  opt_img.r_st_mas = opt_img.r_stat_mas_new ;
  disp( sprintf( 'New stationary radius set to %3.2f mas', opt_img.r_st_mas ) )
  end

% UNIT CONVERSION. IT INCLUDES QE, MIRROR AREA, AND OPTICAL THROUGHPUT.
opt_img = spectral_flux_density_to_photons_s_f( opt_img, i_lmbd, i_bnd ) ;

  if ( do_cnv_f( fl_dt_out, opt_img ) )
  % Get the scene data and the star flux that will be helpful when adding some effects like the perturbed locus
  [ scn_dt opt_img ] = get_scene_data( opt_img, i_lmbd ) ;

  % Stationary convolution
  % if the scene has an even number of pixels add an array of 0s
  % Sentinel
  sz_1_odd = 1 ;
  sz_2_odd = 1 ;
    if size( scn_dt, 1 ) / 2 == floor( size( scn_dt, 1 ) / 2 )
    sz_1_odd = 0 ;
    scn_dt( end + 1, : ) = 0 ;
    end
    if size( scn_dt, 2 ) / 2 == floor( size( scn_dt, 2 ) / 2 )
    sz_2_odd = 0 ;
    scn_dt( :, end + 1 ) = 0 ;
    end

  % Create a hole for the non-stationary region
  sz_1 = size( scn_dt, 1 ) ;
  sz_2 = size( scn_dt, 2 ) ;
  % These must be integers by construction
  cntr_1 = ( sz_1 + 1 ) / 2 ;
  cntr_2 = ( sz_2 + 1 ) / 2 ;
  % Grid of points defining the non-stationary region
  x_1 = ( -ceil( opt_img.r_st_mas / opt_img.px_scn_mas ) : 1 : ceil( opt_img.r_st_mas / opt_img.px_scn_mas ) ) ;
  x_2 = ( -ceil( opt_img.r_st_mas / opt_img.px_scn_mas ) : 1 : ceil( opt_img.r_st_mas / opt_img.px_scn_mas ) ) ;
  % Temporary array
    if ( ( cntr_1 + x_1( 1 ) < 1 ) || ( cntr_2 + x_2( 1 ) <  1 ) )
    disp( sprintf( '(WARNING) The FOV does not cover the full geometric IWA. Set fov_diam_mas to at least %2.0f. Returning', ceil( ( max( [ abs( x_1( 1 ) ), abs( x_2( 1 ) ) ] ) - 1 ) * 2 * opt_img.px_scn_mas ) + 1 ) ) ;
    scn_dt_cnv_ph_s = 0 ;
    scn_dt_cnv_no_plnt_ph_s = 0 ;
    return
    end
  V_TMP = scn_dt( cntr_1 + x_1( 1 ) : cntr_1 + x_1( end ), cntr_2 + x_2( 1 ) : cntr_2 + x_2( end ) ) ;
  % Meshgrid
  [ X_1 X_2 ] = meshgrid( x_1, x_2 ) ;
  q_tmp = sqrt( X_1.^2 + X_2.^2 ) < ( opt_img.r_st_mas / opt_img.px_scn_mas ) ;
  % Sentinel value used to identify non-stationary pixels (the input scene should not have negative values)
  sntnl_non_st = -100 ; % Don't use zero, since some other values on the scene may be zero, and also because -sentinel is used to mark pixels already convolved in the non-stationary region later on.
  V_TMP( q_tmp ) = sntnl_non_st ; % sentinel value
  % Inserting the non-stationary region and getting its pixel positions
  % In order to save a full copy of the scene array, one can just keep the non-stationary region. That requires introducing now some non-stationary variables to be used later on.
  % This ensures that for every wavelength the PSF gets properly normalized in the next calls
  opt_img.psf_st_nrm = 1 ;
    if strcmp( opt_img.starshade.mode, 'spinning' )
    % In the spinning case, one can get the whole PSF basis already
    [ psf_arry opt_img ] = get_normalized_spinning_psf_array( i_bnd, opt_img ) ;
    % Stationary PSF: the farthest psf (already within the stationary region)
    psf_st = squeeze( psf_arry( end, :, : ) ) ;
    else
    % This case stores the norm of the PSF because it corresponds to the stationary PSF (same location as the spinning one)
    % Compute the stationary PSF if there are any data in the stationary region
    tic
    [ psf_st opt_img ] = get_normalized_psf( i_bnd, opt_img, 0, opt_img.r_st_mas, 1 ) ;
    disp( 'Time spent building the stationary PSF:' )
    toc
    end

  % Storing the FWHM of the stationary PSF
  % Approximate FWHM of the PSF
  fwhm_px_tmp = opt_img.lmbd_tmp_nm * opt_img.units_wavelength_si / opt_img.dmtr_tlscp_m * 180 * 3600e3 / pi / opt_img.px_scn_mas ;
  cntr_psf_st = ( size( psf_st, 1 ) + 1 ) / 2 ;
  psf_st_tmp = psf_st( round( cntr_psf_st - 2 * fwhm_px_tmp ) : round( cntr_psf_st + 2 * fwhm_px_tmp ), round( cntr_psf_st - 2 * fwhm_px_tmp ) : round( cntr_psf_st + 2 * fwhm_px_tmp ) ) ;
  hlf_pwr_px = find( psf_st_tmp >= 0.5 * max( psf_st_tmp( : ) ) ) ;
  psf_blnk = psf_st_tmp * 0 ; psf_blnk( hlf_pwr_px ) = 1 ;
  opt_img.fwhm_mas( i_lmbd, : ) = [ opt_img.lmbd_tmp_nm, 2 * sqrt( sum( psf_blnk( : ) ) / pi ) * opt_img.px_scn_mas ] ;
  % Its array dimensions
  opt_img.psf_arry_px = size( psf_st ) ;

  sz_psf_st_1 = size( psf_st, 1 ) ;
  sz_psf_st_2 = size( psf_st, 2 ) ;
    if ( ( sz_psf_st_1 / 2 == floor( sz_psf_st_1 / 2 ) ) || ( sz_psf_st_2 / 2 == floor( sz_psf_st_2 / 2 ) ) )
    disp( 'One of the sizes of the PSF array is even and both should be odd. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  side_psf_st_1 = ( sz_psf_st_1 - 1 ) / 2 ;
  side_psf_st_2 = ( sz_psf_st_2 - 1 ) / 2 ;
  % Notice that the non-stationary PSF affects values farther than the non-stationary region, although only pixels inside the non-stationary region are a peak for a non-stationary PSF
  % Length must be at least the non-stationary radius plus the length of the sides of the PSF.
  n_non_st_1 = ceil( 1.2 * opt_img.r_st_mas / opt_img.px_scn_mas ) + side_psf_st_1 ;
  n_non_st_2 = ceil( 1.2 * opt_img.r_st_mas / opt_img.px_scn_mas ) + side_psf_st_2 ;
  % Copy for the non-stationary convolution
  idx_non_st_1_1 = cntr_1 - n_non_st_1 ;
  idx_non_st_1_2 = cntr_1 + n_non_st_1 ;
  idx_non_st_2_1 = cntr_2 - n_non_st_2 ;
  idx_non_st_2_2 = cntr_2 + n_non_st_2 ;
  % For small FOV these limits may be off. Padding the original scene with 0s. After convolution, the scene is cut back to its original size (l~340, if ( l_1_1 + l_1_2 + l_2_1 + l_2_2 )).
l_1_1 = 0 ; l_1_2 = 0 ; l_2_1 = 0 ; l_2_2 = 0 ;
    if idx_non_st_1_1 < 1
    l_1_1 = 1 - idx_non_st_1_1 ;
    end
    if idx_non_st_1_2 > sz_1
    l_1_2 = idx_non_st_1_2 - sz_1 ;
    end
    if idx_non_st_2_1 < 1
    l_2_1 = 1 - idx_non_st_2_1 ;
    end
    if idx_non_st_1_2 > sz_2
    l_2_2 = idx_non_st_2_2 - sz_2 ;
    end
    if ( l_1_1 + l_1_2 + l_2_1 + l_2_2 )
    scn_dt_tmp = zeros( sz_1 + l_1_1 + l_1_2, sz_2 + l_2_1 + l_2_2 ) ;
    scn_dt_tmp( l_1_1 + 1 : sz_1 + l_1_1, l_2_1 + 1 : sz_2 + l_2_1 ) = scn_dt ;
    clear scn_dt
    scn_dt = scn_dt_tmp ;
    clear scn_dt_tmp
    sz_1 = size( scn_dt, 1 ) ;
    sz_2 = size( scn_dt, 2 ) ;
    % These must be integers by construction
    cntr_1 = cntr_1 + l_1_1 ;
    cntr_2 = cntr_2 + l_2_1 ;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (SRH. Important remark) Way to prove the orientation of the rotated spinning PSF is the right one.
  if ( 0 )
%  scn_dt( 335 + round( 35 * cosd( 80 ) / 3 ), 335 + round( 35 * sind( 80 ) / 3 ) ) = 111 ;
%  scn_dt( 335 + round( 35 * cosd( 80 ) / 3 ), 335 - round( 35 * sind( 80 ) / 3 ) ) = 112 ;
%  scn_dt( 335 - round( 35 * cosd( 80 ) / 3 ), 335 + round( 35 * sind( 80 ) / 3 ) ) = 113 ;
%  scn_dt( 335 - round( 35 * cosd( 80 ) / 3 ), 335 - round( 35 * sind( 80 ) / 3 ) ) = 114 ;
  end
% Set stops at imrotate and flipud, fliplr and flipud( fliplr( ) )
% Identify which case hast stopped at looking at the value of: scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt(this term)
% plt_mnmx( scn_non_st_dt ). Use the plotting tool on the window to check the value of the pixel.
% plt_mnmx( psf_non_st ) ; to see the original one
% plt_mnmx( psf_non_st_rtd ) ; or with flips, to see the actual one used that will match the orientation corresponding to the identified pixel.
% PS: the conventions in the convolution depend on the orientation used to store the spinning PSF. These are all consistent with the current PSF orientation (asymmetries look horizontal instead of vertical).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  scn_non_st_dt = scn_dt( cntr_1 - n_non_st_1 : cntr_1 + n_non_st_1, cntr_2 - n_non_st_2 : cntr_2 + n_non_st_2 ) ;
  scn_dt_cp = scn_dt ;
  scn_dt( cntr_1 + x_1( 1 ) : cntr_1 + x_1( end ), cntr_2 + x_2( 1 ) : cntr_2 + x_2( end ) ) = V_TMP ;
  % Pixels defining the non-stationary region AND where there's some non zero data
  px_non_st_sntnl = find( scn_dt == sntnl_non_st ) ;
  px_non_st_non_zr = find( scn_dt_cp( px_non_st_sntnl ) ~= 0 ) ;
  px_non_st = px_non_st_sntnl( px_non_st_non_zr ) ;
  clear scn_dt_cp
  % Storing the values of the scene in the non-stationary region
  scn_dt_non_st = scn_dt( px_non_st ) ;
  % Inserting the non-stationary region for convolution
  scn_dt( px_non_st_sntnl ) = 0 ;

  % Convolve with the stationary psf
  disp( 'Convolving the stationary region' )
  % Padding the psf to the same size as scn_dt (PS: fft of 3333x3333, clearly not 2^power, takes 0.14 s. Not caring about it)
  psf_pd = scn_dt * 0 ;
  psf_pd( cntr_1 - side_psf_st_1 : cntr_1 + side_psf_st_1, cntr_2 - side_psf_st_2 : cntr_2 + side_psf_st_2 ) = psf_st ;
  tic
  % Matlab way of doing the convolution and re-centering frequencies
  scn_dt_cnv = ifft2( fft2( psf_pd ) .* fft2( scn_dt ) ) ;
  clear psf_pd 
  scn_dt_cnv = circshift( scn_dt_cnv, [ cntr_1 cntr_2 ] ) ;
  toc

    % With zero padding (it's 15% faster, but less than 3% faster for the whole process and requires more memory
    if ( 0 )
    psf_pd_2 = zeros( 4096, 4096 ) ;
    psf_pd_2( 2048 - cntr_1 - side_psf_st_1 : 2048 + cntr_1 + side_psf_st_1, 2048 - cntr_2 - side_psf_st_2 : 2048 + cntr_2 + side_psf_st_2 ) = psf_st ;
    scn_dt_2 = zeros( 4096, 4096 ) ;
    scn_dt_2( 2048 - ( sz_1 - 1 ) / 2: 2048 + ( sz_1 - 1 ) / 2, 2048 - ( sz_2 - 1 ) / 2 : 2048 + ( sz_2 - 1 ) / 2 ) = scn_dt ;
    cnv_tmp = ifft2( fft2( psf_pd_2 ) .* fft2( scn_dt_2 ) ) ;
    clear psf_pd_2 scn_dt_2
    cnv_tmp = circshift( cnv_tmp, [ 2048 2048 ] ) ;
    scn_dt_cnv = cnv_tmp( 2048 - ( sz_1 - 1 ) / 2: 2048 + ( sz_1 - 1 ) / 2, 2048 - ( sz_2 - 1 ) / 2 : 2048 + ( sz_2 - 1 ) / 2 ) ;
    clear cnv_tmp
    end

  % Convolve in the non-stationary region
  n_non_st = numel( px_non_st ) ;
  % Assign sentinel value
  scn_dt( px_non_st ) = sntnl_non_st ;
  % Array that will accumulate the convolution due to the non-stationary region
  % Length must be at least the non-stationary radius plus the length of the sides of the PSF.
  scn_non_st_cnv = zeros( 2 * n_non_st_1 + 1, 2 * n_non_st_2 + 1 ) ;
  cntr_non_st_1 = n_non_st_1 + 1 ;
  cntr_non_st_2 = n_non_st_2 + 1 ;
  % PSF array template to translate and rotate it
  psf_non_st = zeros( 2 * n_non_st_1 + 1, 2 * n_non_st_2 + 1 ) ;
  disp( 'Convolving the non-stationary region' )
  tic


if ( 0 )
  disp( '*** TEST: convolve_with_...' )
  psf_pd_2 = scn_non_st_dt * 0 ;
  psf_pd_2( cntr_non_st_1 - side_psf_st_1 : cntr_non_st_1 + side_psf_st_1, cntr_non_st_2 - side_psf_st_2 : cntr_non_st_2 + side_psf_st_2 ) = psf_st ;
  tic
  % Matlab way of doing the convolution and re-centering frequencies
  scn_non_st_dt_2 = 0 * scn_non_st_dt ;
  d_tmp = ( sz_1 + 1 ) / 2 - n_non_st_1 - 1 ; 
    for i_px = 1 : n_non_st
    px_tmp = px_non_st( i_px ) ;
    x_px = mod( px_tmp, sz_1 ) ;
    y_px = ceil( px_tmp / sz_1 ) ;
    scn_non_st_dt_2( x_px - d_tmp, y_px - d_tmp ) = scn_non_st_dt( x_px - d_tmp, y_px - d_tmp ) ;
    end
  clear d_tmp
  scn_non_st_cnv = ifft2( fft2( psf_pd_2 ) .* fft2( scn_non_st_dt_2 ) ) ;
  scn_non_st_cnv = circshift( scn_non_st_cnv, [ cntr_non_st_1 cntr_non_st_2 ] ) ;
else
  % Depends whether it is spinning or not
    if strcmp( opt_img.starshade.mode, 'spinning' )
      for i_px = 1 : n_non_st
      px_tmp = px_non_st( i_px ) ;
      % DEC
      x_px = mod( px_tmp, sz_1 ) ;
      % RA
      y_px = ceil( px_tmp / sz_1 ) ; 
      % Start with the first quadrant (not including the center) PS: very fast loop until it finds the pixels in the first quadrant (<0.002 s)
        if ( x_px > cntr_1 ) && ( y_px > cntr_2 )
        % Find the corresponding psf (remember the first index is for the centered PSF, therefore 0 distance, and matlab starts indices at 1)
        idx_psf = 1 + round( sqrt( ( x_px - cntr_1 )^2 + ( y_px - cntr_2 )^2 ) * opt_img.px_scn_mas ) ;
        % Taking into account the spacing of the PSF basis
        idx_psf = round( idx_psf / opt_img.psf_spacing_mas ) ;
        psf_tmp = squeeze( psf_arry( idx_psf, :, : ) ) ;
        % PSF is centered in the PSF data array and now has to be on top of the pixel
        % PSF translation on the 'x' axis (PS: side_psf_st is the same as the non-stationary since the psf_arry has the same dimensions for all distances)
        r_px = round( sqrt( ( x_px - cntr_1 )^2 + ( y_px - cntr_2 )^2 ) ) ;
        % reset its values
        psf_non_st = 0 * psf_non_st ;
        psf_non_st( cntr_non_st_1 - side_psf_st_1 : cntr_non_st_1 + side_psf_st_1, cntr_non_st_2 + r_px - side_psf_st_2 : cntr_non_st_2 + r_px + side_psf_st_2 ) = psf_tmp ;
        % PSF rotation towards the 'y' axis
        psf_non_st_rtd = imrotate( psf_non_st, -180 / pi * atan2( x_px - cntr_1, y_px - cntr_1 ), 'crop', 'bilinear' ) ;
        % Assign the convolution value and coadd (+1 is necessary to have a proper centering, if x_px=cntr_1, the pixel is n_non_st_1+1, which is the center of the subarray scn_non_st_dt, which in turn is the center of scn_dt, as scn_non_st_dt was defined.
        scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( x_px - cntr_1 + n_non_st_1 + 1, y_px - cntr_2 + n_non_st_2 + 1 ) * psf_non_st_rtd ;
        % Assign any value different to the sentinel once done
        scn_dt( x_px, y_px ) = -sntnl_non_st ;
        % Check for flips
        % Flip x axis
          if scn_dt( 2 * cntr_1 - x_px, y_px ) == sntnl_non_st
          scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( 2 * cntr_1 - x_px - cntr_1 + n_non_st_1 + 1, y_px - cntr_2 + n_non_st_2 + 1 ) * flipud( psf_non_st_rtd ) ;
          % Assign any value different to the sentinel once done
          scn_dt( 2 * cntr_1 - x_px, y_px ) = -sntnl_non_st ;
          end
        % Flip y axis
          if scn_dt( x_px, 2 * cntr_2 - y_px ) == sntnl_non_st
          scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( x_px - cntr_1 + n_non_st_1 + 1, 2 * cntr_2 - y_px - cntr_2 + n_non_st_2 + 1 ) * fliplr( psf_non_st_rtd ) ;
          % Assign non-zero value
          scn_dt( x_px, 2 * cntr_2 - y_px ) = -sntnl_non_st ;
          end
        % Flip x and y axis
          if scn_dt( 2 * cntr_1 - x_px, 2 * cntr_2 - y_px ) == sntnl_non_st
          scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( 2 * cntr_1 - x_px - cntr_1 + n_non_st_1 + 1, 2 * cntr_2 - y_px - cntr_2 + n_non_st_2 + 1 ) * flipud( fliplr( psf_non_st_rtd ) ) ;
          % Assign any value different to the sentinel once done
          scn_dt( 2 * cntr_1 - x_px, 2 * cntr_2 - y_px ) = -sntnl_non_st ;
          end
        end % first quadrant (not including the center of the scene)
      end % i_px

    % Remaining pixels: repeat the loop
    px_non_st = find( scn_dt == sntnl_non_st ) ;
    n_non_st = numel( px_non_st ) ;
      for i_px = 1 : n_non_st
      px_tmp = px_non_st( i_px ) ;
      x_px = mod( px_tmp, sz_1 ) ;
      y_px = ceil( px_tmp / sz_1 ) ;
      idx_psf = 1 + round( sqrt( ( x_px - cntr_1 )^2 + ( y_px - cntr_2 )^2 ) * opt_img.px_scn_mas ) ;
      % Taking into account the spacing of the PSF basis
      idx_psf = round( idx_psf / opt_img.psf_spacing_mas ) ;
      psf_tmp = squeeze( psf_arry( idx_psf, :, : ) ) ;
      r_px = sign( x_px ) * round( sqrt( ( x_px - cntr_1 )^2 + ( y_px - cntr_2 )^2 ) ) ;
      % Reset its value to avoid some cross-talk with how other pixels filled out the PSF array
      psf_non_st = 0 * psf_non_st ;
      psf_non_st( cntr_non_st_1 - side_psf_st_1 : cntr_non_st_1 + side_psf_st_1, cntr_non_st_2 + r_px - side_psf_st_2 : cntr_non_st_2 + r_px + side_psf_st_2 ) = psf_tmp ;
      psf_non_st_rtd = imrotate( psf_non_st, -180 / pi * atan2( x_px - cntr_1, y_px - cntr_1 ), 'crop', 'bilinear' ) ;
      scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( x_px - cntr_1 + n_non_st_1 + 1, y_px - cntr_2 + n_non_st_2 + 1 ) * psf_non_st_rtd ;
      end % i_px
    toc
    % End of the convolution for the spinning starshade
    % Non-spinning case
    else
    % Use to track how much pixels (in %) have been convolved. First message when:
    i_prcnt = 10 ;
      for i_px = 1 : n_non_st
      px_tmp = px_non_st( i_px ) ;
      % DEC
      x_mas = round( ( mod( px_tmp, sz_1 ) - cntr_1 ) * opt_img.px_scn_mas ) ;
      % RA
      y_mas = round( ( ceil( px_tmp / sz_1 ) - cntr_2 ) * opt_img.px_scn_mas ) ;
      % Notice the order DEC/RA when bulding the PSF (checked that if the PSF is not re-centered -see mft_shift when re-centering xpp, ypp- the PSF centroid is located at (x_mas,y_mas), which agrees with the fact that (,x,y)=(DEC,RA) on the scene).
      psf_tmp = get_normalized_psf( i_bnd, opt_img, x_mas, y_mas ) ;
      x_px = mod( px_tmp, sz_1 ) - cntr_1 ;
      y_px = ceil( px_tmp / sz_1 ) - cntr_2 ;
      % Shifting the PSF to the location of the pixel under consideration (PS: it's not convenient to avoid this step by turning off the re-centering of the PSF in mft_shift, because one would lose some area around the PSF instead of obtaining a symmetric n_lambda_over_d area around its peak).
      % Reset its value to avoid some cross-talk with how other pixels filled out the PSF array
      psf_non_st = psf_non_st * 0 ;
      % psf_non_st has a big enough size to include the PSF. PS: notice that psf_tmp is already centered at the center of the array
      psf_non_st( cntr_non_st_1 + x_px - side_psf_st_1 : cntr_non_st_1 + x_px + side_psf_st_1, cntr_non_st_2 + y_px - side_psf_st_2 : cntr_non_st_2 + y_px + side_psf_st_2 ) = psf_tmp ;
      % Co-adding the result of the non-stationary convolution
      % Coordinates of the pixel on the original scene
      scn_non_st_cnv = scn_non_st_cnv + scn_non_st_dt( x_px + n_non_st_1 + 1, y_px + n_non_st_2 + 1 ) * psf_non_st ;
      % Tracking how much is computed
        if ( i_px == round( n_non_st / 100 * i_prcnt ) )
        disp( sprintf( '%i%% Done', i_prcnt ) )
        toc
        % Steps of 10%
        i_prcnt = 10 * round( ( i_prcnt + 10 ) / 10 ) ;
        end
      end % i_px
    end % End of the non-spinning Starshade convolution
end % if ( 0 )
  % Coadding the non-stationary results with the stationary ones
  scn_dt_cnv( cntr_1 - n_non_st_1 : cntr_1 + n_non_st_1, cntr_2 - n_non_st_2 : cntr_2 + n_non_st_2 ) = scn_dt_cnv( cntr_1 - n_non_st_1 : cntr_1 + n_non_st_1, cntr_2 - n_non_st_2 : cntr_2 + n_non_st_2 ) + scn_non_st_cnv ;

    % In case the FOV was too small, bring it back to its original size
    if ( l_1_1 + l_1_2 + l_2_1 + l_2_2 )
    scn_dt_cnv = scn_dt_cnv( l_1_1 + 1 : sz_1 - l_1_2, l_2_1 + 1 : sz_2 - l_2_2 ) ;
    % Likewise with psf_pd (used later on in the solar glint code)
      if exist( 'psf_pd', 'var' )
      psf_pd = psf_pd( l_1_1 + 1 : sz_1 - l_1_2, l_2_1 + 1 : sz_2 - l_2_2 ) ;
      end
    end

    % * Add new locus and subtract unperturbed one.
    if ( opt_img.lcs.do )
      if strcmp( opt_img.starshade.mode, 'spinning' )
      psf_0 = squeeze( psf_arry( 1, :, : ) ) ;
      else
      psf_0 = get_normalized_psf( i_bnd, opt_img, 0, 0 ) ;
      end
    scn_dt_cnv = add_new_locus( scn_dt_cnv, psf_0, opt_img, i_lmbd, psf_st ) ;
    end

    % * Add solar glint (get the PSF array if it is not present)
    if ( opt_img.slr_glnt.do )
    tic
    disp( 'Adding the solar glint' )
    scn_dt_cnv = add_solar_glint( scn_dt_cnv, opt_img ) ;
    toc
    end

    % * Add stars (get the PSF array if it is not present)

    % Write out the results as data and as an image
    % Cutout the last column/row if the input scene had an even number of pixels on either each of them (sz_1_odd, sz_2_odd)
    if ~( sz_1_odd )
    scn_dt_cnv = scn_dt_cnv( 1 : end - 1, : ) ;
    end
    if ~( sz_2_odd )
    scn_dt_cnv = scn_dt_cnv( :, 1 : end - 1 ) ;
    end

  % To avoid conflicts if files are reloaded
  opt_img_sv = opt_img ;
    % Storing the result for each wavelength
    if ( opt_img.sv_out_2 )
    save( fl_dt_out, 'scn_dt_cnv', 'opt_img_sv' ) ;
    disp( sprintf( 'Data stored: %s', fl_dt_out ) )
    end
  % If do_cnv_f = 0 then read the file (useful for plotting purposes)
  else % do_cnv_f. Next is done when convolution of diffuse areas is skipped, e.g., while planets are moved and there's no change on the backgropund
  load( fl_dt_out )
  % Do not set r_st_mas to the opt_img_sv value by default since it may be modified at the beginning of this code with r_stat_mas_new
  % This ensures that for every wavelength the PSF gets properly normalized in the next calls
  opt_img.psf_st_nrm = 1 ;
    if strcmp( opt_img.starshade.mode, 'spinning' )
    [ psf_arry opt_img ] = get_normalized_spinning_psf_array( i_bnd, opt_img ) ;
    end
    if strcmp( opt_img.starshade.mode, 'non-spinning' )
    [ psf_st opt_img ] = get_normalized_psf( i_bnd, opt_img, 0, opt_img.r_st_mas, 1 ) ;
    end
  end % do_cnv

% * Add planets (get the PSF array if it is not present) Since PSF size >> pixel, does not really matter the zodi issue, unless it would be >> planet emission
  if ( opt_img.plnt.add.do )
  % Keeping a copy of the scene without planets
  scn_dt_cnv_no_plnt_ph_s = scn_dt_cnv ; % It will be correctly defined with W/m2/um or Jy->photons/sec momentarily (saving one duplicate 2d array)
    if strcmp( opt_img.starshade.mode, 'non-spinning' )
    scn_dt_cnv = add_planets_non_spinning( scn_dt_cnv, psf_st, opt_img, i_lmbd, i_bnd ) ;
    else
    [ scn_dt_cnv opt_img ] = add_planets_spinning( scn_dt_cnv, psf_arry, opt_img, i_lmbd ) ;
    end
  else
  scn_dt_cnv_no_plnt_ph_s = 0 ;
  end

%%%%%%%%%%%%
% Converting into photons per second at the detector 
  if strcmp( lower( opt_img.scn.unts ), 'jy' )
  scn_dt_cnv_ph_s = scn_dt_cnv * opt_img.conv_Jy_to_photons_s_detector( i_lmbd ) ;
  scn_dt_cnv_no_plnt_ph_s = scn_dt_cnv_no_plnt_ph_s * opt_img.conv_Jy_to_photons_s_detector( i_lmbd ) ; % Recall the note above about re-using scn_dt_cnv_no_plnt_ph_s
  end

  if strcmp( lower( opt_img.scn.unts ), 'w/m2/um' )
  scn_dt_cnv_ph_s = scn_dt_cnv * opt_img.conv_W_m2_um_to_photons_s_detector( i_lmbd ) ;
  scn_dt_cnv_no_plnt_ph_s = scn_dt_cnv_no_plnt_ph_s * opt_img.conv_W_m2_um_to_photons_s_detector( i_lmbd ) ; % Recall the note above about re-using scn_dt_cnv_no_plnt_ph_s
  end

% Adding some jitter motion of the telescope, if specified. Modeled as a 2-D Gaussian with sigma=RMS
  if ( opt_img.jttr.do )
  disp( sprintf( 'Applying a jitter motion of RMS %2.2f mas', opt_img.jttr.rms_mas ) )
  scn_dt_cnv_ph_s = imgaussfilt( scn_dt_cnv_ph_s, opt_img.jttr.rms_mas / opt_img.px_scn_mas ) ;
  scn_dt_cnv_no_plnt_ph_s = imgaussfilt( scn_dt_cnv_no_plnt_ph_s, opt_img.jttr.rms_mas / opt_img.px_scn_mas ) ;
  end

% Converting to camera pixels if necessary

% Reversing the scene shift if it was applied 
  if ~isnan( opt_img.scn.ra_shft_px )
  scn_dt_cnv_ph_s = fraccircshift( scn_dt_cnv_ph_s, [ -opt_img.scn.dc_shft_px -opt_img.scn.ra_shft_px ] ) ;
  scn_dt_cnv_no_plnt_ph_s = fraccircshift( scn_dt_cnv_no_plnt_ph_s, [ -opt_img.scn.dc_shft_px -opt_img.scn.ra_shft_px ] ) ; 
    if ( i_lmbd == numel( opt_img.lmbd_arry_img_nm ) )
    disp( 'Lateral shift on the scene reversed.' )
    end
  end

  if ( opt_img.px_scn_mas ~= opt_img.px_cmr_mas )
  % The pixels that will cover the FOV on the camera. PS: forcing the output image to have an odd number of pixels
  opt_img.px_cmr_fov = 2 * ceil( size( scn_dt_cnv, 1 ) * opt_img.px_scn_mas / opt_img.px_cmr_mas / 2 ) + 1 ; % PS: stored as part of the options to parse into the cube construction for the PSF around each planet
  scn_dt_cnv_ph_s = imresize( scn_dt_cnv_ph_s, [ opt_img.px_cmr_fov opt_img.px_cmr_fov ], 'bilinear' ) ;
  % imresize needs to be rescaled to preserve the total number of counts: bilinear averages
  scn_dt_cnv_ph_s = scn_dt_cnv_ph_s * ( size( scn_dt_cnv, 1 ) / opt_img.px_cmr_fov )^2 ;
  % Data without planets
    if sum( sum( scn_dt_cnv_no_plnt_ph_s ) )
    scn_dt_cnv_no_plnt_ph_s = imresize( scn_dt_cnv_no_plnt_ph_s, [ opt_img.px_cmr_fov opt_img.px_cmr_fov ], 'bilinear' ) * ( size( scn_dt_cnv, 1 ) / opt_img.px_cmr_fov )^2 ;
    end
  end

% Data cubes
  if ( opt_img.cb.do )
  n_bck = numel( opt_img.cb.bckgrnd ) ;
    for i_bck = 1 : n_bck
    dt_tmp = scn_dt_cnv_ph_s + ( opt_img.cb.bckgrnd( i_bck ) - 1 ) * scn_dt_cnv_no_plnt_ph_s ;
    % Saving some memory. In most applications, 7 digit precision is enough
      if ( opt_img.sngl_prcsn )
      dt_tmp = single( dt_tmp ) ;
      end
    % For the very first call, erase any previous work
    opt_img.cb.strt = 0 ;
      if ( i_lmbd ==1 ) && ( i_bck == 1 )
      opt_img.cb.strt = 1 ;
      end
    fl_nm_cb = write_starshade_cube_file( dt_tmp, i_lmbd, opt_img, i_bck ) ;
    opt_img.cb.fl_nm{ i_bck } = fl_nm_cb ;
    end
  end

%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
% If convolution has to be computed/repeated
function do_cnv = do_cnv_f( fl_dt_out, opt_img )
do_cnv = 1 ; % By default, do the convolution
  if exist( fl_dt_out, 'file' ) == 2
    if ( opt_img.redo_2 )
    delete( fl_dt_out ) ;
    else
    do_cnv = 0 ;
    disp( sprintf( 'File %s exists. Not convolving the scene again', fl_dt_out ) )
    end
  end

%
% Conversion from spectral flux density (aka spectral irradiance) to ph/s (including other intermediate unit conversions)
% IT INCLUDES QE, MIRROR AREA, AND OPTICAL THROUGHPUT.
function opt = spectral_flux_density_to_photons_s_f( opt, i_lmbd, i_bnd )

% Delta_lambda (useful for unit conversion from Jy to photons)
      if i_lmbd == 1 % There's no before
        % At least some bandwidth
        if numel( opt.lmbd_arry_img_nm ) == 1
        opt.delta_lmbd_nm = opt.lmbd_img_2_nm( i_bnd ) - opt.lmbd_arry_img_nm( i_lmbd ) ;
          if ~( opt.delta_lmbd_nm )
          disp( 'WARNING: Both the lower and upper limits of the band have the same wavelength. Set different values. Returning.' )
          scn_dt_cnv_ph_s = 0 ;
          scn_dt_cnv_no_plnt_ph_s = 0 ;
          return
          end
        else
        opt.delta_lmbd_nm = opt.lmbd_arry_img_nm( i_lmbd + 1 ) - opt.lmbd_arry_img_nm( i_lmbd ) ;
        end
      end
      if ( i_lmbd > 1 ) && ( i_lmbd < numel( opt.lmbd_arry_img_nm ) ) 
      opt.delta_lmbd_nm = opt.lmbd_arry_img_nm( i_lmbd ) - opt.lmbd_arry_img_nm( i_lmbd - 1 ) ;
      end
      if ( i_lmbd == numel( opt.lmbd_arry_img_nm ) ) % The next wavelngth is the end of the band (see n_lmbd in sister_imaging_band.m)
        if numel( opt.lmbd_arry_img_nm ) == 1
        opt.delta_lmbd_nm = opt.lmbd_img_2_nm( i_bnd ) - opt.lmbd_arry_img_nm( i_lmbd ) ;
        else
        opt.delta_lmbd_nm = opt.lmbd_arry_img_nm( i_lmbd  ) - opt.lmbd_arry_img_nm( i_lmbd - 1 ) ;
        end
      end

% Jy are units of flux density in frequency per unit area (SI:1e-26 W/m^2/Hz): F_nu. First thing is to convert it into
% a differential of flux: dF = F_nu dnu. We work with wavelength and dnu=-c/lambda^2*dlambda. Obviating the minus sign (because we will integrate from shorter to longer wavelength, which is the reverse of lower to higher frequency):
% dF = c/lambda^2*F_nu*dlambda. This will have units of power per unit area.
% We are not having differential steps of wavelength, but finite sums. The smaller the value of delta lambda scene, the more accurate. 
int_lmbd_nm = opt.lmbd_tmp_nm : opt.lmbd_tmp_nm + opt.delta_lmbd_nm ;
lmbd_tmp_eq_m = int_lmbd_nm * opt.units_wavelength_si ; % m
% PS: physconst('LightSpeed')=299792458 m
% NB: notice this will be for the delta lambda in use
opt.conv_Jy_to_W_m2_arry = 1e-26 * 299792458 ./ lmbd_tmp_eq_m.^2 * opt.delta_lmbd_nm * opt.units_wavelength_si ;
% Deriving the mean value for the conversion between both units of spectral flux density
opt.conv_Jy_to_W_m2_um( i_lmbd ) = mean( opt.conv_Jy_to_W_m2_arry ) / opt.delta_lmbd_nm / opt.units_wavelength_si / 1e6 ;
% Conversion from Jy to photons/s/m2 (h=6.62607004e-34). Taking the mean value. That is accurate if the flux density in Jy is constant across the wavelength bin width. 
% NB: notice this will be for the delta lambda in use (divided by delta lambda will give Jy->photons/s/m^2/nm)
opt.conv_Jy_to_photons_s_m2( i_lmbd ) = mean( opt.conv_Jy_to_W_m2_arry / ( 6.62607004e-34 * 299792458 ) .* lmbd_tmp_eq_m ) ;

% Conversion from W/m2/um to photons/s/m2 (notice this will be for the delta lambda in use)
opt.conv_W_m2_um_to_photons_s_m2( i_lmbd ) = 1e6 * mean( lmbd_tmp_eq_m / ( 6.62607004e-34 * 299792458 ) ) * opt.delta_lmbd_nm * opt.units_wavelength_si ;

% Collecting area, QE and sky transmittance (ground telescope)
% Telescope area
  if isfield( opt, 'secondary_size' )
  opt.ratio_obscured_area = opt.secondary_size^2 ;
  end
  if strfind( opt.starshade.nominal_filename, 'NI2' )
  % From the WFIRST-S pupil (m2=fitsread('input_scenes/locus/in/pupil_D1Kpix_2048.fits'); ( ( numel( find( m2(:) <= 0.5 ) ) - ( 2048*2048 - pi * 1001^2 / 4 ) ) ) / pi / (1001*1001) * 4 = 0.174104209179090. PS: with ==0, it gives 0.1654 )
  opt.ratio_obscured_area = 0.1741 ; % That is, 0.1024=0.32^2, and 0.0717 due to the struts of the secondary.
    if ~isfield( opt, 'qe' )
    % Data from Patrick Morrissey (email 02/05/18)
    fl_qe_ccd201_20 = sprintf( '%sinput_scenes/EMCCD201_20_QE.txt', sister_installation_path() ) ;
      if ( exist( fl_qe_ccd201_20, 'file' ) == 2 )
      qe_tmp = load( fl_qe_ccd201_20 ) ;
      lmbd_qe_nm = qe_tmp( :, 1 ) * 1e3 ;
      qe_ccd201_20 = qe_tmp( :, 2 ) ;
      clear qe_tmp 
      else % Backwards compatibility (before 04/09/20)
      lmbd_qe_nm = [ 400, 500, 600, 700, 800, 900, 1000 ] ; % nm
      qe_ccd201_20 = [ 41, 81, 86, 73, 47, 18, 2 ] / 100 ;
      % Storing the data in an ASCII file
      qe_tmp = [ lmbd_qe_nm / 1e3 ; qe_ccd201_20 ] ;
      fl_id = fopen( fl_qe_ccd201_20, 'w' ) ;
      fprintf( fl_id, '%1.3f %1.3f\n', qe_tmp ) ;
      fclose( fl_id ) ;
      clear qe_tmp
      end
    opt.qe_lambda( i_lmbd ) = mean( interp1( lmbd_qe_nm, qe_ccd201_20, lmbd_tmp_eq_m * 1e9, 'spline' ) ) * opt.cbe_qe_fct ;
    end
  % PS:  Looking at Andrew Romero-Wolf's compilation of effects (05/29/18 Starshade Rendezvous Mission Study), this value can be 0.29 or 0.17 if some throughput factors related with the modified Teledyne-e2v CCD201-20 electron multiplying CCD (EMCCD) are included. In sister.m the default value for any occulter with NI2 is set to 0.23 (the mean value as a compromise between known and unknown).
  end % NI2

% For any other occulter (assuming constant optical throughput and QE, unless there's a file with wavelength dependent values, see the example of EMCCD201_20)
  if ~isfield( opt, 'optical_throughput' )
    if isfield( opt, 'optical_throughput_lambda' )
    disp( sprintf( 'Optical throughput set to %1.3f', opt.optical_throughput_lambda( i_lmbd ) ) )
    else
    disp( 'Set opt.optical_throughput (Obscuration, QE, and PSF core are accounted for separately by the simulation itself or other parameters of the simulations). Stopped' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  else
    if ~ischar( opt.optical_throughput )
    opt.optical_throughput_lambda( i_lmbd ) = opt.optical_throughput * opt.cbe_optical_throughput_fct ;
    else
    % In case a file is provided 
    throughput = load( opt.optical_throughput ) ; % Assuming wavelength is in the first column
    % Guessing whether it is nm or micron (throughput(:,1) is the wavelength)
      if ( mean( throughput( :, 1 ) ) > 100 ) && ( mean( throughput( :, 1 ) ) < 10000 )
      fct_wl = 1 ; % nm
      end
      if ( mean( throughput( :, 1 ) ) > 1e-3 ) && ( mean( throughput( :, 1 ) ) < 10 )
      fct_wl = 1e3 ; % nm
      end
    % Check
      if ( opt.lmbd_tmp_nm < min( throughput( :, 1 ) * fct_wl ) )
      disp( sprintf( 'The minimum wavelength available in the optical throughput file %s is %.1f nm, whereas the imaging simulation is for %.2f. Please set the range of the imaging simulation within the range of the optical throughput file, or modify the file. Stopped', opt.optical_throughput, min( throughput( :, 1 ) * fct_wl ), opt.lmbd_tmp_nm ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    opt.optical_throughput_lambda( i_lmbd ) = mean( interp1( throughput( :, 1 ) * fct_wl, throughput( :, 2 ),  lmbd_tmp_eq_m * 1e9, 'spline' ) ) * opt.cbe_optical_throughput_fct ;
    end
  end

  if ( opt.cbe_optical_throughput_fct == 1 )
  disp( sprintf( 'Optical throughput set to %1.3f', opt.optical_throughput_lambda( i_lmbd ) ) )
  else
  disp( sprintf( 'Optical throughput set to %1.3f (CBE factor of %1.3f included)', opt.optical_throughput_lambda( i_lmbd ), opt.cbe_optical_throughput_fct ) )
  end

 % For now assuming constant QE, unless there's a prescription like with the modified Teledyne-e2v CCD201-20 electron multiplying CCD (EMCCD) above.
  if ~isfield( opt, 'qe' ) 
    if ~isfield( opt, 'qe_lambda' )
    disp( 'WARNING: please set opt.qe in the configuration file, or set a rule in Jy_to_photons_s_f.m in convolve_with_one_wavelength.m. Stopped' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  else
    if ~ischar( opt.qe )
    opt.qe_lambda( i_lmbd ) = opt.qe * opt.cbe_qe_fct ;
    else
    % In case a file is provided
    qe_cmr = load( opt.qe ) ; % Assuming wavelength in microns and first column
    % Guessing whether it is nm or micron (qe_cmr(:,1) is the wavelength)
      if ( mean( qe_cmr( :, 1 ) ) > 100 ) && ( mean( qe_cmr( :, 1 ) ) < 10000 )
      fct_wl = 1 ; % nm
      end
      if ( mean( qe_cmr( :, 1 ) ) > 1e-3 ) && ( mean( qe_cmr( :, 1 ) ) < 10 )
      fct_wl = 1e3 ; % micrometer
      end
    % Check
      if ( opt.lmbd_tmp_nm < min( qe_cmr( :, 1 ) * fct_wl ) )
      disp( sprintf( 'The minimum wavelength available in the QE file %s is %.1f nm, whereas the imaging simulation is for %.2f. Please set the range of the imaging simulation within the range of the QE file, or modify the file. Stopped', opt.qe, min( qe_cmr( :, 1 ) * fct_wl ), opt.lmbd_tmp_nm ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    opt.qe_lambda( i_lmbd ) = mean( interp1( qe_cmr( :, 1 ) * fct_wl, qe_cmr( :, 2 ), lmbd_tmp_eq_m * 1e9, 'spline' ) ) * opt.cbe_qe_fct ;
    end
  end

  if ( opt.cbe_qe_fct == 1 )
  disp( sprintf( 'Camera QE set to %1.3f', opt.qe_lambda( i_lmbd ) ) )
  else
  disp( sprintf( 'Camera QE set to %1.3f (CBE factor of %1.3f included)', opt.qe_lambda( i_lmbd ), opt.cbe_qe_fct ) )
  end

% Check
  if ~isfield( opt, 'ratio_obscured_area' )
  disp( 'Set opt.ratio_obscured_area. It is the ratio between the areas of the primary mirror and the area blocked by the secondary and struts, if any. If opt.pupil_filename is ''ideal'', one can set opt.secondary_size instead. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% In the case of ground telescopes, the atmosphere absorps some of the emissions
% The files were produced by Stefan Kimensberger (sic) "the ESO sky model is a product of a theoretical calculation cross calibrated with thousands of spectra. The files here are based on online version release 2.0.6, based on the papers Noll et al. (2012, A&A 543, A92) and Jones et al. (2013, A&A 560, A91). This version does not yet include the most recent minor revision of theAerosol extinction curve by Jones et al 2019.
% See also: https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  if ( opt.grnd_tlscp.do )
  % Reading the data (using the method from https://www.mathworks.com/matlabcentral/answers/452816-read-inconsistent-ascii-file-to-matrix, Jan)
  sky_tmp_1 = fileread( sprintf( '%ssky_transmittance_%s_moon_z45.dat', opt.grnd_tlscp.dr, opt.grnd_tlscp.mn_phs ) ) ;
  sky_tmp_2 = strsplit( sky_tmp_1, char( 10 ) ) ;
  i_tmp_2 = 1 ;
    for i_tmp = 1 : numel( sky_tmp_2 )
      if ~isempty( sky_tmp_2{ i_tmp } ) && any( sky_tmp_2{ i_tmp }( 1 ) == '1234567890-.' )
      sky_tmp_3 = sscanf( sky_tmp_2{ i_tmp }, '%g %g' ) ;
      sky_lmbd_nm_arry( i_tmp_2 ) = sky_tmp_3( 1 ) ;
      sky_trnsmttnc_arry( i_tmp_2 ) = sky_tmp_3( 2 ) ;
      i_tmp_2 = i_tmp_2 + 1 ;
      end
    end
  % Integrating across the wavelength slice (recall the data are per micrometer)
  int_lmbd_nm = opt.lmbd_tmp_nm : opt.lmbd_tmp_nm + opt.delta_lmbd_nm ;
  % Check of consistency
    if ( int_lmbd_nm( 1 ) < sky_lmbd_nm_arry( 1 ) ) || ( int_lmbd_nm( end ) > sky_lmbd_nm_arry( end ) )
    disp( 'The imaging band is not covered by the sky model data. Check the wavelength range. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  % Clipping with the imaging band (in case some user sets opt.delta_lambda_imaging_nm too large during some tests)
  int_lmbd_nm = int_lmbd_nm( int_lmbd_nm >= opt.lmbd_img_1_nm & int_lmbd_nm <= opt.lmbd_img_2_nm ) ;
  opt.grnd_tlscp.sky_trnsmttnc( i_lmbd ) = mean( interp1( sky_lmbd_nm_arry, sky_trnsmttnc_arry, int_lmbd_nm, 'spline' ) ) ;
  else
  opt.grnd_tlscp.sky_trnsmttnc( i_lmbd ) = 1 ;
  end

opt.optical_clear_area_m2 = ( 1 - opt.ratio_obscured_area ) * pi * opt.dmtr_tlscp_m^2 / 4 ;
% This conversion takes into account the collecting area:
opt.conv_Jy_to_photons_s_telescope( i_lmbd ) = opt.conv_Jy_to_photons_s_m2( i_lmbd ) * opt.optical_clear_area_m2 ;
% Useful when dealing with units of ph/s/nm/arcsec^2
opt.conv_Jy_to_photons_s_m2_nm( i_lmbd ) = opt.conv_Jy_to_photons_s_m2( i_lmbd ) / opt.delta_lmbd_nm ;
% Similar conversions for W/m2/um
opt.conv_W_m2_um_to_photons_s_telescope( i_lmbd ) = opt.conv_W_m2_um_to_photons_s_m2( i_lmbd ) * opt.optical_clear_area_m2 ;
opt.conv_W_m2_um_to_photons_s_m2_nm( i_lmbd ) = opt.conv_W_m2_um_to_photons_s_m2( i_lmbd ) / opt.delta_lmbd_nm ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual conversion from Jy to photons at the detector (Recall opt_img.grnd_tlscp.sky_trnsmttnc=1; unless a ground telescope is used) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.conv_Jy_to_photons_s_detector( i_lmbd ) = opt.conv_Jy_to_photons_s_telescope( i_lmbd ) * opt.optical_throughput_lambda( i_lmbd ) * opt.qe_lambda( i_lmbd ) * opt.grnd_tlscp.sky_trnsmttnc( i_lmbd ) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual conversion from W/m2/um to photons at the detector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.conv_W_m2_um_to_photons_s_detector( i_lmbd ) = opt.conv_W_m2_um_to_photons_s_telescope( i_lmbd ) * opt.optical_throughput_lambda( i_lmbd ) * opt.qe_lambda( i_lmbd ) * opt.grnd_tlscp.sky_trnsmttnc( i_lmbd ) ;

% Getting the scene data for a given wavelength
function [ scn_dt opt_img ] = get_scene_data( opt_img, i_lmbd )
% Sentinel
scn_dt = NaN ;

% Different cases: from file or created
  if ~( opt_img.scn.do )
    if numel( strfind( opt_img.scn.nm, 'modern_cube_zodi1inc' ) )
    scn_dt = read_haystacks_project( opt_img ) ;
    % Checking the scene data have an odd number of pixels on each dimension
      if size( scn_dt, 1 ) ~= 2 * floor( size( scn_dt, 1 ) / 2 ) + 1
      disp( 'External data for the scene must be an array or a cube array with an odd number of elements. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
      if ~isfield( opt_img.scn, 'fov_diam_mas' )
      opt_img.scn.fov_diam_mas = size( scn_dt, 1 ) * opt_img.px_scn_mas ;
      else
        if ( opt_img.scn.fov_diam_mas > size( scn_dt, 1 ) * opt_img.px_scn_mas )
        disp( sprintf( 'WARNING: the external data have a FOV of %ix%i mas, smaller than %ix%i, the one set in the configuration file. Continuing and setting the FOV to be the one in the external data. Check if this is fine, or modify your choices in the configuration file.', size( scn_dt, 1 ) * opt_img.px_scn_mas, size( scn_dt, 1 ) * opt_img.px_scn_mas, opt_img.scn.fov_diam_mas, opt_img.scn.fov_diam_mas ) )
        opt_img.scn.fov_diam_mas = size( scn_dt, 1 ) * opt_img.px_scn_mas ;
        else
        % Getting the FOV with data
        scn_cntr = ( size( scn_dt, 1 ) + 1 ) / 2 ;
        hlf_fov = floor( opt_img.scn.fov_diam_mas / opt_img.px_scn_mas / 2 ) ;
        scn_dt = scn_dt( scn_cntr - hlf_fov : scn_cntr + hlf_fov, scn_cntr - hlf_fov : scn_cntr + hlf_fov ) ;
        end
      end
    end
  end % Scene from file

% Getting the star flux
opt_img = get_star_flux( scn_dt, opt_img, i_lmbd ) ;

% Creating a scene from scratch
  if ( opt_img.scn.do )
  opt_img.n_px_scn = 2 * floor( opt_img.scn.fov_diam_mas / opt_img.px_scn_mas / 2 ) + 1 ;
  cntr_scn = ( opt_img.n_px_scn + 1 ) / 2 ;
  scn_dt = zeros( opt_img.n_px_scn, opt_img.n_px_scn );
  scn_dt( cntr_scn, cntr_scn ) = opt_img.str.flx( i_lmbd ) ;
  % Add exozodi (extragalactic and stars background is added later on)
    if ( opt_img.zd.fct )
    scn_dt = add_exozodi( scn_dt, opt_img, i_lmbd ) ;
    end
  % Add exo-Kuiper 
    if ( opt_img.kpr.do )
    scn_dt = scn_dt + add_exokuiper_belt( opt_img, i_lmbd ) ;
    end
  end

% Adding local zodiacal light (it follows the solar spectrum, STScI, TIR-CRDS-2015-01.pdf, R. Diaz, 04/01/2015. http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/reference-data-for-calibration-and-tools/documentation/_documents/TIR-CRDS-2015-01.pdf)
% In order to speed up the simulations with mutiple epochs, SISTER derives a template value for the convolution of the local zodiacal light. This is done by choosing the surface brightness to have a 0th magnitude and in sister_imaging_band.m is re-scaled to its actual value.
  if ( opt_img.lcl_zd.do ) && ( opt_img.lcl_zd.do_rf )
  % Reminder: the flux from the Sun at 0 mag is being estimated in get_star_flux.m below
  scn_dt = scn_dt + opt_img.str.flx_sun_0mg * ( opt_img.px_scn_mas / 1000 )^2 ;
  disp( 'Adding the local zodiacal light...' )
  end

% Adding the Earth's sky model for ground telescopes
% Basically, we include in SISTER three scenarios for the sky radiance and transmittance depending on the moon phase. Each of them, includes scattered Moonlight, local zodiacal light, scattered starlight, emisison lines of Upper Atmosphere, airglow (residual continuum), and Molecular Emission of Lower Atmosphere (the latter essentially zero below 1.7 micron though).
% The files were produced by Stefan Kimensberger (sic) "the ESO sky model is a product of a theoretical calculation cross calibrated with thousands of spectra. The files here are based on online version release 2.0.6, based on the papers Noll et al. (2012, A&A 543, A92) and Jones et al. (2013, A&A 560, A91). This version does not yet include the most recent minor revision of theAerosol extinction curve by Jones et al 2019.
% See also: https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  if ( opt_img.grnd_tlscp.do )
  % Reading the data (units: photons/second/m^2/micrometer/arcsec^2) (using the method from https://www.mathworks.com/matlabcentral/answers/452816-read-inconsistent-ascii-file-to-matrix, Jan)
  sky_tmp_1 = fileread( sprintf( '%ssky_radiance_%s_moon_z45.dat', opt_img.grnd_tlscp.dr, opt_img.grnd_tlscp.mn_phs ) ) ;
  sky_tmp_2 = strsplit( sky_tmp_1, char( 10 ) ) ;
  i_tmp_2 = 1 ;
    for i_tmp = 1 : numel( sky_tmp_2 )
      if ~isempty( sky_tmp_2{ i_tmp } ) && any( sky_tmp_2{ i_tmp }( 1 ) == '1234567890-.' )
      sky_tmp_3 = sscanf( sky_tmp_2{ i_tmp }, '%g %g' ) ;
      sky_lmbd_nm_arry( i_tmp_2 ) = sky_tmp_3( 1 ) ;
      sky_ph_s_m2_micrometer_arcsec2_arry( i_tmp_2 ) = sky_tmp_3( 2 ) ;
      i_tmp_2 = i_tmp_2 + 1 ;
      end
    end
  % Integrating across the wavelength slice (recall the data are per micrometer)
  int_lmbd_nm = opt_img.lmbd_tmp_nm : opt_img.lmbd_tmp_nm + opt_img.delta_lmbd_nm ;
  % Check of consistency
    if ( int_lmbd_nm( 1 ) < sky_lmbd_nm_arry( 1 ) ) || ( int_lmbd_nm( end ) > sky_lmbd_nm_arry( end ) )
    disp( 'The imaging band is not covered by the sky model data. Check the wavelength range. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end 
  % Clipping with the imaging band (in case some user sets opt.delta_lambda_imaging_nm too large during some tests)
  int_lmbd_nm = int_lmbd_nm( int_lmbd_nm >= opt_img.lmbd_img_1_nm & int_lmbd_nm <= opt_img.lmbd_img_2_nm ) ;
  % Mean value of the sky emission. Later on this value is multiplied by the width of the wavelength slice, removing the 1/nm units. Notice that the files provided have 1/micrometer and we change them to 1 /nm with the 1/1000 factor.
  opt_img.grnd_tlscp.sky_ph_s_m2_nm_arcsec2 = mean( interp1( sky_lmbd_nm_arry, sky_ph_s_m2_micrometer_arcsec2_arry, int_lmbd_nm, 'spline' ) ) / 1000 ;
  % Bringing it back to Jy or W/m2/nm
    if strcmp( lower( opt_img.scn.unts ), 'jy' )
    opt_img.grnd_tlscp.sky_Jy = opt_img.grnd_tlscp.sky_ph_s_m2_nm_arcsec2 * ( opt_img.px_scn_mas / 1000 )^2 / opt_img.conv_Jy_to_photons_s_m2_nm( i_lmbd ) ;
    scn_dt = scn_dt + opt_img.grnd_tlscp.sky_Jy ;
    end
    if strcmp( lower( opt_img.scn.unts ), 'w/m2/um' )
    opt_img.grnd_tlscp.sky_W_m2_um = opt_img.grnd_tlscp.sky_ph_s_m2_nm_arcsec2 * ( opt_img.px_scn_mas / 1000 )^2 / opt_img.conv_W_m2_um_to_photons_s_m2_nm( i_lmbd ) ;
    scn_dt = scn_dt + opt_img.grnd_tlscp.sky_W_m2_um ;
    end
  end

% Remove planets if required (very fast, no need to stored scenes without planets if they had them on them and one wants to erase them to simulate orbital motion)
  if ( opt_img.plnt.rmv.do )
  [ scn_dt  opt_img ] = remove_planets_from_scene( scn_dt, opt_img ) ;
  end

% Add extragalactic background (needs full convolution as the scene. Other effects which only require some pixels are deal with after the long convolution, because they can re-use this result)
  if ( opt_img.bckgrnd.add )
  scn_dt = add_extragalactic_background( scn_dt, opt_img, i_lmbd ) ;
  end

% Checking there was an assignment
  if isnan( scn_dt )
  disp( 'The scene data could not be loaded. Stopped.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  % Applying a lateral shift if there's any (only one may be set up: either mas or meters. Also if the shift on RA or DEC is set, both must be set, see sister_imaging.m)
  opt_img.scn.ra_shft_px = NaN ;
  opt_img.scn.dc_shft_px = NaN ;
  % Getting the equivalent displacement in meters and mas
    if ~isnan( opt_img.scn.ra_shft_mas )
    opt_img.scn.ra_shft_m = opt_img.scn.ra_shft_mas / 180 * pi / 3600e3 * opt_img.dst_strshd_tlscp_m ;
    opt_img.scn.dc_shft_m = opt_img.scn.dc_shft_mas / 180 * pi / 3600e3 * opt_img.dst_strshd_tlscp_m ;
    end
    if ~isnan( opt_img.scn.ra_shft_m )
    opt_img.scn.ra_shft_mas = opt_img.scn.ra_shft_m / opt_img.dst_strshd_tlscp_m * 180 / pi * 3600e3 ;
    opt_img.scn.dc_shft_mas = opt_img.scn.dc_shft_m / opt_img.dst_strshd_tlscp_m * 180 / pi * 3600e3 ;
    opt_img.scn.ra_shft_px = opt_img.scn.ra_shft_mas / opt_img.px_scn_mas ; % It might be a fraction of a pixel for extended objects. Planets are rounded to nearest pixel.
    opt_img.scn.dc_shft_px = opt_img.scn.dc_shft_mas / opt_img.px_scn_mas ;
    end

  % Shift the scene if needs be. The shift is circular. In principle, the shifts are going to be a few mas maximum and the scene at the edges of the scene are not used for analysis.
    if ~isnan( opt_img.scn.ra_shft_px )
    scn_dt = fraccircshift( scn_dt, [ opt_img.scn.dc_shft_px opt_img.scn.ra_shft_px ] ) ;
      if ( i_lmbd == 1 )
      disp( sprintf( 'Scene laterally displaced by %.1f mas, %.1f meters, %i pixels in the horizontal direction and %.1f mas, %.1f meters, %i pixels in the vertical direction', opt_img.scn.ra_shft_mas, opt_img.scn.ra_shft_m, opt_img.scn.ra_shft_px, opt_img.scn.dc_shft_mas, opt_img.scn.dc_shft_m, opt_img.scn.dc_shft_px ) )
      end
    end

% Getting the flux of the star
function opt = get_star_flux( scn_dt, opt, i_lmbd )
% It brings star_flux_0mag
  if strcmp( lower( opt.scn.unts ), 'jy' ) 
  load( [ opt.scn_dr 'HST_stars_0mag_flux_jy.mat' ] ) ;
  end
  if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
  load( [ opt.scn_dr 'HST_stars_0mag_flux_w_m2_um.mat' ] ) ;
  end
% Deriving an interpolated/extrapolated star flux for all types, but the Sun, and F type (only f5v is available in SISTER). Recall that the available spectra in the code are: 'a0v', 'a5v', 'f5v', 'g0v', 'g5v' (closer to sun), 'k0v', 'k5v', 'm0v', 'm5v'
  % Sun's spectrum is dealt separately
  if ~strcmp( opt.str.tp, 'sun' ) 
    if ~strcmp( opt.str.tp( 1 ), 'f' )
    str_tp_lw_flx = eval( [ 'star_flux_0mag.' opt.str.tp( 1 ) '0v' ] ) ;
    str_tp_up_flx = eval( [ 'star_flux_0mag.' opt.str.tp( 1 ) '5v' ] ) ;
    % Stellar spectral types range from 0 to 9
    str_flx_arry = str_tp_lw_flx + str2num( opt.str.tp( 2 ) ) * ( str_tp_up_flx - str_tp_lw_flx ) / 5 ;
    else
    opt.str.tp = 'f5v' ;
    str_flx_arry = eval( [ 'star_flux_0mag.' opt.str.tp ] ) ;
    end
  end  % Interpolate for all types except the Sun and F stars

% Sun's apparent spectral irradiance 
% Default case: C.N.A. Willmer, Astrophysical Journal Supplements, 2018, 236, 47, or http://mips.as.arizona.edu/~cnaw/sun.html
sun_W_m2_micron = load( [ opt.scn_dr 'sun_composite_spectrum_willmer_2018_W_m2_micron.txt' ] ) ;

% In case the AM0 (2000) is used (https://rredc.nrel.gov/solar//spectra/am0/ASTM2000.html (updated 2014))
  if ( opt.str.am0 )
  sun_W_m2_micron = load( [ opt.scn_dr 'AM0_2000_solar_spectrum_W_m2_micron.txt' ] ) ;
  disp( 'WARNING: loading the Sun''s irradiance from AM0 (2000) instead of Willmer (2018)' )
  end

% In case the WMO (1985) is used (https://www.nrel.gov/grid/solar-resource/spectra-wehrli.html)
  if ( opt.str.wmo )
  sun_W_m2_micron = load( [ opt.scn_dr 'wmo_solar_spectrum_W_m2_micron.txt' ] ) ;
  disp( 'WARNING: loading the Sun''s irradiance from WMO (1985) instead of Willmer (2018)' )
  end
% Bringing the irradiance from apparent to 0th V mag
sun_flx_0th_mag = sun_W_m2_micron( :, 2 ) * 10^(-26.76/2.5) ; % Absolute magnitude of the Sun from http://mips.as.arizona.edu/~cnaw/sun.html)
sun_lmbd_nm = sun_W_m2_micron( :, 1 ) * 1e3 ;
  if strcmp( lower( opt.scn.unts ), 'jy' )
  sun_flx_0th_mag = sun_flx_0th_mag * 1e26 / 299792458 .* (sun_lmbd_nm * 1e-9).^2 * 1e6 ; 
  end

  if strcmp( opt.str.tp, 'sun' )
  star_flux_0mag.lambda_nm = sun_lmbd_nm ;
  str_flx_arry = sun_flx_0th_mag ;
  end

% Getting the mean flux in the small wavelength interval (1 nm spacing)
int_lmbd_nm = opt.lmbd_tmp_nm : opt.lmbd_tmp_nm + opt.delta_lmbd_nm ;
% Clipping the range between the available wavelengths
int_lmbd_nm = int_lmbd_nm( int_lmbd_nm >= min( star_flux_0mag.lambda_nm ) & int_lmbd_nm <= max( star_flux_0mag.lambda_nm ) ) ;
% Clipping with the imaging band (in case some user sets opt,delta_lambda_imaging_nm too large during some tests)
int_lmbd_nm = int_lmbd_nm( int_lmbd_nm >= opt.lmbd_img_1_nm & int_lmbd_nm <= opt.lmbd_img_2_nm ) ;
  if ( opt.scn.do )
  % Mean value in the region of interest (linear interpolation because absorption lines are not smooth enough for higher order interpolation to be more accurate)
  opt.str.flx( i_lmbd ) = mean( interp1( star_flux_0mag.lambda_nm, str_flx_arry, int_lmbd_nm ) ) * 10^( -opt.str.mg / 2.5 ) ;
  % If the Haystacks scenes are used, we get the star flux from the scenes. It should be at the center of the scene.
  else
    if numel( findstr( 'modern_cube_zodi', opt.scn.nm ) ) 
    % Haystacks cubes are given in Jy
      if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
      disp( 'Haystacks data are provided in units of Jy, but the simulation is set to units of W/m^2/um. Please, check the setup is correct. For instance, set opt.scene.units=''Jy''.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    cntr_scn_1 = floor( ( size( scn_dt, 1 ) + 1 ) / 2 ) ;
    cntr_scn_2 = floor( ( size( scn_dt, 2 ) + 1 ) / 2 ) ;
    % Just in case the scenes had some even number of pixels
    l = 1 ;
    opt.str.flx( i_lmbd ) = max( max( scn_dt( cntr_scn_1 - l : cntr_scn_1 + l, cntr_scn_2 - l : cntr_scn_2 + l ) ) ) ;
    opt.str.mg = 4.81 ; % Absolute magnitude of the Sun (http://mips.as.arizona.edu/~cnaw/sun.html)
    end
  end

% It's also valuable to store the ratio with the Sun's flux at that wavelength (it is also used when deriving the local zodiacal light)
opt.str.flx_rt_sun = opt.str.flx( i_lmbd ) / ( sum( interp1( sun_lmbd_nm, sun_flx_0th_mag, int_lmbd_nm ) ) / numel( int_lmbd_nm ) * 10^( -opt.str.mg / 2.5 ) ) ;
% And also n *Jy* wrt 500 nm because the simple exozodi model uses a reference zodi computed at 500 nm
opt.str.flx_rt_sun_500nm_jy = opt.str.flx( i_lmbd ) / ( mean( interp1( sun_lmbd_nm, sun_flx_0th_mag, 470 : 530 ) ) * 10^( -opt.str.mg / 2.5 ) ) ;
% if the spectral irradiance is in W/m^2/micro-meter, translate to Jy (but for constant factors)
  if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
  opt.str.flx_rt_sun_500nm_jy = mean( interp1( star_flux_0mag.lambda_nm, str_flx_arry .* star_flux_0mag.lambda_nm.^2, int_lmbd_nm ) ) / mean( interp1( sun_lmbd_nm, sun_flx_0th_mag .* sun_lmbd_nm.^2, 470 : 530 ) ) ;
  end
% Finally, we also store the mean flux corresponding to the Sun in the wavelength slice (useful for estimating the local zodical light)
opt.str.flx_sun_0mg = mean( interp1( sun_lmbd_nm, sun_flx_0th_mag, int_lmbd_nm ) ) ;
% Flux for a A0V star with 0th magnitude. Useful to derive magnitudes in the same band.
a0v_flx = eval( [ 'star_flux_0mag.a0v' ] ) ;
% Wavelength array has to be re-loaded because if the star was the Sun, then the alternative files with the Sun's spectrum have been used.
% It brings star_flux_0mag
  if strcmp( lower( opt.scn.unts ), 'jy' )
  load( [ opt.scn_dr 'HST_stars_0mag_flux_jy.mat' ] ) ;
  end
  if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
  load( [ opt.scn_dr 'HST_stars_0mag_flux_w_m2_um.mat' ] ) ;
  end
opt.str.flx_a0v_0mg = mean( interp1( star_flux_0mag.lambda_nm, a0v_flx, int_lmbd_nm ) ) ;
% Flux of a A0V star within the V band (550 nm, with passband 0.16 *550 nm = 506-594 nm.
opt.str.flx_a0v_0mg_V_bnd = mean( interp1( star_flux_0mag.lambda_nm, a0v_flx, 506 : 594 ) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing the planets from a scene that came with planets
% Notice that if a file without planets would get the 'planet removal', it would only affect if there are point-like sources on the scene at those positions. Otehrwise, the result is very similar to what was before. So, it's not dramatic if some file gets through inadvertedly.
function [ scene opt ] = remove_planets_from_scene( scene, opt )
% Planet location: planet_pos_pix_1 and planet_pos_pix_2 (and add str.arc_ra, arc_dec):
% !!!
% These are defined such that scene( planet_pos_pix_1, planet_pos_pix_2 ) = planet (pixel if point-like or peak if extended).
% !!!

% Depends on each case

% Haystacks project case for different sub-cases
  if findstr( 'modern_cube_zodi', opt.scn.nm )
  % It brings planet_pos_pix_1 and _2
  load( [ opt.scn_dr opt.scn.nm '_planet_positions.mat' ] )
  n_plnt = numel( planet_pos_pix_1 ) ;
    if n_plnt ~= numel( planet_pos_pix_2 )
    disp( 'Number of coordinates for planet positions is different. Correct.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    % Warning message in case a file arrives here with the action of removing planets but there is no planet list assigned to it.
    if ~exist( 'planet_pos_pix_1', 'var' )
    disp( 'There should be two variables named planet_pos_pix_1 and planet_pos_pix_2. They are not there. Checked and correct.' )
    end
  % Removing the planets
  % It may happen that some error of +/- 1 might occur if some planet positions were determined in a language that starts with 0 or with 1
  l = 3 ;
    for i_plnt = 1 : n_plnt
    sqr_tmp = scene( planet_pos_pix_1( i_plnt ) - l : planet_pos_pix_1( i_plnt ) + l, planet_pos_pix_2( i_plnt ) - l : planet_pos_pix_2( i_plnt ) + l ) ;
    q_mx = find( sqr_tmp == max( max( sqr_tmp ) ) ) ;
    planet_pos_1 = planet_pos_pix_1( i_plnt ) - l - 1 + ceil( q_mx / ( 2 * l + 1 ) ) ;
    a( i_plnt ) = planet_pos_1 ;
    planet_pos_2 = planet_pos_pix_2( i_plnt ) - l - 1 + mod( q_mx, ( 2 * l + 1 ) ) ;
    b( i_plnt ) = planet_pos_2 ;
    l2 = 1 ;
    % Background estimation
    bckgrnd_tmp = ( sum( sum( scene( planet_pos_1 - l2 : planet_pos_1 + l2, planet_pos_2 - l2 : planet_pos_2 + l2 ) ) ) - scene( planet_pos_1, planet_pos_2 ) ) / ( ( 2 * l2 + 1 ) * ( 2 * l2 + 1 ) - 1 ) ;
    % Planet's flux
    plnt_flx( i_plnt ) = scene( planet_pos_1, planet_pos_2 ) - bckgrnd_tmp ;
    scene( planet_pos_1, planet_pos_2 ) = bckgrnd_tmp ;
    end % i_plnt
  % Keeping track of the coordinates of the planets removed
  opt.plnt.rmv.arc_dc_px = a - ( size( scene, 1 ) + 1 ) / 2 ;
  opt.plnt.rmv.arc_ra_px = b - ( size( scene, 1 ) + 1 ) / 2 ;
  opt.plnt.add.flx_rt = plnt_flx / opt.str.flx( i_lmbd ) ;
  disp( sprintf( '%i planets removed', n_plnt ) ) 
  return
  end

% Add a return to each case. If it arrives here is due to an error in the settings. For instance, remove_planets is et, but there are in fact no planets to remove.
disp( 'The action remove planets was set, but there is no identified method for this type of scene. Stopped.' )
disp( ' ' )
disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop

function scn_dt = add_extragalactic_background( scn_dt, opt, i_lmbd )
% Units are Jy
% Pixel pitch in Haystacks extragalactic background fields: 10 mas (FITS Header, or PASP paper. As of 12/12/18, the online notebook had 9 mas).
px_hystcks_bckgrnd_mas = 10 ;

% Particular feature from the Haystacks project (Background data have a different array size than the planetary files, and its number of elements is even, 3600x3600, instead of odd as in its scenes)
  if strcmp( opt.bckgrnd.img, 'GALAXIES_10lat' )
  % Making sure there's enough background data for the FOV
  scn_dt_sz_1 = size( scn_dt, 1 ) ;
  scn_dt_sz_2 = size( scn_dt, 2 ) ;
    if ( scn_dt_sz_1 * opt.px_scn_mas > 3600 * px_hystcks_bckgrnd_mas ) || ( scn_dt_sz_2 * opt.px_scn_mas > 3600 * px_hystcks_bckgrnd_mas )
    disp( 'There is not enough data in the extragalactic background field to cover the FOV of the scene. No extragalactic background data added.' )
    return
    end
  % Following the same steps as in get_scene_data
  % Getting the list of wavelengths in the Haystacks project (it brings lambda_array_scene)
  load( [ opt.scn_dr 'modern_cube_zodi1inc60dist10_misc.mat' ] ) ; % All have the same wavelength range
  % Each Haystacks file has 64 scenes
  n_scn_fl = 64 ;
  % Label of each FITS file
  lbl_scn = { '0.30-0.37um.fits', '0.37-0.46um.fits', '0.46-0.57um.fits', '0.57-0.70um.fits', '0.70-0.87um.fits', '0.87-1.07um.fits', '1.07-1.32um.fits', '1.32-1.63um.fits', '1.64-2.02um.fits', '2.02-2.50um.fits' } ;
  n_scn_lbl = numel( lbl_scn ) ;
  % Finding which files to open
  d_l = lambda_array_scene * 1000 - opt.lmbd_tmp_nm ;
  q_l = find( abs( d_l ) == min( abs( d_l ) ) ) ;
  d_l_q = d_l( q_l ) ;
    % If it's exactly the same as one of the wavelength in the scene data, it can be read off directly
    if ~d_l_q
    idx_fl = 1 + floor( ( q_l - 1 ) / n_scn_fl ) ; % Matlab indices are 1 ... N, for a N-element array
    idx_scn = mod( q_l, n_scn_fl ) ;
    % Reading off the data
    scn_dt_tmp = fitsread( [ opt.scn_dr 'GALAXIES_10lat' lbl_scn{ idx_fl } ], 'Image', idx_scn ) ; 
    else % Interpolation
      if d_l_q > 0
      q_l_1 = q_l - 1 ;
        % Extremely rare case where the imaging wavelength is less than the minimum scene wavelength
        if ~q_l_1
        q_l_1 = q_l ;
        end
      q_l_2 = q_l ;
      else
      q_l_1 = q_l ;
      q_l_2 = q_l + 1 ;
        % Extremely rare case where the imaging wavelength is greater than the maximum scene wavelength
        if ( q_l_2 > ( n_scn_fl * n_scn_lbl ) )
        q_l_2 = q_l ;
        end
      end
    idx_fl_1 = 1 + floor( ( q_l_1 - 1 ) / n_scn_fl ) ; % Matlab indices are 1 ... N for an N-element array
    idx_scn_1 = 1 + mod( q_l_1 - 1, n_scn_fl ) ; % Matlab indices are 1 ... N
    scn_dt_tmp = fitsread( [ opt.scn_dr 'GALAXIES_10lat' lbl_scn{ idx_fl_1 } ], 'Image', idx_scn_1 ) ;
    % Most common case. If q_l_1 == q_l_2, then no extrapolation, just identify the scene with the imaging wavelength. This would fail if someone images a scene at a wavelength well beyond the available simulated wavelength in the scenes, which in itself, is a total failure and should not happen.
      if ( q_l_1 ~= q_l_2 )
      idx_fl_2 = 1 + floor( ( q_l_2 - 1 ) / n_scn_fl ) ;
      idx_scn_2 = mod( q_l_2, n_scn_fl ) ;
      scn_dt_tmp_2 = fitsread( [ opt.scn_dr 'GALAXIES_10lat' lbl_scn{ idx_fl_2 } ], 'Image', idx_scn_2 ) ;
      % Linear interpolation
      scn_dt_tmp = ( ( lambda_array_scene( q_l_2 ) * 1000 - opt.lmbd_tmp_nm ) * scn_dt_tmp + ( opt.lmbd_tmp_nm - lambda_array_scene( q_l_1 ) * 1000 ) * scn_dt_tmp_2 ) / 1000 / ( lambda_array_scene( q_l_2 ) - lambda_array_scene( q_l_1 ) ) ;
      end
    end    

  % Recentering the extragalactic background field at the new location
  scn_dt_tmp = fraccircshift( scn_dt_tmp, [ opt.bckgrnd.arc_dc_mas / px_hystcks_bckgrnd_mas, opt.bckgrnd.arc_ra_mas / px_hystcks_bckgrnd_mas ] ) ;

   % Changing the units if necessary (the extragalactic background provided by the Haystacks project is in Jy)
    if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
    scn_dt_tmp = scn_dt_tmp * opt.conv_Jy_to_W_m2_um( i_lmbd ) ;
    end

  % According to the Haystacks project documentation the extragalactic background images have a pixel size of 10 mas, so we need to interpolate to get the scene resolution size (https://github.com/mjrfringes/haystacks/blob/master/HaystacksProcessing.ipynb. Notice that the example in the notebook implies 9 mas for the pixel pitch, but Maxime Rizzo confirmed it should be 10 (email: 12/12/18, Pixel scale on extragalactic background field from Haystacks))
  cntr_bck = size( scn_dt_tmp, 1 ) / 2 ; % recall the extragalactic background scenes are 3600x3600 pixels
  hlf_bck_1 = floor( scn_dt_sz_1 * opt.px_scn_mas / px_hystcks_bckgrnd_mas / 2 ) ;
  hlf_bck_2 = floor( scn_dt_sz_2 * opt.px_scn_mas / px_hystcks_bckgrnd_mas / 2 ) ;
  scn_dt_tmp = scn_dt_tmp( cntr_bck - hlf_bck_1 : cntr_bck + hlf_bck_1, cntr_bck - hlf_bck_2 : cntr_bck + hlf_bck_2 ) ;
  % Factor to use in imresize (recall imresize needs to be compensated in intensity)
  scn_dt_tmp = imresize( scn_dt_tmp, [ scn_dt_sz_1, scn_dt_sz_1 ], 'bilinear' ) * ( ( 2 * hlf_bck_1 + 1 ) / scn_dt_sz_1 )^2 ;
  end % Haystacks extragalactic background: GALAXIES_10lat...

% Adding the scene and the extragalactic background
scn_dt = scn_dt + scn_dt_tmp ;

% Building the effective, normalized PSF for the imaging wavelength and for the spinning Starshade case
function [ psf_array opt_img ] = get_normalized_spinning_psf_array( i_bnd, opt_img )

% In all the applications considered so far, the PSF files correspond to the unperturbed occulter case, except for the star (on-axis). The case of a perturbed starshade is dealt with separately later on with run_new_locus, so all the PSF files here can be assumed to be the unperturbed starshade, until something changes in the future (perturbed starshade+off axis star due to misalignment). 

% Array of available wavelength values in the PSF basis (this is very fast, no issues if it's done for each wavelength. It's clearer to keep it here)
opt_img.lmbd_arry_psf_nm = opt_img.lmbd_psf_1_nm( i_bnd ) : opt_img.dlt_lmbd_psf_nm : opt_img.lmbd_psf_2_nm( i_bnd ) ;
  if opt_img.lmbd_arry_psf_nm( end ) ~= opt_img.lmbd_psf_2_nm( i_bnd )
  opt_img.lmbd_arry_psf_nm( end + 1 ) = opt_img.lmbd_arry_psf_nm( end ) + opt_img.dlt_lmbd_psf_nm ;
  end
% Subdirectory where the imaging basis is found
  if ~isfield( opt_img, 'sbdr_psf' )
  opt_img.sbdr_psf = get_subdir_psf( opt_img ) ;
  end

% Opening as many files as necessary for the interpolation
d_lmbd = opt_img.lmbd_arry_psf_nm - opt_img.lmbd_tmp_nm ;
q_fnd_0 = find( d_lmbd == 0 ) ;
% Basic geometric parameters
[ dummy_dst_m geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt_img ) ;
  if ~isfield( opt_img, 'geo_iwa_mas' )
  opt_img.geo_iwa_mas = geo_iwa_mas ;
  end

  if ~numel( q_fnd_0 )
  % Linear interpolation
  q_fnd = find( abs( d_lmbd ) == min( abs( d_lmbd ) ) ) ;
  q_fnd = q_fnd( 1 ) ;
    if d_lmbd( q_fnd ) < 0
    idx_1 = q_fnd( 1 ) ;
    else
    idx_1 = q_fnd( 1 ) - 1 ;
    end
  idx_2 = idx_1 + 1 ;
  % Loading the two corresponding files (sister_basis.m)
  % Name of the first filename
  fl_nm_1_0 = sprintf( '%s%sstarshade_spinning_psf_Nx%i', opt_img.psf_dr, opt_img.sbdr_psf, opt_img.nx_pupil_pix ) ;
    if isfield( opt_img, 'pupil_filename' )
      if strcmp( opt_img.pupil_filename, '0' ) == 1
      fl_nm_1_0 = [ fl_nm_1_0 '_ideal' ] ;
      end
    end
  fl_nm_1 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%04d_nm.mat', fl_nm_1_0, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( idx_1 ) ) ;
  fl_nm_1_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%04d_nm.mat', fl_nm_1_0, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( idx_1 ) ) ;
    % For delta_lambda_nm that is not an integer
    if opt_img.lmbd_arry_psf_nm( idx_1 ) ~= round( opt_img.lmbd_arry_psf_nm( idx_1 ) )
    fl_nm_1 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%3.1f_nm.mat', fl_nm_1_0, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( idx_1 ) ) ;
    fl_nm_1_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%3.1f_nm.mat', fl_nm_1_0, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( idx_1 ) ) ;
    end

    if ~strcmp( opt_img.starshade.nominal_filename, 'NI2' )
    fl_nm_1 = [ fl_nm_1( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    fl_nm_1_1mas = [ fl_nm_1_1mas( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    end
  
  % Name of the second filename
  fl_nm_2_0 = sprintf( '%s%sstarshade_spinning_psf_Nx%i', opt_img.psf_dr, opt_img.sbdr_psf, opt_img.nx_pupil_pix ) ;
    if isfield( opt_img, 'pupil_filename' )
      if strcmp( opt_img.pupil_filename, '0' ) == 1
      fl_nm_2_0 = [ fl_nm_2_0 '_ideal' ] ;
      end
    end
  fl_nm_2 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%04d_nm.mat', fl_nm_2_0, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( idx_2 ) ) ;
  fl_nm_2_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%04d_nm.mat', fl_nm_2_0, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( idx_2 ) ) ;
    % For delta_lambda_nm that is not an integer
    if opt_img.lmbd_arry_psf_nm( idx_2 ) ~= round( opt_img.lmbd_arry_psf_nm( idx_2 ) )
    fl_nm_2 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%3.1f_nm.mat', fl_nm_2_0, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( idx_2 ) ) ;
    fl_nm_2_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%3.1f_nm.mat', fl_nm_2_0, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( idx_2 ) ) ;
    end
  
    if ~strcmp( opt_img.starshade.nominal_filename, 'NI2' )
    fl_nm_2 = [ fl_nm_2( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    fl_nm_2_1mas = [ fl_nm_2_1mas( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    end

   % Advise to create the imaging PSF basis if it does not exist and neither does the one with 1 mas pixel size
    if exist( fl_nm_1, 'file' ) ~= 2
    % Creating the PSF basis with this pixel size if the reference case of 1 mas already exists
    % Notice how I set 1 / opt_img.px_scn_mas, imposing the 1 mas as the reference pixel size of the PSF.
      if exist( fl_nm_1_1mas, 'file' ) == 2
      disp( sprintf( 'The PSF file %s does not exist. Creating it.', fl_nm_1 ) )
      load( fl_nm_1_1mas )
      t_strt = tic ;
        for i_psf = 1 : size( starshade_spinning_psf_array, 1 )
        starshade_spinning_psf_array_for_scene( i_psf, :, : ) = psf_resize_odd( squeeze( starshade_spinning_psf_array( i_psf, :, : ) ), 1 / opt_img.px_scn_mas ) ;
        end
      clear starshade_spinning_psf_array
      starshade_spinning_psf_array = starshade_spinning_psf_array_for_scene ;
      save( fl_nm_1, 'starshade_spinning_psf_array' ) ;
      disp( sprintf( 'PSF file %s stored (%3f s)', fl_nm_1, toc( t_strt ) ) ) ;
      else
      lbl_tmp = '' ;
        if ( opt_img.px_psf_mas ~= 1 )
        lbl_tmp = ' (neither does it with 1 mas pixel size)' ;
        end
      disp( sprintf( 'The spinning imaging basis %s does not exist%s. Please, create either of these first with sister_basis. Stopped (type dbquit all to leave the debugging mode).', fl_nm_1, lbl_tmp ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  load( fl_nm_1 )
  starshade_spinning_psf_array_1 = starshade_spinning_psf_array ;
  idx_psf_1 = ( size( starshade_spinning_psf_array_1, 2 ) - 1 ) / 2 ;

  % Advise to create the imaging PSF basis if it does not exist and neither does the one with 1 mas pixel size
    if exist( fl_nm_2, 'file' ) ~= 2
    % Creating the PSF basis with this pixel size if the reference case of 1 mas already exists
    % Notice how I set 1 / opt_img.px_scn_mas, imposing the 1 mas as the reference pixel size of the PSF.
      if exist( fl_nm_2_1mas, 'file' ) == 2
      disp( sprintf( 'The PSF file %s does not exist. Creating it.', fl_nm_2 ) )
      load( fl_nm_2_1mas )
      t_strt = tic ;
        for i_psf = 1 : size( starshade_spinning_psf_array, 1 )
        starshade_spinning_psf_array_for_scene( i_psf, :, : ) = psf_resize_odd( squeeze( starshade_spinning_psf_array( i_psf, :, : ) ), 1 / opt_img.px_scn_mas ) ;
        end
      clear starshade_spinning_psf_array
      starshade_spinning_psf_array = starshade_spinning_psf_array_for_scene ;
      save( fl_nm_2, 'starshade_spinning_psf_array' ) ;
      disp( sprintf( 'PSF file %s stored (%3f s)', fl_nm_2, toc( t_strt ) ) ) ;
      else
      disp( sprintf( 'The spinning imaging basis %s does not exist (neither does it with 1 mas pixel size). Please, create either of these first with sister_basis. Stopped (type dbquit all to leave the debugging mode).', fl_nm_2 ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  load( fl_nm_2 )
  cntr_psf_2 = ( size( starshade_spinning_psf_array, 2 ) - 1 ) / 2  + 1 ;
  dlt_lmbd_tmp = d_lmbd( idx_2 ) - d_lmbd( idx_1 ) ;
  % It should be positive
    if ( dlt_lmbd_tmp < 0 )
    disp( 'The amplitude of wavelength for the interpolation should be positive. Stopped' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  % Linear interpolation
  starshade_spinning_psf_array_1 = ( d_lmbd( idx_2 ) * starshade_spinning_psf_array_1 - d_lmbd( idx_1 ) * starshade_spinning_psf_array( :, cntr_psf_2 - idx_psf_1 : cntr_psf_2 + idx_psf_1, cntr_psf_2 - idx_psf_1 : cntr_psf_2 + idx_psf_1 ) ) / dlt_lmbd_tmp ;
  starshade_spinning_psf_array = starshade_spinning_psf_array_1 ;
  else
  % The wavelength is one of the PSF array. No interpolation is necessary
  fl_nm = sprintf( '%s%sstarshade_spinning_psf_Nx%i', opt_img.psf_dr, opt_img.sbdr_psf, opt_img.nx_pupil_pix ) ;
    if isfield( opt_img, 'pupil_filename' )
      if strcmp( opt_img.pupil_filename, '0' ) == 1
      fl_nm = [ fl_nm '_ideal' ] ;
      end
    end
  fl_nm_end = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%04d_nm.mat', fl_nm, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
  fl_nm_end_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%04d_nm.mat', fl_nm, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
    % For delta_lambda_nm that is not an integer
    if opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ~= round( opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) )    
    fl_nm_end = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%3.1f_nm.mat', fl_nm, geo_iwa_mas, opt_img.px_psf_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
    fl_nm_end_1mas = sprintf( '%s_geo_iwa_%03.2f_mas_pix_1_mas_%3.1f_nm.mat', fl_nm, geo_iwa_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
    end

    if ~strcmp( opt_img.starshade.nominal_filename, 'NI2' )
    fl_nm_end = [ fl_nm_end( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    fl_nm_end_1mas = [ fl_nm_end_1mas( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
    end

  % Advise to create the imaging basis if it does not exist (because these are many files)
    if exist( fl_nm_end, 'file' ) ~= 2
    % Creating the PSF basis with this pixel size if the reference case of 1 mas already exists
    % Notice how I set 1 / opt_img.px_scn_mas, imposing the 1 mas as the reference pixel size of the PSF.
      if exist( fl_nm_end_1mas, 'file' ) == 2
      disp( sprintf( 'The PSF file %s does not exist. Creating it.', fl_nm_end ) )
      load( fl_nm_end_1mas )
      t_strt = tic ;
        for i_psf = 1 : size( starshade_spinning_psf_array, 1 )
        starshade_spinning_psf_array_for_scene( i_psf, :, : ) = psf_resize_odd( squeeze( starshade_spinning_psf_array( i_psf, :, : ) ), 1 / opt_img.px_scn_mas ) ;
        end
      clear starshade_spinning_psf_array
      starshade_spinning_psf_array = starshade_spinning_psf_array_for_scene ;
      save( fl_nm_end, 'starshade_spinning_psf_array' ) ;
      disp( sprintf( 'PSF file %s stored (%3f s)', fl_nm_end, toc( t_strt ) ) ) ;
      else
      lbl_tmp = '' ;
        if ( opt_img.px_psf_mas ~= 1 )
        lbl_tmp = ' (neither does it with 1 mas pixel size)' ;
        end
      disp( sprintf( 'The spinning imaging basis %s does not exist%s. Please, create either of these first with sister_basis. Stopped (type dbquit all to leave the debugging mode).', fl_nm_end, lbl_tmp ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  load( fl_nm_end )
  end

  % If the pixel size of the PSF is different to that of the scene, the PSF has to be reduced to the resolution of the scene
  if ( opt_img.px_psf_mas ~= opt_img.px_scn_mas )    
  % Only checking if there's such basis in the case that the imaging is set for the same set of wavelengths as the basis (it means it was a basic configuration, instead of some exploration case)  
  do_rsz = 1 ;
    if numel( q_fnd_0 )
    fl_nm_2 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%04d_nm.mat', fl_nm, geo_iwa_mas, opt_img.px_scn_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
      % For delta_lambda_nm that is not an integer
      if opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ~= round( opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) )
      fl_nm_2 = sprintf( '%s_geo_iwa_%03.2f_mas_pix_%i_mas_%3.1f_nm.mat', fl_nm, geo_iwa_mas, opt_img.px_scn_mas, opt_img.lmbd_arry_psf_nm( q_fnd_0( 1 ) ) ) ;
      end

      if ~strcmp( opt_img.starshade.nominal_filename, 'NI2' )
      fl_nm_2 = [ fl_nm_2( 1 : end - 4 ) '_' opt_img.starshade.nominal_filename '.mat' ] ;
      end
      if exist( fl_nm_2, 'file' ) == 2
      load( fl_nm_2 ) 
      do_rsz = 0 ;
      end
    end
    if ( do_rsz )
    t_strt = tic ;
    disp( sprintf( 'The pixel size of the scene (%i mas) is different to that of the PSF basis (%i mas). Creating a version of the PSF basis with the same pixel size as the scene.', opt_img.px_scn_mas, opt_img.px_psf_mas ) ) 
      for i_psf = 1 : size( starshade_spinning_psf_array, 1 )
      starshade_spinning_psf_array_for_scene( i_psf, :, : ) = psf_resize_odd( squeeze( starshade_spinning_psf_array( i_psf, :, : ) ), opt_img.px_psf_mas / opt_img.px_scn_mas ) ;
      end
    clear starshade_spinning_psf_array
    starshade_spinning_psf_array = starshade_spinning_psf_array_for_scene ;
    % Storing the new basis if the set of wavelengths is the same as the PSF basis stored
      if numel( q_fnd_0 )
      save( fl_nm_2, 'starshade_spinning_psf_array' ) ;
      disp( sprintf( 'PSF file %s stored', fl_nm_2 ) ) ;
      toc( t_strt )
      end
    end
  end % opt_img.px_psf_mas ~= opt_img.px_scn_mas

% Normalizing the PSF
psf_st_nrm = sum( sum( squeeze( starshade_spinning_psf_array( end, :, : ) ) ) ) ;
psf_array = starshade_spinning_psf_array / psf_st_nrm ;
opt_img.psf_st_nrm = psf_st_nrm ; % used in run_new_locus when dealing with a perturbed locus
% Storing the transmission curve
i_lmbd = find( opt_img.lmbd_tmp_nm == opt_img.lmbd_arry_img_nm ) ;
opt_img.transmission_curve( i_lmbd, : ) = sum( sum( psf_array, 2 ), 3 ) ;
% Distance where it becomes 50%
r_0p1mas =  0 : 0.1 : opt_img.r_st_mas ;
% This is the spacing between 2 PSFs in mas
psf_stp_mas = opt_img.r_st_mas / ( numel( opt_img.transmission_curve( i_lmbd, : ) ) - 1 ) ;
trns_tmp = interp1( 0 : psf_stp_mas  : opt_img.r_st_mas, opt_img.transmission_curve( i_lmbd, : ), r_0p1mas ) ;
q_tmp = find( abs( trns_tmp - 0.5 ) == min( abs( trns_tmp - 0.5 ) ) ) ;
opt_img.half_transmission_mas( i_lmbd ) = r_0p1mas( q_tmp( 1 ) ) ;

% Building the effective, normalized PSF for the imaging wavelength and for the non-spinning Starshade case
function [ psf_nrm opt ] = get_normalized_psf( i_bnd, opt, x_mas, y_mas, store_pupil_psf_peak )

  if ~exist( 'store_pupil_psf_peak', 'var' )
  store_pupil_psf_peak = 0 ;
  end

% In all the applications so far, the PSF files correspond to the unperturbed occulter case. The case of peturbed Starshade has appreciable effects on the star (on-axis source) and that is dealt with run_new_locus later on.

% Array of wavelength available in the PSF (this is very fast, no issues if it's done for each wavelength. It's clearer to keep it here). PS: it has to be the same as in the files with the imaging basis information. PPS: it is by construction, when creating the basis and by the filename structure.
opt.lmbd_arry_psf_nm = opt.lmbd_psf_1_nm( i_bnd ) : opt.dlt_lmbd_psf_nm : opt.lmbd_psf_2_nm( i_bnd ) ;
  if opt.lmbd_arry_psf_nm( end ) ~= opt.lmbd_psf_2_nm( i_bnd )
  opt.lmbd_arry_psf_nm( end + 1 ) = opt.lmbd_arry_psf_nm( end ) + opt.dlt_lmbd_psf_nm ;
  end

% Considering as many cases as necessary for the interpolation
d_lmbd = opt.lmbd_arry_psf_nm - opt.lmbd_tmp_nm ;
q_fnd = find( d_lmbd == 0 ) ;
% Basic geometric parameters
[ dummy_dst_m geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt, 1 ) ;
  if ~isfield( opt, 'geo_iwa_mas' )
  opt.geo_iwa_mas = geo_iwa_mas ;
  end      
  if ~numel( q_fnd )
  % Linear interpolation
  q_fnd = find( abs( d_lmbd ) == min( abs( d_lmbd ) ) ) ;
  q_fnd = q_fnd( 1 ) ;
    if d_lmbd( q_fnd ) < 0
    idx_1 = q_fnd( 1 ) ;
    else
    idx_1 = q_fnd( 1 ) - 1 ;
    end
  idx_2 = idx_1 + 1 ;
  % Creating the two corresponding PSF
  % First PSF
  psf_1 = create_normalized_psf_from_E_field( i_bnd, opt, opt.lmbd_arry_psf_nm( idx_1 ), x_mas, y_mas, 0 ) ;
  idx_psf_1 = ( size( psf_1, 2 ) - 1 ) / 2 ;
  % Second PSF
  psf_2 = create_normalized_psf_from_E_field( i_bnd, opt, opt.lmbd_arry_psf_nm( idx_2 ), x_mas, y_mas, 0 ) ;
  sz_psf_2 = size( psf_2, 1 ) ;
  cntr_psf_2 = ( sz_psf_2 + 1 ) / 2 ;
  dlt_lmbd_tmp = d_lmbd( idx_2 ) - d_lmbd( idx_1 ) ;
  % It should be positive
    if ( dlt_lmbd_tmp < 0 )
    disp( 'The amplitude of wavelength for the interpolation should be positive. Stopped' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % Linear interpolation
   psf_nrm = ( d_lmbd( idx_2 ) * psf_1 - d_lmbd( idx_1 ) * psf_2( cntr_psf_2 - idx_psf_1 : cntr_psf_2 + idx_psf_1, cntr_psf_2 - idx_psf_1 : cntr_psf_2 + idx_psf_1 ) ) / dlt_lmbd_tmp ;
  % Only do this once setting store_pupil_psf_peak to 1 for the stationary PSF  
    if ( store_pupil_psf_peak )
    opt.psf_st_nrm = sum( sum( psf_nrm ) ) ;
    psf_nrm = psf_nrm / opt.psf_st_nrm ;
    end
  else
  % The wavelength is one of the PSF array. No interpolation is necessary
  [ psf_nrm opt ] = create_normalized_psf_from_E_field( i_bnd, opt, opt.lmbd_tmp_nm, x_mas, y_mas, store_pupil_psf_peak ) ;
  end
  
  % If the pixel size of the PSF is different to that of the scene, the PSF has to be reduced to the resolution of the scene
  if ( opt.px_psf_mas ~= opt.px_scn_mas )
    if ( opt.verbose )
    tic
    disp( sprintf( 'The pixel size of the scene (%i mas) is different to that of the PSF basis (%i mas). The process could be sped up if the same pixel size as the scene is used.' , opt.px_scn_mas, opt.px_psf_mas ) )
    end
  psf_nrm = psf_resize_odd( psf_nrm, opt.px_psf_mas / opt.px_scn_mas ) ;
    if ( opt.verbose )
    toc
    end
  end

% Function to build the PSF from the stored E fields before the pupil, and normalized to the peak of the stationary PSF
function [ psf_nrm opt ] = create_normalized_psf_from_E_field( i_bnd, opt, lmbd_nm, x_mas, y_mas, store_pupil_psf_peak )
% For the case of an ideal pupil (i.e., without struts) we can exploit the circular symmetry:
  if isfield( opt, 'pupil_filename' )
    if strcmp( opt.pupil_filename, '0' )
    % Reducing the position to that within the first petal (assumed to be aligned along x -matlab- axis, consistent with starshade_spinning_local in sister_basis when dealing with the radial PSF)
    r_tmp_mas = round( sqrt( ( x_mas^2 + y_mas^2 ) ) ) ;
    alph_tmp_rd = atan2( y_mas, x_mas ) ;
    % Angle corresponding to a petal
    alph_ptl_rd = 2 * pi  / opt.n_ptl ;
    % Corresponding angle within one petal
    alph_tmp_2_rd = alph_tmp_rd - floor( ( alph_tmp_rd  + alph_ptl_rd  / 2 ) / alph_ptl_rd ) * alph_ptl_rd ;
    % Corresponding x, y positions
    x_mas = r_tmp_mas * cos( alph_tmp_2_rd ) ;
    y_mas = r_tmp_mas * sin( alph_tmp_2_rd ) ;
    % The PSF is rotated back after being created
    end
  end
% Getting the filename (Quadrant symmetry for the E fields).
opt_2 = opt ;
opt_2.x_source_mas = abs( opt.px_psf_mas * round( x_mas / opt.px_psf_mas ) ) ; % Notice the reduction to the pixel size of the PSF
opt_2.y_source_mas = abs( opt.px_psf_mas * round( y_mas / opt.px_psf_mas ) ) ;
opt_2.lambda_1_nm = opt.lmbd_psf_1_nm( i_bnd ) ;
opt_2.lambda_2_nm = opt.lmbd_psf_2_nm( i_bnd ) ;
opt_2.path_occulter = [ opt.scn_dr 'locus/in/' ] ;
  if ~isfield( opt, 'sbdr_psf' )
  opt.sbdr_psf = get_subdir_psf( opt ) ;
  end
opt_2.save_path = [ opt.psf_dr opt.sbdr_psf ] ;
opt_2 = get_default_options( opt_2 ) ;
fl_psf = [ opt_2.save_path opt_2.save_filename '.mat' ] ;
  if exist( fl_psf, 'file' ) ~= 2
  disp( sprintf( 'The file %s does not exist. Creating it. * Notice that sister_basis is designed to create a full basis * This is meant for one or few PSF.', fl_psf ) )
  opt_2.save = 1 ;
  opt_2.make_occulter_name = opt.starshade.nominal_filename ;
  makeStarshadeImage( opt_2 ) ;
  end
% It brings UTotL, pupil and opt_make, if necessary
load( [ opt.psf_dr opt.sbdr_psf opt_2.save_filename '.mat' ] ) 
% Select the single wavelength (float precision rounding errors require some rounding. In this case we set coincidence within 1 nm)
q_wl = find( round( lambdaIn * 1e9 ) == round( lmbd_nm ) ) ;
  if numel( q_wl ) ~= 1
  disp( sprintf( 'Error: the wavelength of the PSF %04.4f nm has no correspondence in the imaging basis, and it must. Stopped', lmbd_nm ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
UTot = squeeze( UTotL( :, :, q_wl ) ) ;
clear UTotL
lambdaIn = lmbd_nm / 1e9 ;
% Erasing the hot pixels (if any). PS: set to 1 iteration. It should be enough. If there's any non-sense PSF then the number of iterations needs be increased. The appearance of hot pixels is something that should be removed when the core script makeStarshadeImage.m derives the electric fields at the pupil plane (E. Cady will spend some time if it turns out to be an issue).
UTot = erase_hot_pixels( UTot, 1, 0 ) ;
% Re-orienting the E fields (consistent with sister_basis)
  if ( x_mas < 0 )
  UTot = flipud( UTot ) ;
  end
  if ( y_mas < 0 )
  UTot = fliplr( UTot ) ;
  end
% Here, I extract what is necessary from E. Cady's code:
opt.lD2mas = lmbd_nm / 1e9 / opt.dmtr_tlscp_m * 180 / pi * 3600 * 1000 ;
% Number of pixels around the PSF center to be used (from analysis on simulations with NI2, 6.5 l/D keeps relative errors below 1e-3 in accordance with the choice of opt.n_lambda_over_d=7 for the NI2 PSF basis. Every PSF basis has opt.n_lambda_over_d chosen to provide some target precision).
opt.n_px_psf = ceil( opt.n_lambda_over_d * opt.lD2mas / opt.px_psf_mas ) ;
Nx_img = 2 * opt.n_px_psf + 1 ;
imagePlaneDiameterInLambdaOverD = Nx_img * opt.px_psf_mas / opt.lD2mas ;
imagePlane = mft_shift( imagePlaneDiameterInLambdaOverD, Nx_img, UTot .* pupil, x_mas, y_mas, opt.px_psf_mas ) ;
int_src = abs( imagePlane ).^2 ; 

if ( 0 )
% Proper matlab ordering
% Computing the centroid 
[ y_tmp x_tmp ] = meshgrid( 1 : Nx_img, 1 : Nx_img ) ;
x_cntrd = sum( sum( x_tmp .* int_src ) ) / sum( sum( int_src ) ) ; 
y_cntrd = sum( sum( y_tmp .* int_src ) ) / sum( sum( int_src ) ) ;
% Computing the pixel where the maximum intensity happens
q_mx = find( int_src == max( int_src( : ) ) ) ;
y_mx = ceil( q_mx / size( int_src, 1 ) ) ;
x_mx = mod( q_mx, size( int_src, 1 ) ) ;
% Finding a new centroid in the vicinity of the PSF's peak
l = 20 ;
psf_tmp = int_src( x_mx - l : x_mx + l, y_mx - l : y_mx + l ) ;
[ y_tmp x_tmp ] = meshgrid( - l : l, -l : l ) ;
y_mx_cntrd = y_mx + sum( sum( y_tmp .* psf_tmp ) ) / sum( sum( psf_tmp ) ) ;
x_mx_cntrd = x_mx + sum( sum( x_tmp .* psf_tmp ) ) / sum( sum( psf_tmp ) ) ;
disp( 'Centroid' )
disp( sprintf( '%3.3f, %3.3f', x_cntrd - ( Nx_img + 1 ) / 2, y_cntrd  - ( Nx_img + 1 ) / 2 ) ) ;
disp( 'Max' )
disp( sprintf( '%3.3f, %3.3f', x_mx - ( Nx_img + 1 ) / 2, y_mx - ( Nx_img + 1 ) / 2 ) ) ;
disp( 'Max, centroid' )
disp( sprintf( '%3.3f, %3.3f', x_mx_cntrd - ( Nx_img + 1 ) / 2, y_mx_cntrd - ( Nx_img + 1 ) / 2 ) ) ;
disp( 'Initial-Max. Centroid' )
disp( sprintf( '%3.3f, %3.3f', x_mx_cntrd - round( x_mas / opt.px_psf_mas ) - ( Nx_img + 1 ) / 2, y_mx_cntrd - round( y_mas / opt.px_psf_mas ) - ( Nx_img + 1 ) / 2 ) ) ;
disp( '' )
%disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
end

% Rotating the PSF back to its position in the case of an ideal pupil, since we have used the reduction to within a petal above
  if isfield( opt, 'pupil_filename' )
    if strcmp( opt.pupil_filename, '0' )
    % Checked that rotating by a + angle in imrotate it does rotate the PSF to the right place
    int_src = imrotate( int_src, 180 / pi * floor( ( alph_tmp_rd  + alph_ptl_rd  / 2 ) / alph_ptl_rd ) * alph_ptl_rd, 'crop', 'bilinear' ) ;
    end
  end

% Only do this once setting store_pupil_psf_peak to 1 for the stationary PSF
  if ( store_pupil_psf_peak )
  opt.psf_st_nrm = sum( sum( int_src ) ) ;
  end
% Normalizing the PSF to the peak of the stationary PSF (that's why the previous psf_peak is computed only once, when the switch is set 1 for the case of the stationary PSF)
psf_nrm = int_src / opt.psf_st_nrm ;

% Reading the cubes from the Haystacks project: https://asd.gsfc.nasa.gov/projects/haystacks/haystacks.html and https://arxiv.org/abs/1710.06328
function scn_dt = read_haystacks_project( opt )
% Getting the list of wavelengths in the Haystacks project (it brings lambda_array_scene)
load( [ opt.scn_dr 'modern_cube_zodi1inc60dist10_misc.mat' ] ) ; % All Haystacks scenes have the same wavelength range
% Each Haystacks file has 64 scenes
n_scn_fl = 64 ;
% Label of each FITS file
lbl_scn = { '_0.30-0.37um.fits', '_0.37-0.46um.fits', '_0.46-0.57um.fits', '_0.57-0.70um.fits', '_0.70-0.87um.fits', '_0.87-1.07um.fits', '_1.07-1.32um.fits', '_1.32-1.63um.fits', '_1.64-2.02um.fits', '_2.02-2.50um.fits' } ;
n_scn_lbl = numel( lbl_scn ) ;
% Finding which files to open
d_l = lambda_array_scene * 1000 - opt.lmbd_tmp_nm ;
q_l = find( abs( d_l ) == min( abs( d_l ) ) ) ;
d_l_q = d_l( q_l ) ;
  % If it's exactly the same as one of the wavelength in the scene data, it can be read off directly
  if ~d_l_q
  idx_fl = 1 + floor( ( q_l - 1 ) / n_scn_fl ) ; % Matlab indices are 1 ... N, for a N-element array
  idx_scn = mod( q_l, n_scn_fl ) ;
  % Reading off the data
   tic
scn_dt = fitsread( [ opt.scn_dr opt.scn.nm lbl_scn{ idx_fl } ], 'Image', idx_scn ) ;
   disp( 'Time reading the scene data' ) ; toc % Sub-second
    if ( 1 ) % It works
    % Check of consistency
    fts_info = fitsinfo( [ opt.scn_dr opt.scn.nm lbl_scn{ idx_fl } ] ) ;
    lmbd_info = fts_info.Image( idx_scn ).Keywords( 9, 2 ) ;
    % disp( 'Time reading the FITS info' ) ; toc % Sub-second
      if ( abs( lmbd_info{ 1 } * 1000 - opt.lmbd_tmp_nm ) > 0.5 )
      disp( 'The wavelength from the scene to be opened and the one read seem to be different. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end % if (0/1)
  else % Case of a wavelength in the image list that is not exactly a wavelength of the scene data
    if d_l_q > 0
    q_l_1 = q_l - 1 ;
    % Extremely rare case where the imaging wavelength is less than the minimum scene wavelength
      if ~q_l_1
      q_l_1 = q_l ;
      end
    q_l_2 = q_l ;
    else
    q_l_1 = q_l ;
    q_l_2 = q_l + 1 ;
      % Extrenely rare case where the imaging wavelength is greater than the maximum scene wavelength
      if ( q_l_2 > ( n_scn_fl * n_scn_lbl ) )
      q_l_2 = q_l ;
      end
    end
  idx_fl_1 = 1 + floor( ( q_l_1 - 1 ) / n_scn_fl ) ; % Matlab indices are 1 ... N for an N-element array
  idx_scn_1 = 1 + mod( q_l_1 - 1, n_scn_fl ) ; % Matlab indices are 1 ... N
  scn_dt = fitsread( [ opt.scn_dr opt.scn.nm lbl_scn{ idx_fl_1 } ], 'Image', idx_scn_1 ) ;
  % Most common case. If q_l_1 == q_l_2, then no extrapolation, just identify the scene with the imaging wavelength. This would fail if someone images a scene at a wavelength well beyond the available simulated wavelength in the scenes, which in itself, is a total failure and should not happen.
    if ( q_l_1 ~= q_l_2 )
    idx_fl_2 = 1 + floor( ( q_l_2 - 1 ) / n_scn_fl ) ;
    idx_scn_2 = mod( q_l_2, n_scn_fl ) ;
    scn_dt_2 = fitsread( [ opt.scn_dr opt.scn.nm lbl_scn{ idx_fl_2 } ], 'Image', idx_scn_2 ) ;
    % Linear interpolation
    scn_dt = ( ( lambda_array_scene( q_l_2 ) * 1000 - opt.lmbd_tmp_nm ) * scn_dt + ( opt.lmbd_tmp_nm - lambda_array_scene( q_l_1 ) * 1000 ) * scn_dt_2 ) / 1000 / ( lambda_array_scene( q_l_2 ) - lambda_array_scene( q_l_1 ) ) ;
      % Check of consistency
      if ( lambda_array_scene( q_l_2 ) - lambda_array_scene( q_l_1 ) ) <= 0
      disp( 'The different of wavelength is expected to be positive and it is not. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    else
    disp( sprintf( '*** Odd case: the wavelength required to be imaged is outside the range of the simulated scenes. Imaging wavelength is %4.2d nm', opt.lmbd_tmp_nm ) )
    end
  end % Imaging wavelength different to scene wavelength

% Adding the planets
% 1) The case of a spinning Starshade
function [ scn_dt_cnv opt ] = add_planets_spinning( scn_dt_cnv, psf_array, opt, i_lmbd )

% Special case when planets are removed and added at the same positions to simulate a background estimation
  if ( opt.plnt.rmv.do == -1 )
  opt.plnt.add.arc_dc_px = opt.plnt.rmv.arc_dc_px ;
  opt.plnt.add.arc_ra_px = opt.plnt.rmv.arc_ra_px ;
  opt.plnt.add.arc_dc_mas = opt.plnt.rmv.arc_dc_px * opt.px_scn_mas ;
  opt.plnt.add.arc_ra_mas = opt.plnt.add.arc_ra_px * opt.px_scn_mas ;
  end

% Getting the number of planets
n_pl = numel( opt.plnt.add.arc_dc_px ) ;
disp( sprintf( 'Adding %i planets', n_pl ) )

  for i_pl = 1 : n_pl
  % NB: 2e-3 seconds spent to do this (fast)
  [ plnt_psf opt ] = get_planet_spinning_psf( opt, psf_array, scn_dt_cnv, i_pl, i_lmbd ) ;
  % In case the planet does not fall within the FOV (see get_planet_non_spinning_psf.m and get_planet_spinning_psf.m below)
    if isnan( sum( plnt_psf( : ) ) )
    continue
    end
  scn_dt_cnv( opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ), opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) ) = ...
  scn_dt_cnv( opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ), opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) ) + ...
  opt.str.flx( i_lmbd ) * opt.plnt.add.flx_rt( i_pl, i_lmbd ) * plnt_psf ;


  % Storing the indices for later use (e.g., for the data cubes)
    if ( opt.cb.plnt.do )
    opt.cb.plnt.n_plnt_psf_1( i_pl, i_lmbd ) = size( plnt_psf, 1 ) ;
    opt.cb.plnt.n_plnt_psf_2( i_pl, i_lmbd ) = size( plnt_psf, 2 ) ;
    opt.cb.plnt.plnt_psf( i_pl, i_lmbd, 1 : opt.cb.plnt.n_plnt_psf_1( i_pl, i_lmbd ) , 1 : opt.cb.plnt.n_plnt_psf_2( i_pl, i_lmbd ) ) = plnt_psf ;
    end
  end

% Case of a non-spinning Starshade
function scn_dt_cnv = add_planets_non_spinning( scn_dt_cnv, psf_st, opt, i_lmbd, i_bnd )

% Special case when planets are removed and added at the same positions to simulate a background estimation
  if ( opt.plnt.rmv.do == -1 )
  opt.plnt.add.arc_dc_px = opt.plnt.rmv.arc_dc_px ;
  opt.plnt.add.arc_ra_px = opt.plnt.rmv.arc_ra_px ;
  opt.plnt.add.arc_dc_mas = opt.plnt.rmv.arc_dc_px * opt.px_scn_mas ;
  opt.plnt.add.arc_ra_mas = opt.plnt.add.arc_ra_px * opt.px_scn_mas ;
  end

% Getting the number of planets
n_pl = numel( opt.plnt.add.arc_dc_px ) ;
disp( sprintf( 'Adding %i planets', n_pl ) )

  for i_pl = 1 : n_pl
  [ plnt_psf opt ] = get_planet_non_spinning_psf( opt, psf_st, scn_dt_cnv, i_pl, i_lmbd, i_bnd ) ;
  scn_dt_cnv( opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ), opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) ) = ...
  scn_dt_cnv( opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ), opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) ) + ...
  opt.str.flx( i_lmbd ) * opt.plnt.add.flx_rt( i_pl, i_lmbd ) * plnt_psf ;

   % Storing the indices for later use (e.g., for the data cubes)
    if ( opt.cb.plnt.do )
    opt.cb.plnt.n_plnt_psf_1( i_pl, i_lmbd ) = size( plnt_psf, 1 ) ;
    opt.cb.plnt.n_plnt_psf_2( i_pl, i_lmbd ) = size( plnt_psf, 2 ) ;
    opt.cb.plnt.plnt_psf( i_pl, i_lmbd, 1 : opt.cb.plnt.n_plnt_psf_1( i_pl, i_lmbd ) , 1 : opt.cb.plnt.n_plnt_psf_2( i_pl, i_lmbd ) ) = plnt_psf ;
    end
  end

% Add a new locus coming from some perturbation
function scn_dt_cnv = add_new_locus( scn_dt_cnv, psf_0, opt, i_lmbd, psf_st )
opt.lcs.sz_psf_rf = size( psf_0, 1 ) ; % Assuming a square
opt.lcs.px_scn_mas = opt.px_scn_mas ; % This is the pixel scale of the nominal PSF
opt.lcs.starshade = opt.starshade ;
opt.lcs.px_psf_mas = opt.px_psf_mas ;
opt.lcs.geo_iwa_mas = opt.geo_iwa_mas ; 
opt.lcs.diam_img_mas = size( scn_dt_cnv, 1 ) * opt.px_scn_mas ;
opt.lcs.px_res_mas = min( opt.px_psf_mas, 1 ) ; % At most 1 mas when dealing with non-ideal starshades
opt.lcs.pupil_filename = opt.pupil_filename ;
opt.lcs.nx_pupil_pix = opt.nx_pupil_pix ;
opt.lcs.secondary_size = opt.secondary_size ;
opt.lcs.lmbd_tmp_nm = opt.lmbd_tmp_nm ;
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  opt.lcs.transmission_curve_tmp = opt.transmission_curve( find( opt.lmbd_tmp_nm == opt.lmbd_arry_img_nm ), : ) ;
  opt.lcs.half_transmission_mas = opt.half_transmission_mas( end ) ; % Only parse the current value
  end
opt.lcs.r_st_mas = opt.r_st_mas ;
opt.lcs.au2km = opt.au2km ;
opt.lcs.au2mas = opt.au2mas ;
% Using the xVals, yVals of the unperturbed locus. Respecting the filename convention
opt.lcs.unprtrbd = sprintf( '%s', opt.starshade.nominal_filename ) ;
% Dummy (necessary to pass get_default_options)
  if isfield( opt, 'diameter_telescope_m' )
  opt.lcs.diameter_telescope_m = opt.diameter_telescope_m ;
  else
  opt.lcs.diameter_telescope_m = opt.dmtr_tlscp_m ;
  end
  if isfield( opt, 'distance_starshade_telescope_m' )
  opt.lcs.distance_starshade_telescope_m = opt.distance_starshade_telescope_m ; 
  else
  opt.lcs.distance_starshade_telescope_m = opt.dst_strshd_tlscp_m ;
  end
disp( sprintf( 'Considering a new locus: %s', opt.lcs.prtrbd ) )
% Getting the proper normalization for the perturbed PSF
opt.lcs.mx_psf_nmnl = max( psf_0( : ) ) ;
% Generating the non-ideal PSF response. The setup for the perturbed locus is different than for the unperturbed one, in order to increase the spatial extension of the PSF since the perturbed PSF is brighter at larger distances than the unperturbed case. We calibrate the new perturbed on-axis response, keeping track of the stationary PSF generated in the same setup.
psf_cntr_prt = run_new_locus( opt.lcs ) ;
% After some testing, it's more accurate to use imresize to resize the image than the method from the mmf code in SITER. From perturbed pixel scale to scene pixel scale
sz_prt_scn = 2 * floor( size( psf_cntr_prt, 1 ) * opt.lcs.px_res_mas / opt.px_scn_mas / 2 ) + 1 ;
psf_cntr_prt = imresize( psf_cntr_prt, [ sz_prt_scn sz_prt_scn ] ) ;
% PS: run_new_locus already normalizes conveniently the PSF response on-axis. No need to re-scale imresize now.
% Subtracting the non-ideal PSF and adding the perturbed one
cntr_scn_dt_1 = ( size( scn_dt_cnv, 1 ) + 1 ) / 2 ;
cntr_scn_dt_2 = ( size( scn_dt_cnv, 2 ) + 1 ) / 2 ;

% Subtracting the reference PSF
hlf_psf_0_1 = ( size( psf_0, 1 ) - 1 ) / 2 ;
hlf_psf_0_2 = ( size( psf_0, 2 ) - 1 ) / 2 ;
idx_1_1 = cntr_scn_dt_1 - hlf_psf_0_1 ;
% Reference PSF within FOV?
mssg_tmp = 0 ;
  if ( idx_1_1 < 1 ), idx_1_1 = 1 ; mssg_tmp = 1 ; end
idx_1_2 = cntr_scn_dt_1 + hlf_psf_0_1 ;
  if ( idx_1_2 > size( scn_dt_cnv, 1 ) ), idx_1_2 = size( scn_dt_cnv, 1 ) ; mssg_tmp = 1 ; end
idx_2_1 = cntr_scn_dt_2 - hlf_psf_0_2 ;
  if ( idx_2_1 < 1 ), idx_2_1 = 1 ; mssg_tmp = 1 ; end
idx_2_2 = cntr_scn_dt_2 + hlf_psf_0_2 ;
   if ( idx_2_2 > size( scn_dt_cnv, 2 ) ), idx_2_2 = size( scn_dt_cnv, 2 ) ; mssg_tmp = 1 ; end

% If the PSF fits inside the FOV
  if ~( mssg_tmp )
  scn_dt_cnv( idx_1_1 : idx_1_2, idx_2_1 : idx_2_2 ) = scn_dt_cnv( idx_1_1 : idx_1_2, idx_2_1 : idx_2_2 ) - opt.str.flx( i_lmbd ) * psf_0( 1 : idx_1_2 - idx_1_1 + 1, 1 : idx_2_2 - idx_2_1 + 1 ) ;
  else
  % Warning message if the FOV is less than the extension of the PSF (lambda_over_d)
  disp( sprintf( 'WARNING: The Field of View is %04.1f mas. It is smaller than the nominal PSF domain %04.1f mas, which is set by opt.n_lambda_over_d(=%02.1f). Using the portion of the PSF within the FOV. To cover the full extension of the PSF, set opt.scene.fov_diam_mas bigger than or equal to the PSF domain.', opt.scene.fov_diam_mas, size( psf_0, 1 ) * opt.px_scn_mas, opt.n_lambda_over_d ) )
  cntr_psf_0_1 = ( size( psf_0, 1 ) + 1 ) / 2 ;
  cntr_psf_0_2 = ( size( psf_0, 2 ) + 1 ) / 2 ;
  hlf_scn_1 = ( size( scn_dt_cnv, 1 ) - 1 ) / 2 ;
  hlf_scn_2 = ( size( scn_dt_cnv, 2 ) - 1 ) / 2 ;
  scn_dt_cnv = scn_dt_cnv - opt.str.flx( i_lmbd ) * psf_0( cntr_psf_0_1 - hlf_scn_1 : cntr_psf_0_1 + hlf_scn_1, cntr_psf_0_2 - hlf_scn_2 : cntr_psf_0_2 + hlf_scn_2 ) ;
  end

% Adding the non-ideal PSF response
hlf_psf_prt_1 = ( size( psf_cntr_prt, 1 ) - 1 ) / 2 ;
hlf_psf_prt_2 = ( size( psf_cntr_prt, 2 ) - 1 ) / 2 ;
% Non-ideal PSF within FOV?
mssg_tmp = 0 ;
idx_1_1 = cntr_scn_dt_1 - hlf_psf_prt_1 ;
  if ( idx_1_1 < 1 ), idx_1_1 = 1 ; mssg_tmp = 1 ; end
idx_1_2 = cntr_scn_dt_1 + hlf_psf_prt_1 ;
  if ( idx_1_2 > size( scn_dt_cnv, 1 ) ), idx_1_2 = size( scn_dt_cnv, 1 ) ; mssg_tmp = 1 ; end
idx_2_1 = cntr_scn_dt_2 - hlf_psf_prt_2 ;
  if ( idx_2_1 < 1 ), idx_2_1 = 1 ; mssg_tmp = 1 ; end
idx_2_2 = cntr_scn_dt_2 + hlf_psf_prt_2 ;
   if ( idx_2_2 > size( scn_dt_cnv, 2 ) ), idx_2_2 = size( scn_dt_cnv, 2 ) ; mssg_tmp = 1 ; end

% If the PSF fits inside the FOV
  if ~( mssg_tmp )
  scn_dt_cnv( idx_1_1 : idx_1_2, idx_2_1 : idx_2_2 ) = scn_dt_cnv( idx_1_1 : idx_1_2, idx_2_1 : idx_2_2 ) + opt.str.flx( i_lmbd ) * psf_cntr_prt( 1 : idx_1_2 - idx_1_1 + 1, 1 : idx_2_2 - idx_2_1 + 1 ) ;
  else
  disp( sprintf( 'WARNING: The Field of View is %04.1f mas. It is smaller than the non-ideal PSF domain %04.1f mas, which is set by opt.n_lambda_over_d(=%02.1f). The non-ideal PSF extension has been increased a few times with respect the nominal PSF, because the starlight is brighter (see run_new_locus, opt.times_nominal). Using the portion of the PSF within the FOV. To cover the full extension of the PSF set opt.scene.fov_diam_mas bigger than or equal to the PSF domain.', opt.scene.fov_diam_mas, size( psf_cntr_prt, 1 ) * opt.px_scn_mas, opt.n_lambda_over_d ) )
  cntr_psf_prt_1 = ( size( psf_cntr_prt, 1 ) + 1 ) / 2 ;
  cntr_psf_prt_2 = ( size( psf_cntr_prt, 2 ) + 1 ) / 2 ;
  hlf_scn_1 = ( size( scn_dt_cnv, 1 ) - 1 ) / 2 ;
  hlf_scn_2 = ( size( scn_dt_cnv, 2 ) - 1 ) / 2 ;
  scn_dt_cnv = scn_dt_cnv + opt.str.flx( i_lmbd ) * psf_cntr_prt( cntr_psf_prt_1 - hlf_scn_1 : cntr_psf_prt_1 + hlf_scn_1, cntr_psf_prt_2 - hlf_scn_2 : cntr_psf_prt_2 + hlf_scn_2 ) ;
  end

% Adding some solar glint
function scn_dt_cnv = add_solar_glint( scn_dt_cnv, opt ) 

% Getting the apparent flux of the sun 
opt_tmp = opt ; opt_tmp.str.tp = 'sun' ; opt_tmp.str.mg = -26.76 ; opt_tmp = get_star_flux( scn_dt_cnv, opt_tmp, 1 ) ; % Absolute magnitude of the Sun from http://mips.as.arizona.edu/~cnaw/sun.html
% Only picking the first index because get_star_flux is run anew everytime for each new wavelength
sun_app_flx = opt_tmp.str.flx( 1 ) ; clear opt_tmp 

% Creating the array of relative fluxes (solar glint)
[ fluxgridst, fluxgrid ] = glint_distrib4imaging( opt, size( scn_dt_cnv ) ) ;
% Constructing the PSF from the Starshade
psf_ss = get_psf_from_starshade( opt ) ;
% Padding the PSF to the size of the convolved data
psf_ss_pd = 0 * scn_dt_cnv ;
cntr_cnv = ( size( scn_dt_cnv, 1 ) + 1 ) / 2 ;
hlf_psf_ss = ( size( psf_ss, 1 ) -1 ) / 2 ;
% It may happen that the PSF from the starshade is larger than the actual FOV. In that case, SISTER trims the PSF and displays a warning message.
  if ( cntr_cnv - hlf_psf_ss < 1 )
  nw_fov_mas = ceil( size( psf_ss, 1 ) * opt.px_scn_mas * max( opt.lmbd_arry_img_nm ) / opt.lmbd_tmp_nm ) ;
size( psf_ss, 1 ) * opt.px_scn_mas
max( opt.lmbd_arry_img_nm ) / opt.lmbd_tmp_nm
  psf_ss = psf_ss( 2 - cntr_cnv + hlf_psf_ss : end + cntr_cnv - hlf_psf_ss - 1, 2 - cntr_cnv + hlf_psf_ss : end + cntr_cnv - hlf_psf_ss - 1 ) ;
  hlf_psf_ss = ( size( psf_ss, 1 ) -1 ) / 2 ;
  disp( sprintf( 'WARNING: the PSF for the solar glint has been re-sized to the size of the FOV of the scene. If you want to cover the entire PSF from the solar glint, please, set opt.scene.fov_diam_mas=%4.0f mas or to a larger value.', nw_fov_mas ) )
  end
psf_ss_pd( cntr_cnv - hlf_psf_ss : cntr_cnv + hlf_psf_ss, cntr_cnv - hlf_psf_ss : cntr_cnv + hlf_psf_ss ) = psf_ss ;
% Convolving it with the SS PSF
ipsf_ss_pd = fft2( psf_ss_pd ) ;
  if ( opt.slr_glnt.stlth )
  ifluxgrid_tmp = fft2( fluxgridst ) ;
  else
  ifluxgrid_tmp = fft2( fluxgrid ) ;
  end
slr_glnt_fct = fftshift( ifft2( ifluxgrid_tmp .* ipsf_ss_pd ) ) ;
% Adding the absolute value of the solar glint
scn_dt_cnv = scn_dt_cnv + sun_app_flx * slr_glnt_fct ;

% Creating an exo-zodi for the case when the scene is built from some user specs
function scn_dt = add_exozodi( scn_dt, opt, i_lmbd )

% Temporary hack for the WFIRST Exoplanet Data Challenge
  if strcmp( lower( opt.run ), 'dc2' )
  % Realistic model created with Dustmap and Zodipic by Neil Zimmerman. Average over the 425-552 nm band at 448.5 nm (to be used with the 425-552 nm pass band). Intensity is 5 zodi and the pixel scale is 5.02 mas. Units are Jy. 
  zd_dc2 = fitsread( '/Users/srhildeb/Desktop/Starshade/WFIRST_DC2/DC2_101019/finalexozodi_SS_5zodi.fits' ) ; 
  % Choosing the corresponding epoch
  zd_dc2 = squeeze( zd_dc2( :, :, opt.dc2.i_rw ) ) ;
    % For the test image only
      if isfield( opt.dc2, 'test' )
      zd_dc2 = imrotate( mean( zd_dc2, 3 ), 34 ) ;
      % 06/24/19. 5 zodi, 5 mas pixel scale.
        if opt.dc2.test == 2
        zd_dc2 = fitsread( '/Users/srhildeb/Desktop/Starshade/WFIRST_DC2/DC2_062419/zodipic_ipachack_wav488nm_05zodi.fits' ) ;
        end
      end
  disp( '*** Exozodiacal model for the DC2 read.' )
  % Bringing it to the pixel scale of the simulated scene, with odd numbers. The total flux has to be the same, so imresize needs be scaled
  fct_rsz_zd = ( 2 * floor( 5 * size( zd_dc2, 1 ) / opt.px_scn_mas / 2 ) + 1 ) / size( zd_dc2, 1 ) ;
  zd_dc2 = imresize( zd_dc2, fct_rsz_zd, 'bilinear' ) / fct_rsz_zd^2 ;
  % Changing the units if necessary
    if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
    zd_dc2 = zd_dc2 * opt.conv_Jy_to_W_m2_um( i_lmbd ) ;
    end
  % Inserting the exozodiacal emission into the scene
  sz_scn_dt = size( scn_dt, 1 ) ;
  sz_dc2 = size( zd_dc2, 1 ) ;
    if ( sz_scn_dt < sz_dc2 )
    cntr_zd_dc2 = ( size( zd_dc2, 1 ) + 1 ) / 2 ;
    hlf_scn_dt = ( sz_scn_dt - 1 ) / 2 ;
    scn_dt = scn_dt + zd_dc2( cntr_zd_dc2 - hlf_scn_dt : cntr_zd_dc2 + hlf_scn_dt, cntr_zd_dc2 - hlf_scn_dt : cntr_zd_dc2 + hlf_scn_dt ) ;
    else
    cntr_scn_dt = ( size( scn_dt, 1 ) + 1 ) / 2 ;
    hlf_zd_dc2 = ( sz_dc2 - 1 ) / 2 ;
    scn_dt( cntr_scn_dt - hlf_zd_dc2 : cntr_scn_dt + hlf_zd_dc2, cntr_scn_dt - hlf_zd_dc2 : cntr_scn_dt + hlf_zd_dc2 ) = scn_dt( cntr_scn_dt - hlf_zd_dc2 : cntr_scn_dt + hlf_zd_dc2, cntr_scn_dt - hlf_zd_dc2 : cntr_scn_dt + hlf_zd_dc2 ) + zd_dc2 ;
    end
  return
  end

% Re-scaling, rotating, inclining a reference exozodiacal template.
% Exozodi file generated by B. Mennesson on 02/22/18 (zodipic,fnu,13,0.5,pixnum=128,zodis=3,inclination=60,positionangle=135,distance=5)
% It is 128x128 with 13 mas/pix. Units are Jy.
zd_bm = fitsread( [ opt.scn_dr 'ezflux.fits' ] ) ;
% Get the total flux for a 1 solar zodiacal dust emission at 1 pc and 500 nm (for a Sun-like star)
zd_flx_jy = sum( sum( zd_bm ) ) / 3 * 5^2 ;
% Get a clearly centered zodi
sz_zd_nw = ( 2 * floor( 13 / opt.px_scn_mas * 128 / 2 ) + 1 ) ;
zd_nw = imrotate( imresize( zd_bm, [ sz_zd_nw, sz_zd_nw ], 'bilinear' ), -135, 'bilinear', 'crop' ) ;
cntr_zd_nw = ( size( zd_nw, 1 ) + 1 ) / 2 ;
hlf_zd_nw = ( size( zd_nw, 1 ) - 1 ) / 2 ;
% Get a cut through the center
zd_rd = squeeze( zd_nw( cntr_zd_nw, cntr_zd_nw : end ) ) ;
% Get a characteristic size
d_zd = abs( zd_rd - max( zd_rd ) / 2 ) ;
% The AU have to be computed for the distance at which the zodi was simulated: 5 pc
opt.zd_org_hlf_au = find( d_zd == min( d_zd ) ) * opt.px_scn_mas / 200 ;
% Streching the radial profile of the zodi
fct_tmp = opt.zd.rd_au / opt.zd_org_hlf_au ;
n_zd_rd = numel( zd_rd ) ;
zd_rd = interp1( 1 : n_zd_rd, zd_rd, 1 : 1 / fct_tmp : n_zd_rd ) ; 
% Creating an ellipse with the radial zodi profile (just 0.03 sec, fast enough to avoid meshgrid, etc ...)
n_zd_rd = numel( zd_rd ) ;
% If the exozodi is going to be bigger than the scene, trim it
hlf_scn = ( size( scn_dt, 1 ) - 1 ) / 2 ;
fct_trim = 1 ;
  if n_zd_rd > hlf_scn
  n_zd_rd = hlf_scn ;
  % Estimating the part of the exozodi intensity that will be missing in the final image from the whole exozodi created
  fct_trim = sum( zd_rd( 1 : n_zd_rd ) ) / sum( zd_rd ) ;
  end
zd_2d = zeros( 2 * n_zd_rd + 1, 2 * n_zd_rd + 1 ) ;
% Exozodi inclination in rad
inc_rd = opt.zd.inc_dg * pi / 180 ;
  for i_1 = 1 : 2 * n_zd_rd + 1
    for i_2 = 1 : 2 * n_zd_rd + 1
    r_tmp = round( sqrt( ( i_1 - 1 - n_zd_rd )^2 + ( i_2 - 1 - n_zd_rd )^2 ) ) ;
    % Center left with contiguous value (fine because that's where the star will be )
      if r_tmp == 0
      r_tmp = 1 ;
      end
      if ( r_tmp <= n_zd_rd )
      idx_1 = round( n_zd_rd + 1 + ( i_1 - n_zd_rd - 1 ) * cos( inc_rd ) ) ; 
      zd_2d( idx_1, i_2 ) = zd_2d( idx_1, i_2 ) + zd_rd( r_tmp ) ;
      end
    end
  end

% Rotate to PA (https://en.wikipedia.org/wiki/Position_angle)
% Since the images are in the
%     |---E-->
%     N
%     |
%     v
% orientation, and PA goes from N to east, imrotate does the correct rotation 
zd_2d = imrotate( zd_2d, opt.zd.pa_dg, 'crop', 'bilinear' ) ;

% Normalize with zd_flx_jy and apply x solar-zodi and distance effect. Still at 500 nm for a Sun-like star
fct_recal = zd_flx_jy * opt.zd.fct / opt.str.dst_pc^2 * cos( inc_rd ) * fct_trim^2 ;
% Add the normalization wrt the star being used
fct_recal = fct_recal * opt.str.flx_rt_sun_500nm_jy ;
% Changing the units if necessary
  if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
  fct_recal = fct_recal * opt.conv_Jy_to_W_m2_um( i_lmbd ) ;
  end
zd_2d = zd_2d * fct_recal / sum( sum( zd_2d ) ) ;
% Adding the exozodi to the scene data
cntr_scn = ( size( scn_dt, 1 ) + 1 ) / 2 ;
hlf_zd_2d = n_zd_rd ;
idx_l_1 = cntr_scn - hlf_zd_2d ;
  if idx_l_1 < 1
  idx_l_1 = 1 ;
  end
idx_r_1 = cntr_scn + hlf_zd_2d ;
  if  idx_r_1 > size( scn_dt, 1 )
  idx_r_1 = size( scn_dt, 1 ) ;
  end
idx_l_2 = cntr_scn - hlf_zd_2d ;
  if idx_l_2 < 1
  idx_l_2 = 1 ;
  end
idx_r_2 = cntr_scn + hlf_zd_2d ;
  if  idx_r_2 > size( scn_dt, 2 ) ;
  idx_r_2 = size( scn_dt, 2 ) ; ;
  end
% Odd numbers
idx_l_1 = 2 * floor( idx_l_1 / 2 ) + 1 ;
idx_r_1 = 2 * floor( idx_r_1 / 2 ) + 1 ;
idx_l_2 = 2 * floor( idx_l_2 / 2 ) + 1 ;
idx_r_2 = 2 * floor( idx_r_2 / 2 ) + 1 ;
hlf_1 = ( idx_r_1 - idx_l_1 ) / 2 ;
hlf_2 = ( idx_r_2 - idx_l_2 ) / 2 ;
scn_dt( idx_l_1 : idx_r_1, idx_l_2 : idx_r_2 ) = scn_dt( idx_l_1 : idx_r_1, idx_l_2 : idx_r_2 ) + zd_2d( n_zd_rd + 1 - hlf_1 : n_zd_rd + 1 + hlf_1, n_zd_rd + 1 - hlf_2 : n_zd_rd + 1 + hlf_2 ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tweaking the Kuiper belt %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xkpr_blt = add_exokuiper_belt( opt, i_lmbd ) ; % Do not output opt from this function. It changes some values locally that should not be propagated.

% Opening the Haystacks cube with zero inclination. Units are Jy.
opt.scn.nm = 'modern_cube_zodi1inc0dist10' ;
scn_ext = read_haystacks_project( opt ) ;

% Removing the Sun's flux at the center
cntr_scn_ext = ( size( scn_ext, 1 ) + 1 ) / 2 ;
scn_ext( cntr_scn_ext, cntr_scn_ext ) = ( sum(sum( scn_ext( cntr_scn_ext - 2 : cntr_scn_ext + 2, cntr_scn_ext - 2 : cntr_scn_ext + 2 ) ) ) - scn_ext( cntr_scn_ext, cntr_scn_ext ) ) / 24 ;

% Haystacks pixel scale is 3 mas
px_hystcks_mas = 3 ;
% Haystacks image is at 10 pc
au2mas_hystcks = 1000 / 10 ;

% Extending the Haystacks data when using imrotate (some big enough value )
sz_scn_2 = 2 * floor( opt.scn.fov_diam_mas / px_hystcks_mas / 2 ) + 5 ; % pixels
cntr_scn_2 = ( sz_scn_2 + 1 ) / 2 ;
scn_2 = zeros( sz_scn_2, sz_scn_2 ) ;
hlf_scn = ( size( scn_ext, 1 ) - 1 ) / 2 ;
% If the FOV of the simulation is larger than the Haystacks FOV
  if ( ( cntr_scn_2 - hlf_scn ) > 0 )
  scn_2( cntr_scn_2 - hlf_scn : cntr_scn_2 + hlf_scn, cntr_scn_2 - hlf_scn : cntr_scn_2 + hlf_scn ) = scn_ext ;
  % Filling the new area with the edge value. This is fine, since it will be cropped later on.
  scn_2( 1 : cntr_scn_2 - hlf_scn - 1, : ) = scn_2( cntr_scn_2 - hlf_scn, cntr_scn_2 ) ;
  scn_2( cntr_scn_2 + hlf_scn + 1 : end, : ) = scn_2( cntr_scn_2 + hlf_scn, cntr_scn_2 ) ;
  scn_2( :, 1 : cntr_scn_2 - hlf_scn - 1 ) = scn_2( cntr_scn_2, cntr_scn_2 - hlf_scn ) ;
  scn_2( :, cntr_scn_2 + hlf_scn + 1 : end ) = scn_2( cntr_scn_2, cntr_scn_2 + hlf_scn ) ;
  % Redefining the variable
  scn_ext = scn_2 ;
  clear scn_2
  end

% Granularity
% pixels for the boxes
l_bx = 17 ;
% Assigning random values (Choosing the greater array to cover the case when the scene's FOV<Haystacks' FOV)
sz_kpr = max( [ sz_scn_2, size( scn_ext, 1 ) ] ) ;
cntr_xkpr = ( sz_kpr + 1 ) / 2 ;
sz_rndm = floor( sz_kpr / l_bx ) ;
% Using the same random realization
rng( 1 ) ;
scn_rndm = ( 0.5 + 0.5 * rand( sz_rndm, sz_rndm ) ) / 1.5 ; % Keeping the maximum value at 1
scn_rndm_fct = imresize( scn_rndm, ( sz_kpr + 0.5 ) / sz_rndm, 'bilinear' ) ;
% Some roundings add one pixel
  if ( size( scn_rndm_fct, 1 ) > sz_kpr )
  scn_rndm_fct = scn_rndm_fct( 2 : end, 2 : end ) ;
  end

xkpr_blt = scn_rndm_fct .* scn_ext ;
%clear scn_ext
% Constructing the Kuiper belt
x_ax = 1 : sz_kpr ;
[ x_px y_px ] = meshgrid( x_ax, x_ax ) ;
% Circular contours (original data with 0 inclination)
crc_xy = sqrt( ( x_px - cntr_xkpr ).^2 + ( ( y_px - cntr_xkpr ) ).^2 ) ;
% Building a radial-azimuthal scaling function
% The argument of the hyperbolic function changes depending on the disk location (azimuthal angle)
cntr_xkpr_px = ( opt.kpr.ext_au + opt.kpr.int_au ) / 2 * au2mas_hystcks / px_hystcks_mas ; % Final target size
wdth_kpr_px = ( opt.kpr.ext_au - opt.kpr.int_au ) / 2 * au2mas_hystcks / px_hystcks_mas ;
arg_hp = 3 * ( crc_xy - cntr_xkpr_px ) / wdth_kpr_px ;
% Some forward scattering
% Azimuthal angle (notice we take into account the ring is projected, so we remove the y cos factor to derive the true theta angle)
% The direction of the center of the foward scattering is zero degrees (N) and rotates as the PA parameter of the exozodiacal emission.
tht_xy = pi - flipud( abs( atan2( x_px - cntr_xkpr, y_px - cntr_xkpr ) ) ) ; % Centering angle=0 in the vertical direction, and then increasing towards pi in both directions
% Some amplitude for the main bulk of the exo-Kuiper emission (choosing a Gaussian profile)
sgm_tht = pi / 2.5 ;
% Forward scattering modifies the brighness between the relative factor of opt.kpr.mn_scttrng and opt.kpr.mx_scttrng.
cff_1 =  exp( pi^2 / 2 / sgm_tht^2 ) ;
cff_2 = ( cff_1 * opt.kpr.mn_scttrng - opt.kpr.mx_scttrng ) / ( cff_1 - 1 ) ;
cff_3 = cff_1 * ( opt.kpr.mx_scttrng - opt.kpr.mn_scttrng ) / ( cff_1 - 1 ) ;
fwd_scttrng = cff_2 + cff_3 * exp( -0.5 * ( tht_xy / sgm_tht ).^2 ) ;
% Additional intensity factor to apply to the Haystacks images (double tanh on top of each other)
fct_kpr = 1 - abs( tanh( 1.2 * arg_hp ) ).^(30)  + 1 - abs( tanh( 3 * arg_hp ) ).^(60);
% Normalized to 1
fct_kpr = fct_kpr / max( fct_kpr( : ) ) ;

% Additional factor to reduce the dust outside the belt (final target size)
fct_kpr_ext = 0 * fct_kpr +  1 ;
fct_kpr_ext( find( crc_xy > ( opt.kpr.ext_au * au2mas_hystcks / px_hystcks_mas ) ) ) = 0.2 ;
% Similarly for the interior
fct_kpr_int = 0 * fct_kpr +  1 ;
fct_kpr_int( find( crc_xy < ( opt.kpr.int_au * au2mas_hystcks / px_hystcks_mas ) ) ) = 0.2 ;
% Combined
fct_kpr_edg = fct_kpr_int .* fct_kpr_ext ;

% Multiplying by the factors for forward scattering and
xkpr_blt = xkpr_blt .* fct_kpr .* fct_kpr_edg .* fwd_scttrng ;

% Inclining the exo-Kuiper disk: elliptical contours
xkpr_blt_tmp = 0 * xkpr_blt ;
hlf_kpr = ( sz_kpr + 1 ) / 2 ; % Must be an integer
% This loop takes ~0.16 seconds (<<full simulation)
  for i_x = 1 : sz_kpr
  xkpr_blt_tmp( hlf_kpr + round( ( i_x - hlf_kpr ) * cos( pi / 180 * opt.zd.inc_dg ) ), : ) = xkpr_blt( i_x, : ) ;
  end
clear xkpr_blt
xkpr_blt = xkpr_blt_tmp ;
clear xkpr_blt_tmp

% Rotating it to the final configuration
xkpr_blt = imrotate( xkpr_blt, 180 + opt.zd.pa_dg, 'crop', 'bilinear' ) ;
% Final intensity (same relative factor wrt the Solar System as the interior exozodiacal emission)
xkpr_blt = xkpr_blt * opt.zd.fct ;

% Changing the distance (we will impose that the pixel scale is the same before/after the distance change, so that the brightness density per pixel is the same. Therefore, imresize, does not need any correction factor). The Haystacks data are at 10 pc
sz_dstnc = 2 * floor( size( xkpr_blt, 1 ) * 10 / opt.str.dst_pc / 2 ) + 1 ;
xkpr_blt = imresize( xkpr_blt, [ sz_dstnc sz_dstnc ] ) ;
% Changing the units if necessary
  if strcmp( lower( opt.scn.unts ), 'w/m2/um' )
  xkpr_blt = xkpr_blt * opt.conv_Jy_to_W_m2_um( i_lmbd ) ;
  end
% Cutting-out the FOV
hlf_fov_px = ( opt.n_px_scn - 1 ) / 2 ;
cntr_xkpr_blt = ( size( xkpr_blt, 1 ) + 1 ) / 2 ;
  if ( cntr_xkpr_blt - hlf_fov_px > 0 )
  xkpr_blt_tmp = xkpr_blt( cntr_xkpr_blt - hlf_fov_px : cntr_xkpr_blt + hlf_fov_px, cntr_xkpr_blt - hlf_fov_px : cntr_xkpr_blt + hlf_fov_px ) ;
  clear xkpr_blt
  xkpr_blt = xkpr_blt_tmp ;
  clear xkpr_blt_tmp
  else
  disp( 'The construction of the exo-Kuiper belt requires a size for the input scene that does not seem to be correct. Stopped.' )
  make_a_stop
  end

% Here it shows it is 499.449 nm
% fts_info = fitsinfo( sprintf( '%sinput_scenes/modern_cube_zodi1inc60dist10_0.46-0.57um.fits', opt.rt_dr ) ) ;
% lmbd_tmp = fts_info.Image( 25 ).Keywords( 9, 2 ) ;

% Plot that shows Neptune and the roll-off of the Kuiper dust (not dramatic)
% plot( ((2000-1667):(3333-1667))*3/100, scn_ext(1673,2000:end)) ; xlabel( 'AU' ) ; ylabel( 'Jy' ) ; title( '         Horizontal cut Haystacks solar model 500 nm (10 pc, 60 deg inc)', 'FontSize', 10 ) ; ylim( [ 1.85e-11, 2.15e-11 ] ) ;
% text( 32, 2.1e-11, sprintf( 'Neptune=%2.2d Jy', max( scn_ext( 1673, 2000:end) ) ) )
% text( 32, 2.08e-11, sprintf( 'Sun=%2.2d Jy', scn_ext(1667,1667) ) )

% Constructing the PSF from the Starshade
% Written by S.B. Shaklan. Adapted to SISTER by S.R. Hildebrandt to consider any telescope size.
function psf_ss = get_psf_from_starshade( opt ) 
% parabolic sag in meters.  1 um for remote occulter. For 39 m aperture, this is just 1 um sag.
sag = 0.5 * ( opt.dmtr_tlscp_m / 2 )^2 / opt.dst_strshd_tlscp_m ;  
% Conversion from mas (seen at the telescope) to m (on the starshade) 
mas2m = 1 / 3600 / 1e3 / 180 * pi * opt.dst_strshd_tlscp_m ;
% PSF needs to have enough pixels for the calculations. 
% r_out is the radius of the aperture in some new pixels
r_px_out = 20 ; % Any value to sample appropriately the primary and secondary apertures
% New pixel size
nw_px_mas = opt.dmtr_tlscp_m / mas2m / r_px_out ;
% Ratio wrt the actual pixels sampling the PSF in SISTER (the pxiel scale of the PSF used in the convolution witht he scene is opt.px_scn_mas)
px_rt = nw_px_mas / opt.px_scn_mas ;
% Size of a new array that is similar in absolute size to the current PSF (odd number of elements)
sz_pxpsf = 2 * floor( round( opt.psf_arry_px( 1 ) / px_rt ) / 2 ) + 1 ;
[ x, y ] = meshgrid( 1 : sz_pxpsf, 1 : sz_pxpsf ) ;
x = x - ( x( end ) - 1 ) / 2 ;
y = y - ( y( end ) - 1 ) / 2 ;  
r_px = sqrt( x.^2 + y.^2 ) ;
% meters of sag.  Goes to value sag at edge of aperture
sagval = sag * ( r_px * nw_px_mas * mas2m / opt.dmtr_tlscp_m ).^2 ; 
sagfun = exp( -2 * pi * i * sagval / ( opt.lmbd_tmp_nm * 1e-9 ) ) ;

apfun = 0 * sagfun ;
% Primary (r_px_out) & secondary (r_px_in)
r_px_in = r_px_out * opt.secondary_size ;
% Primary
apfun( r_px <= r_px_out ) = 1 ;
% Secondary
apfun( r_px <= r_px_in ) = 0 ;
apfun_sag = apfun .* sagfun ;
% Field 
psf_fld = fftshift( fft2( apfun_sag ) ) ;
% Intensity
psf_ss = abs( psf_fld ).^2 ;
% The fft2 changes the meaning of the scale of the pixels, so we compare with an ideal aperture
apfun( r_px <= r_px_out ) = 1 ;
% Field for the ideal aperture
psf_fld = fftshift( fft2( apfun ) ) ;
% Intensity
psf_ideal = abs( psf_fld ).^2 ;
% Estimating the FWHM
cntr_psf = ( size( psf_ideal, 1 ) + 1 ) / 2 ;
psf_1d = squeeze( psf_ss( cntr_psf, : ) ) ; 
psf_1d = psf_1d / max( psf_1d ) ; 
% Increasing the smaple points
idx_tmp = 1 : 1/1000 : numel( psf_1d ) ; 
psf_1d_intrp = interp1( 1 : numel( psf_1d ), psf_1d, idx_tmp ) ; 
d_fwhm = abs( psf_1d_intrp - 0.5 ); 
n_idx = numel( idx_tmp ) ;
q = find( abs( d_fwhm( 1 : round( n_idx / 2 ) ) ) == min( d_fwhm( 1 : round( n_idx / 2 ) ) ) ) ;
q2 =find( abs( d_fwhm( round( n_idx / 2 ) : end ) ) == min( d_fwhm( round( n_idx / 2 :end ) ) ) ) ;
fwhm_ideal_mas = ( round( n_idx / 2 ) + q2 - 1 - q ) / 1000 * opt.px_scn_mas ;
% FWHM of an Airy PSF: 1.02 l/D
fwhm_airy_mas = 1.02 * opt.lmbd_tmp_nm * 1e-9 / opt.dmtr_tlscp_m * 180 * 3600e3 / pi ;
% Re-size to match both ideal cases (odd number of pixels). The final pixel scale is opt.px_scn_mas
sz_end = 2 * floor( fwhm_airy_mas / fwhm_ideal_mas * size( psf_ss, 1 ) / 2 ) + 1 ;
psf_ss = imresize( psf_ss, [ sz_end, sz_end ], 'bilinear' ) ;

% Padding or triming to match the PSF used to convolve in SISTER
  if ( size( psf_ss, 1 ) < opt.psf_arry_px( 1 ) )
  psf_ss_tmp = zeros( opt.psf_arry_px ) ;
  cntr_psf_sstr = ( ( opt.psf_arry_px( 1 ) + 1 ) / 2 ) ;
  hlf_psf_ss = ( size( psf_ss, 1 ) - 1 ) / 2 ;
  psf_ss_tmp( cntr_psf_sstr - hlf_psf_ss : cntr_psf_sstr + hlf_psf_ss, cntr_psf_sstr - hlf_psf_ss : cntr_psf_sstr + hlf_psf_ss ) = psf_ss ;
  psf_ss = psf_ss_tmp ; clear psf_ss_tmp
  else
  cntr_psf_ss = ( size( psf_ss, 1 ) + 1 ) / 2 ;
  hlf_psf_sstr = ( opt.psf_arry_px( 1 ) - 1 ) / 2 ;
  psf_ss = psf_ss( cntr_psf_ss - hlf_psf_sstr : cntr_psf_ss + hlf_psf_sstr, cntr_psf_ss - hlf_psf_sstr : cntr_psf_ss + hlf_psf_sstr ) ;
  end

% Normalize after padding (in the case is trimed down, it will take into account the intensity loss)
psf_ss = psf_ss / sum( psf_ss( : ) ) ;

% Planet's PSF for a spinning Starshade
function [ plnt_psf opt ] = get_planet_spinning_psf( opt, psf_array, scn_dt_cnv, i_pl, i_lmbd )

% The stationary PSF
psf_st = squeeze( psf_array( end, :, : ) ) ;
sz_st_1 = size( psf_st, 1 ) ;
sz_st_2 = size( psf_st, 2 ) ;
hlf_st_1 = ( sz_st_1 - 1 ) / 2 ;
hlf_st_2 = ( sz_st_2 - 1 ) / 2 ;
cntr_st_1 = ( sz_st_1 + 1 ) / 2 ;
cntr_st_2 = ( sz_st_2 + 1 ) / 2 ;

% Center of the scene
cntr_scn_dt_1 = ( size( scn_dt_cnv, 1 ) + 1 ) / 2 ;
cntr_scn_dt_2 = ( size( scn_dt_cnv, 2 ) + 1 ) / 2 ;
% Half width of the scene
hlf_scn_dt_1 = ( size( scn_dt_cnv, 1 ) - 1 ) / 2 ;
hlf_scn_dt_2 = ( size( scn_dt_cnv, 2 ) - 1 ) / 2 ;

% Choosing the corresponding PSF
pl_1_px = opt.plnt.add.arc_dc_px( i_pl ) ;
% Lateral shift if any (cubes derive the planet's position after reverting the shift. pl_1_px only defines the planet's position here to identify the corresponding PSF)
  if ~isnan( opt.scn.ra_shft_px )
  pl_1_px = pl_1_px + round( opt.scn.dc_shft_px ) ; % Planets are point-like
  end

% Some planets may be outside the FOV
fov_add_1 = 0 ;
  if abs( pl_1_px ) > hlf_scn_dt_1
  disp( sprintf( 'Notice: planet %i is outside the FOV', i_pl ) )
  fov_add_1 = ( abs( pl_1_px ) - hlf_scn_dt_1 ) * opt.px_scn_mas ;
  end
pl_2_px = opt.plnt.add.arc_ra_px( i_pl ) ;
  if ~isnan( opt.scn.ra_shft_px )
  pl_2_px = pl_2_px + round( opt.scn.ra_shft_px ) ; % Planets are point-like
  end

fov_add_2 = 0 ;
  if abs( pl_2_px ) > hlf_scn_dt_2
  disp( sprintf( 'Notice: planet %i is outside the FOV', i_pl ) )
  fov_add_2 = ( abs( pl_2_px ) - hlf_scn_dt_2 ) * opt.px_scn_mas ;
  end
  if ( fov_add_1 ) || ( fov_add_2 )
  disp( sprintf( 'In order to include planet %i the diameter of the FOV should at least be %4.0f mas', i_pl, round( opt.scn.fov_diam_mas + 2 * max( [ fov_add_1, fov_add_2 ] ) ) ) )
  plnt_psf = NaN ;
  return
  end
r_pl_px = sqrt( pl_1_px^2 + pl_2_px^2 ) ;
  % If it is in the stationary region, choose the last element of the PSF array and recall it is assumed to be translationally invariant. Otherwise, perform the corresponding rotation.
  if ( r_pl_px >= ( opt.r_st_mas / opt.px_scn_mas ) )
  psf_rtd = psf_st ;
  else
  % Getting an accurate position in mas (the non-stationary basis is computed every mas)
  r_psf = round( sqrt( opt.plnt.add.arc_dc_mas( i_pl )^2 + opt.plnt.add.arc_ra_mas( i_pl )^2 ) ) ;
  % Some rounding may get it a mas farther, but it would be the same PSF (stationary one)
    if ( r_psf + 1 > size( psf_array, 1 ) )
    r_psf = size( psf_array, 1 ) - 1 ;
    end
  psf_non_st = squeeze( psf_array( r_psf + 1, :, : ) ) ;
  % Rotating the non-stationary PSF to the position of the planet
  % PSF rotation towards the 'y' axis (checked it rotates the PSF to the position of the planet (pl_1_px,pl_2_px) and that the asymmetric structure of the PSF is consistent with the corresponding non-spinning PSF at that location.
  psf_rtd = imrotate( psf_non_st, - 180 / pi * atan2( opt.plnt.add.arc_dc_mas( i_pl ), opt.plnt.add.arc_ra_mas( i_pl ) ), 'crop', 'bilinear' ) ;
  end

% In some cases (small FOV) part of the planet's PSF may be off the scene
idx_pl_1_1 = cntr_scn_dt_1 + pl_1_px - hlf_st_1 ;
idx_pl_1_2 = cntr_scn_dt_1 + pl_1_px + hlf_st_1 ;
idx_pl_2_1 = cntr_scn_dt_2 + pl_2_px - hlf_st_2 ;
idx_pl_2_2 = cntr_scn_dt_2 + pl_2_px + hlf_st_2 ;
l_1_1 = 0 ; l_1_2 = 0 ; l_2_1 = 0 ; l_2_2 = 0 ;
sz_1 = size( scn_dt_cnv, 1 ) ;
sz_2 = size( scn_dt_cnv, 2 ) ;
  if idx_pl_1_1 < 1
  l_1_1 = 1 - idx_pl_1_1 ;
  idx_pl_1_1 = 1 ;
  end
  if idx_pl_1_2 > sz_1
  l_1_2 = idx_pl_1_2 - sz_1 ;
  idx_pl_1_2 = sz_1 ;
  end
  if idx_pl_2_1 < 1
  l_2_1 = 1 - idx_pl_2_1 ;
  idx_pl_2_1 = 1 ;
  end
  if idx_pl_2_2 > sz_2
  l_2_2 = idx_pl_2_2 - sz_2 ;
  idx_pl_2_2 = sz_2 ;
  end
% PSF, centered on the planet and cropped to the scene FOV 
plnt_psf = psf_rtd( l_1_1 + 1 : end - l_1_2, l_2_1 + 1 : end - l_2_2 ) ;
opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) = idx_pl_1_1 ; opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ) = idx_pl_1_2 ; opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) = idx_pl_2_1 ; opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) = idx_pl_2_2 ;
opt.cb.plnt.l_1_1( i_pl, i_lmbd ) = l_1_1 ; opt.cb.plnt.l_1_2( i_pl, i_lmbd ) = l_1_2 ; opt.cb.plnt.l_2_1( i_pl, i_lmbd ) = l_2_1 ; opt.cb.plnt.l_2_2( i_pl, i_lmbd ) = l_2_2 ;
% Saving some memory. In most applications, 7 digit precision is enough
  if ( opt.sngl_prcsn )
  plnt_psf = single( plnt_psf ) ;
  end

% Planet's PSF for a non-spinning Starshade
function [ plnt_psf opt ] = get_planet_non_spinning_psf( opt, psf_st, scn_dt_cnv, i_pl, i_lmbd, i_bnd )

% The stationary PSF
sz_st_1 = size( psf_st, 1 ) ;
sz_st_2 = size( psf_st, 2 ) ;
hlf_st_1 = ( sz_st_1 - 1 ) / 2 ;
hlf_st_2 = ( sz_st_2 - 1 ) / 2 ;
cntr_st_1 = ( sz_st_1 + 1 ) / 2 ;
cntr_st_2 = ( sz_st_2 + 1 ) / 2 ;

% Center of the scene
cntr_scn_dt_1 = ( size( scn_dt_cnv, 1 ) + 1 ) / 2 ;
cntr_scn_dt_2 = ( size( scn_dt_cnv, 2 ) + 1 ) / 2 ;
% Half width of the scene
hlf_scn_dt_1 = ( size( scn_dt_cnv, 1 ) - 1 ) / 2 ;
hlf_scn_dt_2 = ( size( scn_dt_cnv, 2 ) - 1 ) / 2 ;

% Choosing the corresponding PSF
pl_1_px = opt.plnt.add.arc_dc_px( i_pl ) ;
% Some planets may be outside the FOV
fov_add_1 = 0 ;
  if pl_1_px > hlf_scn_dt_1
  disp( sprintf( 'Notice: planet %i is outside the FOV', i_pl ) )
  fov_add_1 = ( pl_1_px - hlf_scn_dt_1 ) * opt.px_scn_mas ;
  end
pl_2_px = opt.plnt.add.arc_ra_px( i_pl ) ;
fov_add_2 = 0 ;
  if pl_2_px > hlf_scn_dt_2
  disp( sprintf( 'Notice: planet %i is outside the FOV', i_pl ) )
  fov_add_2 = ( pl_2_px - hlf_scn_dt_2 ) * opt.px_scn_mas ;
  end
  if ( fov_add_1 ) || ( fov_add_2 )
  disp( sprintf( 'In order to include planet %i the diameter of the FOV should be at least %4.0f mas', i_pl, round( opt.scn.fov_diam_mas + 2 * max( [ fov_add_1, fov_add_2 ] ) ) ) )
  plnt_psf = NaN ;
  return
  end
r_pl_px = sqrt( pl_1_px^2 + pl_2_px^2 ) ;
  % If it is in the stationary region, use the stationary PSF. Otherwise, build the corresponding PSF
  if ( r_pl_px >= ( opt.r_st_mas / opt.px_scn_mas ) )
  psf_tmp = psf_st ;
  else
  % Notice the x/y vs. RA/DEC order for the PSF basis
  psf_tmp = get_normalized_psf( i_bnd, opt, round( opt.plnt.add.arc_dc_mas( i_pl ) ), round( opt.plnt.add.arc_ra_mas( i_pl ) ) ) ;
  end

% In some cases (small FOV) part of the planet's PSF may be off the scene
idx_pl_1_1 = cntr_scn_dt_1 + pl_1_px - hlf_st_1 ;
idx_pl_1_2 = cntr_scn_dt_1 + pl_1_px + hlf_st_1 ;
idx_pl_2_1 = cntr_scn_dt_2 + pl_2_px - hlf_st_2 ;
idx_pl_2_2 = cntr_scn_dt_2 + pl_2_px + hlf_st_2 ;
l_1_1 = 0 ; l_1_2 = 0 ; l_2_1 = 0 ; l_2_2 = 0 ;
sz_1 = size( scn_dt_cnv, 1 ) ;
sz_2 = size( scn_dt_cnv, 2 ) ;
  if idx_pl_1_1 < 1
  l_1_1 = 1 - idx_pl_1_1 ;
  idx_pl_1_1 = 1 ;
  end
  if idx_pl_1_2 > sz_1
  l_1_2 = idx_pl_1_2 - sz_1 ;
  idx_pl_1_2 = sz_1 ;
  end
  if idx_pl_2_1 < 1
  l_2_1 = 1 - idx_pl_2_1 ;
  idx_pl_2_1 = 1 ;
  end
  if idx_pl_2_2 > sz_2
  l_2_2 = idx_pl_2_2 - sz_2 ;
  idx_pl_2_2 = sz_2 ;
  end
plnt_psf = psf_tmp( l_1_1 + 1 : end - l_1_2, l_2_1 + 1 : end - l_2_2 ) ;
opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) = idx_pl_1_1 ; opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ) = idx_pl_1_2 ; opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) = idx_pl_2_1 ; opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) = idx_pl_2_2 ;
opt.cb.plnt.l_1_1( i_pl, i_lmbd ) = l_1_1 ; opt.cb.plnt.l_1_2( i_pl, i_lmbd ) = l_1_2 ; opt.cb.plnt.l_2_1( i_pl, i_lmbd ) = l_2_1 ; opt.cb.plnt.l_2_2( i_pl, i_lmbd ) = l_2_2 ;

% Saving some memory. In most applications, 7 digit precision is enough
  if ( opt.sngl_prcsn )
  plnt_psf = single( plnt_psf ) ;
  end

% Sub-function to write a file with the data per wavelength slice
function fl_nm_cb = write_starshade_cube_file( scn_tmp, i_lmbd, opt, i_bck )

% Ouput location
dr_fts = [ opt.out_dr opt.cb.dr '/' ] ;
  if ~isdir( dr_fts )
  mkdir( dr_fts ) ;
  end
% Name of the output file
opt_cb_lbl = '' ;
  if ~strcmp( opt.cb.lbl, '' )
  opt_cb_lbl = [ '_' opt.cb.lbl ] ;
  end
  if ~( opt.cb.bckgrnd( i_bck ) )
  opt_cb_lbl = [ opt_cb_lbl '_no_background' ] ;
  end
  if ( opt.n_lp > 1 )
  opt_cb_lbl = sprintf( '%s_%i', opt_cb_lbl, opt.i_lp ) ;
  end

fl_nm_cb = sprintf( '%ssister_cube_%s_%s_mas_%s_mas_%i%s', dr_fts, opt.add_label, strrep( sprintf( '%.1f', opt.cb.ra_mas ), '-', 'm' ), strrep( sprintf( '%.1f', opt.cb.dc_mas ), '-', 'm' ), numel( opt.lmbd_arry_img_nm ), opt_cb_lbl ) ;
% If the cubes are extracted about each planet
  if ( opt.cb.plnt.do )
  fl_nm_cb = sprintf( '%ssister_cube_%s_%04i_%04i_%i_planets%s', dr_fts, opt.add_label, opt.lmbd_img_1_nm, opt.lmbd_img_2_nm, numel( opt.lmbd_arry_img_nm ), opt_cb_lbl ) ;
  end
% Storing the filename
opt.cb.fl_nm = fl_nm_cb ;
% For the FITS file data
fl_nm_cb_psf = [ fl_nm_cb '_psf' ] ;

% Erase anything present at the start
  if ( i_lmbd == 1 ) && ( opt.cb.strt )
  fl_cb = dir( [ opt.cb.fl_nm '*' ] ) ;
    if ( numel( fl_cb ) )
    delete( [ opt.cb.fl_nm '*' ] )
    end
  end

% Number of wavelength bins under consideration
n_lmbd = numel( opt.lmbd_arry_img_nm ) ;

% Getting the sub-array of the scene to be stored (recall that an array is ordered as DECxRA in Matlab)

% If the cut is for some specific region (not around each planet)
  if ~( opt.cb.plnt.do )
  sz_tmp_1 = size( scn_tmp, 1 ) ;
  cntr_scn_1 = ( sz_tmp_1 - 1 ) / 2 + 1 ;
  hlf_scn_1 = cntr_scn_1 - 1 ;
  % Proposed margins
  idx_dc_1 = floor( cntr_scn_1 + ( opt.cb.dc_mas  - opt.cb.fov_diam_mas / 2 ) / opt.px_cmr_mas ) ;
    if idx_dc_1 < 1
    idx_dc_1 = 1 ;
    end
  idx_dc_2 = floor( cntr_scn_1 + ( opt.cb.dc_mas  + opt.cb.fov_diam_mas / 2 ) / opt.px_cmr_mas ) ;
    if idx_dc_2 > sz_tmp_1
    idx_dc_2 = sz_tmp_1 ;
    end

  idx_ra_1 = floor( cntr_scn_1 + ( opt.cb.ra_mas  - opt.cb.fov_diam_mas / 2 ) / opt.px_cmr_mas ) ;
    if idx_ra_1 < 1
    idx_ra_1 = 1 ;
    end
  idx_ra_2 = floor( cntr_scn_1 + ( opt.cb.ra_mas  + opt.cb.fov_diam_mas / 2 ) / opt.px_cmr_mas ) ;
    if idx_ra_2 > sz_tmp_1
    idx_ra_2 = sz_tmp_1 ;
    end
  % Reading the existing data
    if ( i_lmbd == 1 )
    sister_cube = scn_tmp( idx_dc_1 : idx_dc_2, idx_ra_1 : idx_ra_2 ) ;
    else
      if ( opt.cb.fits )
      sister_cube = fitsread( [ fl_nm_cb '.fits' ] ) ;
      else
      load( [ fl_nm_cb '.mat' ] )
      end
      if ( i_lmbd == 2 )
      sister_cube_2( 1, :, : ) = sister_cube ;
      sister_cube_2( 2, :, : ) = scn_tmp( idx_dc_1 : idx_dc_2, idx_ra_1 : idx_ra_2 ) ;
      clear sister_cube
      sister_cube = sister_cube_2 ;
      clear sister_cube_2
      else
      sister_cube( i_lmbd, :, : ) = scn_tmp( idx_dc_1 : idx_dc_2, idx_ra_1 : idx_ra_2 ) ;
      end
    end
    % Storing: FITS or Matlab
    if ( opt.cb.fits )
    fitswrite( sister_cube, [ fl_nm_cb '.fits' ] ) ;
      if ( i_lmbd ==1 ) || ( i_lmbd == n_lmbd )
      disp( sprintf( 'Cube data file %s.fits written', fl_nm_cb ) )
      end
    else
    opt_sister = opt ;
    save( [ fl_nm_cb '.mat' ], 'sister_cube', 'opt_sister' ) ;
      if ( i_lmbd ==1 ) || ( i_lmbd == n_lmbd )
      disp( sprintf( 'Cube data file %s.mat written', fl_nm_cb ) )
      end
    end
  end

% If the data are extracted around each planet
  if ( opt.cb.plnt.do )
  sz_tmp_1 = size( scn_tmp, 1 ) ;
  cntr_scn_1 = ( sz_tmp_1 - 1 ) / 2 + 1 ;
  hlf_scn_1 = cntr_scn_1 - 1 ;
  % The data cube will be centered at the planet location +/- the extension, unless some planet has part of the cube area falling off the full FOV
  opt.cube.fov_diam_mas = opt.cb.plnt.n_fwhm * opt.fwhm_mas( i_lmbd, 2 ) ; % The first element is the wavelength
  opt.cb.fov_diam_mas = opt.cube.fov_diam_mas ;
  hlf_cb_plnt_px = ceil( opt.cb.fov_diam_mas / 2 / opt.px_cmr_mas ) ; % Rounded up
    for i_pl = 1 : opt.n_pl
    % Approximate pixel where the planet should lay on. Trying to imitate what imresize may do when changing from scene pixels to camera pixels.
    plnt_dc_px = sign( opt.plnt.add.arc_dc_mas( i_pl ) ) * ceil( abs( opt.plnt.add.arc_dc_mas( i_pl ) ) / opt.px_cmr_mas ) ;
    plnt_ra_px = sign( opt.plnt.add.arc_ra_mas( i_pl ) ) * ceil( abs( opt.plnt.add.arc_ra_mas( i_pl ) ) / opt.px_cmr_mas ) ;
    idx_dc_1 = cntr_scn_1 + plnt_dc_px - hlf_cb_plnt_px ;
      if idx_dc_1 < 1
      disp( sprintf( 'WARNING: part of the area of planet #%i cube lies out off the FOV. Make the FOV larger by %i mas, or reduce the area around the planets. Stopped.', i_pl, ceil( ( 1 - idx_dc_1 ) * opt.px_cmr_mas ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    idx_dc_2 = cntr_scn_1 + plnt_dc_px + hlf_cb_plnt_px ;
      if idx_dc_2 > sz_tmp_1
      disp( sprintf( 'WARNING: part of the area of the planet #%i cube lies out off the FOV. Make the FOV larger by %i mas, or reduce the area around the planets. Stopped.', i_pl, ceil( ( idx_dc_2 - sz_tmp_1 ) * opt.px_cmr_mas ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end

    idx_ra_1 = cntr_scn_1 + plnt_ra_px - hlf_cb_plnt_px ;
      if idx_ra_1 < 1
      disp( sprintf( 'WARNING: part of the area of the planet #%i cube lies out off the FOV. Make the FOV larger by %i mas, or reduce the area around the planets. Stopped.', i_pl, ceil( ( 1 - idx_ra_1 ) * opt.px_cmr_mas ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    idx_ra_2 = cntr_scn_1 + plnt_ra_px + hlf_cb_plnt_px ;
      if idx_ra_2 > sz_tmp_1
       disp( sprintf( 'WARNING: part of the area of the planet #%i cube lies out off the FOV. Make the FOV larger by %i mas, or reduce the area around the planets. Stopped.', i_pl, ceil( ( idx_ra_2 - sz_tmp_1 ) * opt.px_cmr_mas ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    % Keeping track of the camera's pixel indices used to extract the cubes
    opt.cb.plnt.idx_dc_1( i_pl, i_lmbd ) = idx_dc_1 ; opt.cb.plnt.idx_dc_2( i_pl, i_lmbd ) = idx_dc_2 ;
    opt.cb.plnt.idx_ra_1( i_pl, i_lmbd ) = idx_ra_1 ; opt.cb.plnt.idx_ra_2( i_pl, i_lmbd ) = idx_ra_2 ;
    sister_cube_tmp( i_pl, 1 : 2 * hlf_cb_plnt_px + 1, 1 : 2 * hlf_cb_plnt_px + 1 ) = scn_tmp( idx_dc_1 : idx_dc_2, idx_ra_1 : idx_ra_2 ) ;
    % Building the corresponding resized PSF for each particular planet
    % Same size as the scene
    scn_psf = zeros( opt.n_px_scn, opt.n_px_scn ) ;
    % Saving some memory. In most applications, 7 digit precision is enough
      if ( opt.sngl_prcsn )
      scn_psf = single( scn_psf ) ;
      end
    % Inserting the corresponding PSF for the planet (NB: same way as planets were added in add_planets_spinning/non_spinning)
    scn_psf( opt.cb.plnt.idx_pl_1_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_1_2( i_pl, i_lmbd ), opt.cb.plnt.idx_pl_2_1( i_pl, i_lmbd ) : opt.cb.plnt.idx_pl_2_2( i_pl, i_lmbd ) ) = opt.cb.plnt.plnt_psf( i_pl, i_lmbd, 1 : opt.cb.plnt.n_plnt_psf_1( i_pl, i_lmbd ) , 1 : opt.cb.plnt.n_plnt_psf_2( i_pl, i_lmbd ) ) ;
    scn_psf = imresize( scn_psf, [ opt.px_cmr_fov opt.px_cmr_fov ], 'bilinear' ) ;
    % imresize needs to be rescaled: bilinear averages
    scn_psf = scn_psf * ( opt.n_px_scn / opt.px_cmr_fov )^2 ;
    % Adding some jitter motion of the telescope, if specified. Modeled as a 2-D Gaussian with sigma=RMS
      if ( opt.jttr.do )
        if ( i_pl == 1 )
        disp( sprintf( 'Applying a jitter motion of RMS %2.2f mas', opt.jttr.rms_mas ) )
        end
      scn_psf = imgaussfilt( scn_psf, opt.jttr.rms_mas / opt.px_scn_mas ) ;
      end
    psf_cube_tmp( i_pl, 1 : 2 * hlf_cb_plnt_px + 1, 1 : 2 * hlf_cb_plnt_px + 1 ) = scn_psf( idx_dc_1 : idx_dc_2, idx_ra_1 : idx_ra_2 ) ;
    % I/O results
      if ( i_lmbd == 1 )
      sister_cube = sister_cube_tmp ;
      psf_cube = psf_cube_tmp ;
      else
      % Load the results first
        if ( i_pl == 1 )
          if ( opt.cb.fits )
          sister_cube = fitsread( [ fl_nm_cb '.fits' ] ) ;
          psf_cube = fitsread( [ fl_nm_cb_psf '.fits' ] ) ;
          else
          load( [ fl_nm_cb '.mat' ] )
          end
        % For i_lmbd==2, the new dimension is created, so the i_lmbd=1 gets promoted to (1,n_pl,:,:)
          if ( i_lmbd == 2 )
          sister_cube_2( 1, :, :, : ) = sister_cube ;
          sister_cube = sister_cube_2 ; clear sister_cube_2 ;
          psf_cube_2( 1, :, :, : ) = psf_cube ;
          psf_cube = psf_cube_2 ; clear psf_cube_2
          end
        % Re-center the data with the new data cube that will be appended
        sz_cb_tmp = size( sister_cube, 3 ) ; % has to be a square with an odd number of elements
        hlf_cb_tmp = ( sz_cb_tmp - 1 ) / 2 ;
        sister_cube_2( 1 : i_lmbd - 1, 1 : size( sister_cube, 2 ), hlf_cb_plnt_px + 1 - hlf_cb_tmp : hlf_cb_plnt_px + 1 + hlf_cb_tmp, hlf_cb_plnt_px + 1 - hlf_cb_tmp : hlf_cb_plnt_px + 1 + hlf_cb_tmp ) = sister_cube ;
        clear sister_cube
        sister_cube = sister_cube_2 ;
        clear sister_cube_2
        psf_cube_2( 1 : i_lmbd - 1, 1 : size( sister_cube, 2 ), hlf_cb_plnt_px + 1 - hlf_cb_tmp : hlf_cb_plnt_px + 1 + hlf_cb_tmp, hlf_cb_plnt_px + 1 - hlf_cb_tmp : hlf_cb_plnt_px + 1 + hlf_cb_tmp ) = psf_cube ;
        clear psf_cube
        psf_cube = psf_cube_2 ;
        clear psf_cube_2
        end % i_pl==1
      sister_cube( i_lmbd, i_pl, 1 : 2 * hlf_cb_plnt_px + 1, 1 :  2 * hlf_cb_plnt_px + 1 ) = sister_cube_tmp( i_pl, :, : ) ;
      psf_cube( i_lmbd, i_pl, 1 : 2 * hlf_cb_plnt_px + 1, 1 :  2 * hlf_cb_plnt_px + 1 ) = psf_cube_tmp( i_pl, :, : ) ;
      end % i_lmbd>1
    end % i_pl

    % Storing: FITS or Matlab
    if ( opt.cb.fits )
    fitswrite( sister_cube, [ fl_nm_cb '.fits' ] ) ;
    fitswrite( psf_cube, [ fl_nm_cb_psf '.fits' ] ) ;
      if ( i_lmbd ==1 ) || ( i_lmbd == n_lmbd )
      disp( sprintf( 'Cube data file %s.fits written', fl_nm_cb ) )
      end
    else
      if ( i_lmbd < n_lmbd )
      save( [ fl_nm_cb '.mat' ], 'sister_cube', 'psf_cube' ) ;
      else
      opt_sister = opt ;
      clear opt
      % No need to store replicated information
      opt_sister.cb.plnt = rmfield( opt_sister.cb.plnt, 'plnt_psf' ) ;
      save( [ fl_nm_cb '.mat' ], 'sister_cube', 'psf_cube', 'opt_sister' ) ;
      end
      if ( i_lmbd ==1 ) || ( i_lmbd == n_lmbd )
      disp( sprintf( 'Cube data file %s.mat written', fl_nm_cb ) )
      end
    end % Storing FITS or Matlab
  end % opt.cube.planets.do
