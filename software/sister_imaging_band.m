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
function [ scn_dt_cnv_e, scn_dt_cnv_no_plnt_e, scn_ns_e, scn_dt_cnv_lcl_zd_e, opt_img ] = sister_imaging_band( i_bnd, opt_img )
% Script to generate a wavelength slice simulation of SISTER
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Debug mode yes/no (exit back to regular session mode by typing 'dbquit all' if it crashes)
dbstop if error

% Looping over different changes of the main setup (e.g., Keplerian orbits)
% Unless set earlier the loop should be complete
  if ~isfield( opt_img, 'i_lp_1' )
  opt_img.i_lp_1 = 1 ;
  end
  if ~isfield( opt_img, 'i_lp_2' )
  opt_img.i_lp_2 = opt_img.n_lp ;
  end

  for i_lp = opt_img.i_lp_1 : opt_img.i_lp_2
    if ( opt_img.n_scnr == 1 )
    disp( sprintf( 'Simulation %i of %i', i_lp, opt_img.n_lp ) )
    end
  opt_img.i_lp = i_lp ;

  % If Keplerian orbits is set, then some individual changes for the scenes happen
  % Update individual planet additions, and local zodi at each iteration.
    if ( opt_img.kplr.do )
    opt_img = get_kepler_imaging_options( opt_img, i_lp, i_bnd ) ;
    % Adding the change in the extragalactic_background position (defined here because it is for all bands)
      if ( opt_img.bckgrnd.add )
        if ( i_bnd == 1 ) && ( i_lp == 1 )
        bckgrnd_ini_ra_mas = opt_img.bckgrnd.arc_ra_mas ;
        bckgrnd_ini_dc_mas = opt_img.bckgrnd.arc_dc_mas ;
        end
      opt_img.bckgrnd.arc_ra_mas = bckgrnd_ini_ra_mas + round( opt_img.bckgrnd.dlt_ra_mas( i_lp ) ) ;
      opt_img.bckgrnd.arc_dc_mas = bckgrnd_ini_dc_mas + round( opt_img.bckgrnd.dlt_dc_mas( i_lp ) ) ;
      end
    else
      if isfield( opt_img.plnt.add, 'phs_ang_dg' )
      % Assign the flux ratio for this epoch and band
      opt_img.plnt.add.flx_rt = squeeze( opt_img.plnt.add.flx_rt_arry( :, i_bnd, : ) ) ;
      % Particular case of only one planet (matlab selects n_lmbdx1, instead of n_plxn_lmbd.
        if ( opt_img.n_pl == 1 )
        opt_img.plnt.add.flx_rt = opt_img.plnt.add.flx_rt' ;
        end
      end
    end

  % Output filenames
  % Removing the extension
    if strcmp( opt_img.scn.nm( end - 3 : end ), '.mat' )
    opt_img.scn.nm = opt_img.scn.nm( 1 : end - 4 ) ;
    end
  out_dr_img = sprintf( '%simage_data/', opt_img.out_dr ) ;
    if ~isdir( out_dr_img )
    mkdir( out_dr_img )
    end
  % Whether this is from a config file, or not (add_label gets filled out in sister):
    if isfield( opt_img, 'case' )
      if isnumeric( opt_img.case )
      fl_out = sprintf( '%ssister_%i%s', out_dr_img, opt_img.add_label, opt_img.tg_lp ) ;
      end
      if ischar( opt_img.case )
      fl_out = sprintf( '%ssister_run_%s%s', out_dr_img, opt_img.add_label, opt_img.tg_lp ) ; 
      end
    else
    fl_out = sprintf( '%ssister_%s_%s_%s_%s', out_dr_img, opt_img.starshade.mode, opt_img.scn.nm, opt_img.tg_nm, opt_img.tg_lp ) ;
    end
  fl_dt_out = [ fl_out '.mat' ] ; opt_img.fl_dt_out = fl_dt_out ;
  fl_fts_out = [ fl_out '.fits' ] ; opt_img.fl_fts_out = fl_fts_out ;

    if ( do_img_f( fl_out, opt_img ) )
    disp( sprintf( 'Considering the band %i/%i', i_bnd, numel( opt_img.lmbd_img_1_nm ) ) ) 
    % Selecting those wavelengths within the limits of the band (lmbd_arry_img_nm is defined in sister.m)
    q_tmp_1 = ( opt_img.lmbd_arry_img_nm >= opt_img.lmbd_img_1_nm( i_bnd ) ) ;
    q_tmp_2 = ( opt_img.lmbd_arry_img_nm < opt_img.lmbd_img_2_nm( i_bnd ) ) ; % Avoiding choosing the last wavelength exactly as the top one, because that would not leave any wavelength interval when getting the conversion to Jy
    q_tmp = q_tmp_1 .* q_tmp_2 ; 
    idx_in = find( q_tmp == 1 ) ;
    opt_img.lmbd_arry_img_nm = opt_img.lmbd_arry_img_nm( idx_in ) ;
    n_lmbd = numel( idx_in ) ;

    % Taking into account if there's local zodiacal light. The first epoch derives a template that is then normalized to the corresponding value of the surafce brightness.
      if ( i_lp == opt_img.i_lp_1 )
      opt_img.lcl_zd.do_rf = 0 ; 
      % For plotting purposes only
      scn_dt_cnv_lcl_zd_e = NaN ; 
      end

    % Looping over wavelength within the band
      for i_lmbd = 1 : n_lmbd
      opt_img.lmbd_tmp_nm = opt_img.lmbd_arry_img_nm( i_lmbd ) ;
      disp( sprintf( 'Considering the wavelength %04d nm (%i/%i)', opt_img.lmbd_tmp_nm, i_lmbd, n_lmbd ) )
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  Convolution of the scene with the Starshade PSF   %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [ scn_dt_cnv_ph_s_tmp scn_dt_cnv_no_plnt_ph_s_tmp opt_img fl_sngl_wl ] = convolve_with_one_wavelength( i_bnd, i_lmbd, opt_img ) ;

      % Adding the local zodiacal light at each wavelength slice
        if ( opt_img.lcl_zd.do )
          if ( i_lp == opt_img.i_lp_1 )
          % Scaling the convolution of the local zodiacal light by the star magnitude (recall that convolve_local_zodi_reference sets a reference surface brightness of 0 V mag/srcsec^2)
          scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) = convolve_local_zodi_reference( i_bnd, i_lmbd, opt_img ) ;
          scn_dt_cnv_ph_s_tmp = scn_dt_cnv_ph_s_tmp + 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) ;
          scn_dt_cnv_no_plnt_ph_s_tmp = scn_dt_cnv_no_plnt_ph_s_tmp + 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) ;
          else 
          % No need to repeat the convolution for the constant signal associated witht he reference local zodiacal light
          scn_dt_cnv_ph_s_tmp = scn_dt_cnv_ph_s_tmp + 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * squeeze( scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) ) ;
          scn_dt_cnv_no_plnt_ph_s_tmp = scn_dt_cnv_no_plnt_ph_s_tmp + 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * squeeze( scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) ) ;
          end
        end

      % Co-adding wavelength slices
        if ( i_lmbd == 1 )
        scn_dt_cnv_ph_s = scn_dt_cnv_ph_s_tmp ;
        scn_dt_cnv_no_plnt_ph_s = scn_dt_cnv_no_plnt_ph_s_tmp ;
        % For plotting purposes
          if ( opt_img.lcl_zd.do )
          scn_dt_cnv_lcl_zd_tmp = 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * squeeze( scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, 1 ) ) ;
          end
        else
        scn_dt_cnv_ph_s = scn_dt_cnv_ph_s + scn_dt_cnv_ph_s_tmp ;
        scn_dt_cnv_no_plnt_ph_s = scn_dt_cnv_no_plnt_ph_s + scn_dt_cnv_no_plnt_ph_s_tmp ;
          if ( opt_img.lcl_zd.do )
          scn_dt_cnv_lcl_zd_tmp = scn_dt_cnv_lcl_zd_tmp + 10^( -opt_img.lcl_zd.mg_v_arcs2 / 2.5 ) * squeeze( scn_lcl_zd_cnv_nrm_ph_s_rf( :, :, i_lmbd ) ) ;
          end
        end
      end % i_lmbd

    % Total counts per integration time
    scn_dt_cnv_e = opt_img.ns.gn * opt_img.ns.exp_tm_s * scn_dt_cnv_ph_s ;
    scn_dt_cnv_no_plnt_e = opt_img.ns.gn * opt_img.ns.exp_tm_s * scn_dt_cnv_no_plnt_ph_s ;
    % For plotting purposes only
      if ( opt_img.lcl_zd.do )
      scn_dt_cnv_lcl_zd_e = opt_img.ns.gn * opt_img.ns.exp_tm_s * scn_dt_cnv_lcl_zd_tmp ;
      end
    scn_ns_e = 0 ;
    % Instrumental noise
      if ( opt_img.ns.do )
      disp( 'Computing the noise of the image ...' )
      scn_ns_e = generate_noise( scn_dt_cnv_e, opt_img ) ;
      % If noise is on, we assume there's a detetcor and counts are integer values
      scn_dt_cnv_e = round( scn_dt_cnv_e ) ;
      scn_dt_cnv_no_plnt_e = round( scn_dt_cnv_no_plnt_e ) ;
      end

    % Finally, return `average' counts
    scn_ns_e = scn_ns_e / opt_img.ns.gn ;
    scn_dt_cnv_e = scn_dt_cnv_e / opt_img.ns.gn ;
    scn_dt_cnv_no_plnt_e = scn_dt_cnv_no_plnt_e / opt_img.ns.gn ;

    % Storing the results
      if ( opt_img.sv_out )
      % Notice that the FITS file output stores the full simulation only
        if ~( opt_img.ftsio )
        % To avoid any conflicts when re-loading files
        opt_img_sv = opt_img ;
        % Saving some memory. Most applications won't need more than 7 digit precision.
          if ( opt_img.sngl_prcsn )
          scn_ns_e = single( scn_ns_e ) ;
          scn_dt_cnv_e = single( scn_dt_cnv_e ) ;
          scn_dt_cnv_no_plnt_e = single( scn_dt_cnv_no_plnt_e ) ;
          scn_dt_cnv_lcl_zd_e = single( scn_dt_cnv_lcl_zd_e ) ;
          end
        save( fl_dt_out, 'scn_ns_e', 'scn_dt_cnv_e', 'scn_dt_cnv_no_plnt_e', 'scn_dt_cnv_lcl_zd_e', 'opt_img_sv' ) ;
        disp( sprintf( 'Data stored: %s', fl_dt_out ) )
        else
        fitswrite( scn_dt_cnv_e, fl_fts_out ) ;
        disp( sprintf( 'Simulation stored in: %s', fl_fts_out ) )
        disp( 'PS: remember that the FITS file stores the full noiseless simulation.' )
        end
      end
    % If do_cnv = 0 then read the file (useful for plotting purposes if the image data are already on disk)
    else
      if ( opt_img.ftsio )
      scn_dt_cnv_e = fitsread( fl_fts_out ) ;
      scn_dt_cnv_no_plnt_e = 0 ; 
      disp( 'PS: remember that the FITS file stores the full noiseless simulation.' )
      else
      load( fl_dt_out )
      end
      % Re-doing the noise (if sister_imaging_band is run once again using stored products, it updates the noise terms)
      if ( opt_img.ns.do )
      disp( 'Computing the noise of the image ...' )
      % Restore the gain term in the signal
      scn_ns_e = generate_noise( opt_img.ns.gn * scn_dt_cnv_e, opt_img ) ;
      % Return average counts again
      scn_ns_e = scn_ns_e / opt_img.ns.gn ;
        % Storing the results
        if ( opt_img.sv_out )
        % To avoid any conflicts when re-loading files
        opt_img_sv = opt_img ;
        % Notice that the FITS file output stores the full simulation only
          if ~( opt_img.ftsio )
            % Saving some memory
            if ( opt_img.sngl_prcsn )
            scn_ns_e = single( scn_ns_e ) ;
            scn_dt_cnv_e = single( scn_dt_cnv_e ) ;
            scn_dt_cnv_no_plnt_e = single( scn_dt_cnv_no_plnt_e ) ;
            end
          save( fl_dt_out, 'scn_ns_e', 'scn_dt_cnv_e', 'scn_dt_cnv_no_plnt_e', 'opt_img_sv' ) ;
          disp( sprintf( 'Data stored with a new noise realization: %s', fl_dt_out ) )
          end
        end
      else
      scn_ns_e = 0 ;
      end
    end % do_img

  % Some plot 
    if ( opt_img.plt.do )
    opt_img = plot_starshade_figures( scn_dt_cnv_e, scn_dt_cnv_no_plnt_e, scn_dt_cnv_lcl_zd_e, scn_ns_e, opt_img, i_lp, i_bnd ) ;
    end
  end % i_lp

% Removing the single wavelength files after the run is finished (they usually take 600 Mb for each run and may end up being a few Gb unless the user is aware of)
% In the case of Keplerian orbits, recover what the user had set up (see sister.m)
  if ( opt_img.kplr.do )
  opt_img.sv_out_2 = opt_img.save_single_wavelength_user ;
  end
  % This file is the one in convolve_with_one_wavelength.m
  if ~( opt_img.sv_out_2 ) % Equivalent to save_single_wavelength
    if exist( 'fl_sngl_wl', 'var' )
    delete( [ fl_sngl_wl( 1 : max( strfind( fl_sngl_wl, '_' ) ) ) '*.mat' ] )
    end
  end

%%%%%%%%%%%%%%%%%
% Sub-functions %
%%%%%%%%%%%%%%%%%

% Checking if there's some work to do
function do_img = do_img_f( fl_out, opt_img )
do_img = 1 ; % By default, do the imaging
fl_dt_out = [ fl_out '.mat' ] ;
  if ( opt_img.ftsio )
  fl_dt_out = [ fl_out '.fits' ] ;
  end
  if exist( fl_dt_out, 'file' ) == 2
    if ( opt_img.redo )
    delete( [ fl_out '.*' ] ) ; % Also delete any extension
    else
    do_img = 0 ;
    disp( sprintf( 'File %s exists. Loading it.', fl_dt_out ) )
    end
  end

% Defining the local options for the case of Keplerian orbits
function opt  = get_kepler_imaging_options( opt, i_lp, i_bnd )
% Position in mas and conversion to pix here
  if isfield( opt.kplr, 'arc_ra_mas' )
  opt.plnt.add.arc_ra_mas = opt.kplr.arc_ra_mas( :, i_lp ) ;
  opt.plnt.add.arc_dc_mas = opt.kplr.arc_dc_mas( :, i_lp ) ;
  opt.plnt.add.arc_ra_px = round( opt.plnt.add.arc_ra_mas / opt.px_scn_mas ) ;
  opt.plnt.add.arc_dc_px = round( opt.plnt.add.arc_dc_mas / opt.px_scn_mas ) ;
  end
% Assign the flux ratio for this epoch and band
opt.plnt.add.flx_rt = squeeze( opt.plnt.add.flx_rt_arry( :, i_lp, i_bnd, : ) ) ;
% Particular case of only one planet (Matlab selects n_lmbdx1, instead of n_plxn_lmbd.
  if ( opt.n_pl == 1 )
  opt.plnt.add.flx_rt = opt.plnt.add.flx_rt' ;
  end

% Assigning the actual value of the local zodiacal light
  if isfield( opt.str, 'lcl_zd' )
  opt.lcl_zd.mg_v_arcs2 = opt.str.lcl_zd.mg_v_arcs2( i_lp ) ;
    if isnan( opt.lcl_zd.mg_v_arcs2 )
    disp( 'WARNING (NaN): the surface brightness density of the local zodiacal light cannot be derived because the line of sight is too close to the Sun.' )
    end
  end

% Tag name associated with the loop
opt.tg_lp = sprintf( '_%04i', i_lp ) ;
% In case it's been set externally (loop in sister)
  if isfield( opt, 'tg_lp_ext' )
  opt.tg_lp = opt.tg_lp_ext ;
  end
% redo_2
  if isfield( opt, 'redo_2_arry' )
  opt.redo_2 = opt.redo_2_arry( i_lp ) ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolving the local zodi and normalizing it to speed up the convolution process in a Keplerian case %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not output the options from this function or they may overwrite the actual options of the simulation
function scn_lcl_zd_cnv_nrm_ph_s = convolve_local_zodi_reference( i_bnd, i_lmbd, opt )
% Modifying some parameters to only convolve a constant signal
opt.scn.do = 1 ;
opt.slr_glnt.do = 0 ;
opt.lcs.do = 0 ;
opt.str.mg = 40 ; % Significantly less bright than the local zodiacal light, which is brighter than 24.
opt.zd.fct = 0 ;
opt.planets.add.do = 0 ;
opt.kpr.do = 0 ;
opt.bckgrnd.add = 0 ;
opt.lcl_zd.do = 1 ;
opt.lcl_zd.do_rf = 1 ; % It forces the convolution with an equivalent surface brightness equal to 0th magnitude that is the reference convolution.
opt.sv_out_2 = 0 ; % Do not store any of these results on disk
opt.plnt.add.do = 0 ; % Do not add any planets
opt.plnt.rmv.do = 0 ; % Do not remove any planets
opt.cb.do = 0 ; % Do not produce any data cubes
scn_lcl_zd_cnv_nrm_ph_s =  convolve_with_one_wavelength( i_bnd, i_lmbd, opt ) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
