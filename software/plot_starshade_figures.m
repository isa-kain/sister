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
function opt = plot_starshade_figures( scn_dt_cnv_e, scn_dt_cnv_no_plnt_e, scn_dt_cnv_lcl_zd_e, scn_ns_e, opt, i_lp, i_bnd )
% Plotting the starshade simulation
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% standard deviation of the noise in the image (if it is noiseless, it would be 0)
opt.plt.std_ns = std( scn_ns_e( : ) ) ;

% Get the default options for the plots
opt = get_plot_options( opt, i_lp, i_bnd ) ;

% Combinations of signal and noise
fct_ns     = [ 1, 1, 0, 0, 1 ] ;
fct_no_bck = [ 0, -1, 0, -1, 0 ] ;
% The local zodiacal light can be removed to ease some plots that have a changing level, like in the Keplerian orbits.
fct_no_lzl = [ 0, 0, 0, 0, -1 ] ;
spttl_tmp = { '', '(Background subtracted)', '(Noiseless)', '(Only planets)', '' } ;
n_add_lbl = numel( opt.plt.add_plnt_lbl ) ;

% Setting the properties that go to plt_mnmx
opt_plt.px_unt = opt.plt.px_unt ;
opt_plt.pbaspect = opt.plt.pbaspect ;
opt_plt.clrbr = opt.plt.clrbr ;
opt_plt.clrbr_lbl_sz = opt.plt.clrbr_lbl_sz ;
opt_plt.log_10 = opt.plt.log_10 ;
% If the AU axis wants to be plotted +/- instead of -/+
opt_plt.noflip = opt.plt.noflp ;

% Plotting
  if ( opt.n_scnr == 1 ) && ( i_lp == 1 )
  figure( 1 )
  clf
  setwinsize( gcf, 600, 500 )
  end

% Figure counter
i_fg = 1 ;
% For each combination
  for i_cmb = opt.plt.cmb_lst
    % Reminder about FITS file use
    if ( opt.ftsio )
      if ( fct_no_bck( i_cmb ) )
      disp( 'PS: remember that when using FITS files for input/output, the simulation without planets is not stored. Thus, the image combination that subtracts the background will still have the background on it.' )
      end
    end
  % For each zooming
    for i_zm = 1 : numel( opt.plt.zm_lst )
    dt_plt = scn_dt_cnv_e + fct_no_bck( i_cmb ) * scn_dt_cnv_no_plnt_e + fct_ns( i_cmb ) * scn_ns_e ;
    % Removing the local zodiacal light
      if fct_no_lzl( i_cmb )
        if isnan( scn_dt_cnv_lcl_zd_e )
        disp( 'WARNING (NaN): the value of the local zodiacal light to be removed is NaN. Review whether opt.local_zodi.do=1. If so, make sure that opt.local_zodi.mag_v_arcsec2 is not a NaN, due to a line of sight too close to the Sun. If opt.local_zodi.do=0, and you get this warning, it may be because your plot options include the subtraction of an (inexistent) local zodiacal light. Check opt.plot.combination_list' )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop 
        end
      dt_plt = dt_plt + fct_no_lzl( i_cmb ) * scn_dt_cnv_lcl_zd_e ;
      end
    % Same scale for all plots (If it is changed to be estimated at every zoom, which would record the background for the last zoom instead of the first one, then clear up the fields min and max: opt = rmfield( opt.plot, 'min' ) ; opt = rmfield( opt.plot, 'max' ) ;
      if ( i_lp == 1 ) && ( i_zm == 1 ) 
       % Some percentile. Not using prctile from Matlab since it belongs to a Toolbox, and may cause some trouble. This is fast (<0.003 sec)
        if ~fct_no_lzl( i_cmb )
        dt_srt = sort( scn_dt_cnv_e( : ) ) ;
          if ~sum( isnan( scn_dt_cnv_lcl_zd_e( : ) ) )
          opt.plt.bckgrnd_scn = median( scn_dt_cnv_lcl_zd_e( : ) ) ;
          else
          opt.plt.bckgrnd_scn = median( dt_srt ) ;
          end
        else
        dt_srt = sort( scn_dt_cnv_e( : ) + fct_no_lzl( i_cmb ) * scn_dt_cnv_lcl_zd_e( : ) ) ;
        opt.plt.bckgrnd_scn = median( dt_srt ) ;
        end
        if isfield( opt.plot, 'min' )
        mn_plt = opt.plt.bckgrnd_scn + opt.plt.mn ; % Adding the background in case there's the loca zodiacal light, for instance
        % Keeping track of whether the minimum value was set in the configuration file
        mn_st = 1 ;
        else
        mn_plt = opt.plt.bckgrnd_scn - 3 * opt.plt.std_ns ;
         % For images without noise (contrast that shows background and bright features. Not easy for all cases)
          if ~( opt.ns.do )
          mn_plt = opt.plt.bckgrnd_scn ;
          end
        mn_st = 0 ;
        end
      opt.plt.mn = mn_plt ;
      opt.plot.min = mn_plt ;
        if isfield( opt.plot, 'max' )
        mx_plt = opt.plt.bckgrnd_scn + opt.plt.mx ; % Adding the background in case there's the loca zodiacal light, for instance
        % Keeping track of whether the maximum value was set in the configuration file
        mx_st = 1 ;
        else
        mx_plt = opt.plt.bckgrnd_scn + opt.plt.std_ns * opt.plt.n_sgm_ns ;
        % For images without noise
          if ~( opt.ns.do )
          mx_plt = max( dt_srt( 1 : round( numel( dt_srt ) / 100 * 99.9 ) ) ) ;
          end
        mx_st = 1 ;
        end
      opt.plt.mx = mx_plt ;
      opt.plot.max = mx_plt ;
      end
    % Additional factor to better show the subtracted image
    fct_no_bck_mx = 1 ;
    opt_plt.clrbr_lbl = opt.plt.clrbr_unt ;
      if ( fct_no_bck( i_cmb ) == -1 ) && numel( scn_dt_cnv_no_plnt_e ) > 1 && sum( sum( scn_dt_cnv_no_plnt_e ) )
    % If there's background, and both, with and without subtraction, are to be plotted
        if ( numel( find( opt.plt.cmb_lst == 1 ) ) || numel( find( opt.plt.cmb_lst == 3 ) ) ) && fct_no_bck( i_cmb ) , fct_no_bck_mx = 0.5 ; end
        if ( opt.ns.do ) && ~mn_st && ~mx_st
        opt_plt.clrbr_lbl = sprintf( 'e (background-3+%1.0f %s_{N} )', opt.plt.n_sgm_ns * fct_no_bck_mx, '\sigma' ) ;
        else
        opt_plt.clrbr_lbl = 'e' ;
        end
      end
    % If the local zodiacal light has been removed:
      if ( fct_no_lzl( i_cmb ) == -1 )
        if ( opt.ns.do )
        opt_plt.clrbr_lbl = sprintf( 'e (-3+%1.0f %s_{N} )', opt.plt.n_sgm_ns * fct_no_bck_mx, '\sigma' ) ;
        else
        opt_plt.clrbr_lbl = 'e' ;
        end
      end
    opt_plt.mnmx = [ opt.plt.mn, opt.plt.mn + fct_no_bck_mx * ( opt.plt.mx - opt.plt.mn ) ] ;
    opt_plt.zm = opt.plt.zm_lst( i_zm ) ;
     if ~isfield( opt, 'run' )
     clf
     else
    end
    plt_mnmx( dt_plt, opt_plt  ) ;
    % For boundaries (this may give problems when using Nx servers to work remotely. If so, set xy_fg_lm to a big number, xy_fg_lm=1e5;. The issue is that it may plot some ticks of the planets outside the central image when using Kepler).
    h_fg = get( gcf ) ;
    xy_fg_lm = abs( h_fg.CurrentAxes.XLim( 2 ) ) ; % It's always a square
      % Adding the planet labels
      for i_lbl = 1 : n_add_lbl
        if opt.plt.add_plnt_lbl( i_lbl )
        fct_lbl_ps = opt.fwhm_px / i_zm ;
        % Just a warning message if the number of labels is greater than the number of planets (if there were less, the default options would stop the run)
          if numel( opt.plt.plnt_lbl ) ~= opt.n_pl
          disp( sprintf( 'WARNING: Using the first %i labels, since there are %i planets and %i labels', opt.n_pl, opt.n_pl, numel( opt.plt.plnt_lbl ) ) ) 
          end
          for i_pl = 1 : opt.n_pl
          tan_alph_pl = abs( opt.plnt.add.arc_dc_mas( i_pl ) / opt.plnt.add.arc_ra_mas( i_pl ) ) ;
            if isinf( tan_alph_pl ), tan_alph_pl = 1e6 ; end
          ps_plnt_ra = opt.plnt.add.arc_ra_px( i_pl ) + fct_lbl_ps * opt.plt.plnt_lbl_ps_fwhm( i_pl ) / sqrt( 1 + tan_alph_pl^2 ) * sign( opt.plnt.add.arc_ra_mas( i_pl ) ) ;
          ps_plnt_dc = opt.plnt.add.arc_dc_px( i_pl ) + fct_lbl_ps * opt.plt.plnt_lbl_ps_fwhm( i_pl ) * tan_alph_pl / sqrt( 1 + tan_alph_pl^2 ) * sign( opt.plnt.add.arc_dc_mas( i_pl ) ) ;
          ps_plnt_ra_arry( i_pl ) = ps_plnt_ra * opt.px_scn_mas * opt.plt.fct_plt_mas ;
          ps_plnt_dc_arry( i_pl ) = ps_plnt_dc * opt.px_scn_mas * opt.plt.fct_plt_mas ;
          sz_txt = opt.plt.plnt_lbl_sz( i_pl ) * ( 1 + ( opt.plt.zm_lst( i_zm ) - 1 ) / 10 ) ;
          % If it falls within the plotted area
            if ( abs( ps_plnt_ra_arry( i_pl ) ) <= xy_fg_lm ) && ( abs( ps_plnt_dc_arry( i_pl ) ) <= xy_fg_lm )
            text( ps_plnt_ra_arry( i_pl ), ps_plnt_dc_arry( i_pl ), opt.plt.plnt_lbl{ i_pl }, 'Color', opt.plt.plnt_lbl_clr{ i_pl }, 'FontSize', sz_txt )
            end
          end
        end % opt.plt.add_plnt_lbl

        % Adding the shape of the starshade as a circle with full precision
        if ( opt.plt.ss_crcl.do ) || ( opt.plt.ss_ptls.do )
        fl_nm_xyzVals =  sprintf( '%s/locus/in/%s_xyzVals.mat', opt.scene_dir, opt.starshade.nominal_filename ) ;
        % it brings xVals and yVals
          if exist( fl_nm_xyzVals, 'file' ) == 2
          load( fl_nm_xyzVals )
          else
          % Petal locus. This is very fast ~0.01 sec all
          load( sprintf( '%s/locus/in/%s.mat', opt.scene_dir, opt.starshade.nominal_filename ) )
          vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, opt.n_ptl, {}) ;
          tmpxVals = []; tmpyVals = []; tmpzVals = [];
            for j = 1 : opt.n_ptl
            tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)];
            tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)];
            tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)];
            end
          xVals = [tmpxVals tmpxVals(1)]; yVals = [tmpyVals tmpyVals(1)]; zVals = [tmpzVals tmpzVals(1)];
          save( fl_nm_xyzVals, 'xVals', 'yVals', 'zVals' ) ;
          end
        end

        if ( opt.plt.ss_crcl.do )
        hold all
        geo_iwa_mas = sqrt( ( max( xVals( : ) ) - min( xVals( : ) ) )^2 + ( max( yVals( : ) ) - min( yVals( : ) ) )^2 )  / 2 / sqrt( 2 ) / opt.dst_strshd_tlscp_m * 180 / pi * 3600e3 ;
        plot( geo_iwa_mas * opt.plt.fct_plt_mas * cos( 0 : 0.001 : 2 * pi ), geo_iwa_mas * opt.plt.fct_plt_mas * sin( 0 : 0.001 : 2 * pi ), 'Color', opt.plt.ss_crcl.clr, 'LineWidth', opt.plt.ss_crcl.ln_wdth, 'LineStyle', opt.plt.ss_crcl.ln_st )
        hold off
        end

        % Adding the shape of the starshade with its petals
        if ( opt.plt.ss_ptls.do )
        hold all
        % it brings xVals and yVals. Nominal starshade or perturbed one if it is set in the simulation
          if ~( opt.lcs.do )
          load( sprintf( '%s/locus/in/%s.mat', opt.scene_dir, opt.starshade.nominal_filename ) )
          else
          load( sprintf( '%s/locus/in/%s.mat', opt.scene_dir, opt.locus.perturbed_filename ) )
          end
        % Factor to go from meters mas and then to the units in the plot
        fct_m2plt = 180 / pi * 3600e3 / opt.dst_strshd_tlscp_m * opt.plt.fct_plt_mas ;
        % If the scene was displaced laterally, apply the inverse shift to the starshade petals
          if ~isnan( opt.scn.ra_shft_px )
          xVals = xVals - opt.scn.ra_shft_m ;
          yVals = yVals - opt.scn.dc_shft_m ; 
          end
        plot( xVals( 1 : round( 1 / opt.plt.ss_ptls.n_vl ) : end ) * fct_m2plt, yVals( 1 : round( 1 / opt.plt.ss_ptls.n_vl ) : end ) * fct_m2plt, '.', 'Color', opt.plt.ss_ptls.clr )
        hold off
        end

        if isfield( opt.plt, 'x_lbl' )
        xlabel( opt.plt.x_lbl, 'Fontsize', opt.plt.fnt_sz) ;
        ylabel( opt.plt.y_lbl, 'Fontsize', opt.plt.fnt_sz ) ;
        end
    
      % Size of tick marks
      ax = gca ;
      ax.FontSize = opt.plt.fnt_sz ;

        if isfield( opt.plt, 'spttl' )
        spttl = sprintf( '%s %s', opt.plt.spttl, spttl_tmp{ i_cmb } ) ; ;
        else
        spttl = spttl_tmp ;
        end
      
      % Title (double replacement, in case there was a subindex already)
      title( strrep( strrep( opt.plt.ttl, '_', '\_' ), '\_{', '_{' ), 'FontSize', opt.plt.fnt_sz )
      
     
      % Super title
        if ( opt.n_scnr == 1 )
        suptitle( strrep( opt.plt.spttl, '_', '\_' ) ) ;
        end
      % Store image
        if ( opt.sv_plt )
          if opt.plt.add_plnt_lbl( i_lbl )
          fg_plnt_lbl = '';
          else
          fg_plnt_lbl = '_no_label' ;
          end
          % For a single case
          if ( opt.n_scnr == 1 )
            if isfield( opt, 'cs_nm' )
              if numel( opt.add_label )
              fg_nm = sprintf( 'sister_%s_cmb_%i_zm_%i%s%s', opt.add_label, i_cmb, i_zm, fg_plnt_lbl, opt.tg_lp ) ;
              else
              fg_nm = sprintf( 'sister_%s_cmb_%i_zm_%i%s%s', opt.cs_nm, i_cmb, i_zm, fg_plnt_lbl, opt.tg_lp ) ;
              end
            else
            fg_nm = sprintf( '%s_%s_cmb%i_zm%i%s_%s%s', opt.scn.nm, opt.starshade.mode, i_cmb, i_zm, fg_plnt_lbl, opt.tg_nm, opt.tg_lp ) ;
            end
          img = getframe( gcf ) ;
          % Creating the directory if it does not exist
          if ~isdir( opt.fg_dr )
          mkdir( opt.fg_dr ) ;
          end
          % Static image or video
            if ( opt.vd.do )
            % Creating the directory if it does not exist
                if ~isdir( opt.vd_dr )
                mkdir( opt.vd_dr ) ;
                end
            % Write to the GIF File
            fg_nm_vd = sprintf( '%s%s.gif', opt.vd_dr, fg_nm( 1 : end - 5 ) ) ;
            img_2 = frame2im( img );
            [ imind, cm ] = rgb2ind(img_2, 256 ) ;
              if i_lp == 1
              imwrite( imind, cm, fg_nm_vd, 'gif', 'Loopcount', inf ) ;
              else
              imwrite( imind, cm, fg_nm_vd, 'gif', 'WriteMode', 'append', 'DelayTime', opt.vd_dly_s ) ;
              end
            else
            imwrite( img.cdata, [ opt.fg_dr fg_nm '.' opt.frmt_fg ] ) ; 
            end % opt.vd.do
          end % opt.sv_fg
        end % opt.n_scnr == 1
      end % i_lbl
    end % i_zm
  end % i_cmb

% Additional plot with the orbits (once), if there is an orbit
  if ( opt.kplr.do ) && ( i_lp == 1 ) && ( opt.n_pl ) && ( opt.n_scnr == 1 )
  opt_plt.clrbr = 0 ;
  opt_plt.pbaspect = opt.plt.pbaspect ;
  dt_plt = 0 * scn_dt_cnv_e ;
  opt_plt.mnmx = [ 0, 160 ] ;
  n_stmp = numel( opt.kplr.tm_yr_obs ) ;
  figure( i_fg + 1 )
  setwinsize( gcf, 600, 500 )
  % For the different zooming options
    for i_zm = 1 : numel( opt.plt.zm_lst )
    clf
    hold all
    opt_plt.zm = opt.plt.zm_lst( i_zm ) ;
    plt_mnmx( dt_plt, opt_plt ) ;
    % For boundaries (this may give problems when using Nx servers to work remotely. If so, set xy_fg_lm to a big number, xy_fg_lm=1e5;. The issue is that it may plot some ticks of the planets outside the central image when using Kepler).
    h_fg = get( gcf ) ;
    xy_fg_lm = abs( h_fg.CurrentAxes.XLim( 2 ) ) ; % It's always a square
    % Central star
    text( 0, 0, 'o', 'Color', opt.plt.plnt_lbl_clr{ 1 }, 'FontSize', 4 )
    fct_lbl_ps = 1 - (15/75) * ( opt.plt.zm_lst( i_zm  ) - 1 ) / 2.3 ; % To place the labels for the planets
      for i_stmp = 1 : n_stmp
        for i_pl = 1 : opt.n_pl
        % Planet position
        ps_plnt_ra = opt.kplr.arc_ra_mas( i_pl, i_stmp ) * opt.plt.fct_plt_mas ;
        ps_plnt_dc = opt.kplr.arc_dc_mas( i_pl, i_stmp ) * opt.plt.fct_plt_mas ;
        sz_txt = opt.plt.plnt_lbl_sz( i_pl ) * ( 1 + ( opt.plt.zm_lst( i_zm ) - 1 ) / 6 ) / 4 ;
        % If it falls within the plotted area
          if ( abs( ps_plnt_ra ) <= xy_fg_lm ) && ( abs( ps_plnt_dc ) <= xy_fg_lm )
          text( ps_plnt_ra, ps_plnt_dc, '+', 'Color', opt.plt.plnt_lbl_clr{ i_pl }, 'FontSize', sz_txt )
          end
        % Planet label
          if ( i_stmp == 1 )
          tan_alph_pl = opt.plnt.add.arc_dc_mas( i_pl ) / opt.plnt.add.arc_ra_mas( i_pl ) ;
          ps_plnt_ra = opt.plnt.add.arc_ra_px( i_pl ) * opt.px_scn_mas * opt.plt.fct_plt_mas ;
          ps_plnt_dc = opt.plnt.add.arc_dc_px( i_pl ) * opt.px_scn_mas * opt.plt.fct_plt_mas ;
          sz_txt = opt.plt.plnt_lbl_sz( i_pl ) * ( 1 + ( opt.plt.zm_lst( i_zm ) - 1 ) / 6 ) ;
          % If it falls within the plotted area
            if ( abs( ps_plnt_ra ) <= xy_fg_lm ) && ( abs( ps_plnt_dc ) <= xy_fg_lm )
            text( ps_plnt_ra, ps_plnt_dc, opt.plt.plnt_lbl{ i_pl }, 'Color', opt.plt.plnt_lbl_clr{ i_pl }, 'FontSize', sz_txt )
            end
          end % i_stmp
        end % i_pl
      end % i_stmp
    if isfield( opt.plt, 'spttl' )
    spttl = sprintf( '%s: orbits', opt.plt.spttl ) ; ;
    else
    spttl = 'Orbits' ;
    end
    suptitle( strrep( spttl, '_', '\_' ) )
    % Title
    ttl = sprintf( 'YEARS: %4.2f to %4.2f', opt.kplr.tm_yr_obs( 1 ), opt.kplr.tm_yr_obs( end ) ) ;
    title( ttl )
    xlabel( 'AU' )
    ylabel( 'AU' )
    % Store image (not to be deleted: adding orbit_ at the beginning)
      if ( opt.sv_plt )
        if isfield( opt, 'case' )
          if isnumeric( opt.case )
          fg_nm = sprintf( 'orbit_config_%i_cmb_%i_zm_%i%s.%s', opt.case, i_cmb, i_zm, fg_plnt_lbl, opt.frmt_fg ) ;
          end
          if ischar( opt.case )
          fg_nm = sprintf( 'orbit_config_%s_cmb_%i_zm_%i%s.%s', opt.case, i_cmb, i_zm, fg_plnt_lbl, opt.frmt_fg ) ;
          end
        else
          if isfield( opt, 'add_label' )
          fg_nm = [ 'orbit_' opt.scn.nm '_' opt.tg_nm_1 '_' opt.tg_nm_2 '_zm' num2str( i_zm ) '_' opt.add_label '.' opt.frmt_fg ] ;
          else
          fg_nm = [ 'orbit_' opt.scn.nm '_' opt.tg_nm_1 '_' opt.tg_nm_2 '_zm' num2str( i_zm ) '.' opt.frmt_fg ] ;
          end
        end
      img = getframe( gcf ) ;
      % Creating the directory if it does not exist
        if ~isdir( opt.fg_dr )
        mkdir( opt.fg_dr ) ;
        end
      imwrite( img.cdata, [ opt.fg_dr fg_nm ] ) ;
      end
    end % i_zm
  %figure( 1 ) % Return to the figure with the image of the planets
  end % i_lp == 1 (special figure of the orbits)

%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS %
%%%%%%%%%%%%%%%%%

% Default options for plotting
function opt = get_plot_options( opt, i_lp, i_bnd )

% Lower case fields
opt = lower_opt( opt ) ;

% Overplot of the starshade shape
% A) as a circle
  if ~isfield( opt.plot, 'starshade_circle' )
  opt.plt.ss_crcl.do = 0 ; % Not plotted by default
  opt.plot.starshade_circle.do = 0 ;
  else
    if ~isfield( opt.plot.starshade_circle, 'do' )
    opt.plt.ss_crcl.do = 0 ; % Not plotted by default
    opt.plot.starshade_circle.do = 0 ;
    else
    opt.plt.ss_crcl.do = opt.plot.starshade_circle.do ;
    end
    % Color
    if ~isfield( opt.plot.starshade_circle, 'color' )
    opt.plt.ss_crcl.clr = 'w' ; % white
    opt.plot.starshade_circle.color = opt.plt.ss_crcl.clr ; % Assigned to be used in get_running_options.m/sister.m
    else
    opt.plt.ss_crcl.clr = opt.plot.starshade_circle.color ;
    end
    % Line width
    if ~isfield( opt.plot.starshade_circle, 'line_width' )
    opt.plt.ss_crcl.ln_wdth = 1.5 ; % default (fine for dots. For lines 1 may be better)
    opt.plot.starshade_circle.line_width = opt.plt.ss_crcl.ln_wdth ;
    else
    opt.plt.ss_crcl.ln_wdth = opt.plot.starshade_circle.line_width ; % Assigned to be used in get_running_options.m/sister.m
    end
    % Line style
    if ~isfield( opt.plot.starshade_circle, 'line_style' )
    opt.plt.ss_crcl.ln_st = ':' ; % dots
    opt.plot.starshade_circle.line_style = opt.plt.ss_crcl.ln_st ;
    else
    opt.plt.ss_crcl.ln_st = opt.plot.starshade_circle.line_style ; % Assigned to be used in get_running_options.m/sister.m
    end
  end % field( opt.plot, 'starshade_circle' )  

% B) all the petals
  if ~isfield( opt.plot, 'starshade_petals' )
  opt.plt.ss_ptls.do = 0 ; % Not plotted by default
  opt.plot.starshade_petals.do = 0 ;
  else
    if ~isfield( opt.plot.starshade_petals, 'do' )
    opt.plt.ss_ptls.do = 0 ; % Not plotted by default
    opt.plot.starshade_petals.do = 0 ;
    else
    opt.plt.ss_ptls.do = opt.plot.starshade_petals.do ;
    end
    % Density of pixels (too many may interfere with the simulated image itself)
    if ~isfield( opt.plot.starshade_petals, 'n_values' )
    opt.plt.ss_ptls.n_vl = 0.1 ; % Fraction of all the values in xVals, yVals that will be plotted. 
    opt.plot.starshade_petals.n_values = opt.plt.ss_ptls.n_vl ;
    else
    opt.plt.ss_ptls.n_vl = opt.plot.starshade_petals.n_values ;
    end
    % Check of consistency
    if ( opt.plt.ss_ptls.n_vl > 1 )
    disp( sprintf( 'WARNING. Choose a value for opt.plot.starshade_petals.n_values that is <= 1. Your choce was %f. Set it to 1. Continuing.', opt.plt.ss_ptls.n_vl ) ) ;
    end
    % Color
    if ~isfield( opt.plot.starshade_petals, 'color' )
    opt.plt.ss_ptls.clr = 'w' ; % white
    opt.plot.starshade_petals.color = opt.plt.ss_ptls.clr ; % Assigned to be used in get_running_options.m/sister.m
    else
    opt.plt.ss_ptls.clr = opt.plot.starshade_petals.color ;
    end
  end % field( opt.plot, 'starshade_petals' )

% If there are no planets, then set n_pl to zero
  if ~isfield( opt, 'n_pl' )
  opt.n_pl = 0 ;
  end

% Labels for the planets
lbl_plnt_str = { 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'x', 'y', 'z' } ;
  if ~isfield( opt.plot, 'planet_label' )
    for i_pl = 1 : opt.n_pl
    opt.plt.plnt_lbl{ i_pl } = lbl_plnt_str{ i_pl } ;
    end
  else
    if numel( opt.plot.planet_label ) < opt.n_pl
    disp( sprintf( 'There are not enough labels defined (%i) for the number of planets (%i). Please, fix either one.', numel( opt.plt.plnt_lbl ), opt.n_pl ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  opt.plt.plnt_lbl = opt.plot.planet_label ;
  end

% Color of the planet labels
  if ~isfield( opt.plot, 'planet_label_color' )
    for i_pl = 1 : opt.n_pl
    opt.plt.plnt_lbl_clr{ i_pl } = 'w' ;
    end
  else
  opt.plt.plnt_lbl_clr = opt.plot.planet_label_color ;
  end

% Size of the font
  if ~isfield( opt.plot, 'planet_label_fontsize' )
    for i_pl = 1 : opt.n_pl
    opt.plt.plnt_lbl_sz( i_pl ) = 10 ;
    end
  else
  opt.plt.plnt_lbl_sz = opt.plot.planet_label_fontsize ;
  end

% Planet label position: concentric + some pixels
  if ~isfield( opt.plot, 'planet_label_pos_fwhm' )
    for i_pl = 1 : opt.n_pl
    opt.plt.plnt_lbl_ps_fwhm( i_pl ) = 1 ; % times FWHM
    end
  else
    if numel( opt.plot.planet_label_pos_fwhm ) == 1
      for i_pl = 1 : opt.n_pl
      opt.plt.plnt_lbl_ps_fwhm( i_pl ) = opt.plot.planet_label_pos_fwhm ;
      end
    else
    opt.plt.plnt_lbl_ps_fwhm = opt.plot.planet_label_pos_fwhm ;
      if numel( opt.plot.planet_label_pos_fwhm ) == 1
        for i_pl = 1 : opt.n_pl
        opt.plt.plnt_lbl_ps_fwhm( i_pl ) = opt.plot.planet_label_pos_fwhm ;
        end
      elseif numel( opt.plot.planet_label_pos_fwhm ) < opt.n_pl
      disp( sprintf( 'There are not enough elements in opt.plot.planet_label_pos_fwhm (there are %i) for all the planets in the simulation (%i). Please, fix either.', numel( opt.plot.planet_label_pos_fwhm ), opt.n_pl ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  end
% Changing it into pixels
% Longest wavelength
opt.fwhm_px = 1.02 * opt.lmbd_img_2_nm( i_bnd ) * 1e-9 / opt.dmtr_tlscp_m * 180 / pi * 60 * 60 * 1000 / opt.px_scn_mas ;

% Number of zooms
  if ~isfield( opt.plot, 'zoom_list' )
  opt.plt.zm_lst = [ 1 ] ; % It may be multiple values, e.g., [ 1 20 ] and it will plot 2 figures with each zoom in factor.
  else
  opt.plt.zm_lst = opt.plot.zoom_list ;
  end

% Check of consistency
  for i_zm = 1 : numel( opt.plt.zm_lst )
    if opt.plt.zm_lst( i_zm ) < 1
    disp( 'All the values in opt.plot.zoom_list must be greater or equal than 1. Please, fix it.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end

% Min/Max values
  if ~isfield( opt.plot, 'min' )
  opt.plt.mn = 0 ;
  else
  opt.plt.mn = opt.plot.min ;
  end
  if ~isfield( opt.plot, 'max' )
  opt.plt.mx = 0 ;
  else
  opt.plt.mx = opt.plot.max ;
  end

% Some number of times the std of the noise
  if ~isfield( opt.plot, 'n_sigma_noise' )
  opt.plt.n_sgm_ns = 10 ;
  else
  opt.plt.n_sgm_ns = opt.plot.n_sigma_noise ;
  end

% Axis labels
  if ~isfield( opt.plot, 'x_label' )
  opt.plt.x_lbl = 'AU' ;
  else
  opt.plt.x_lbl = opt.plot.x_label ;
  end
  if ~isfield( opt.plot, 'y_label' )
  opt.plt.y_lbl = 'AU' ;
  else
  opt.plt.y_lbl = opt.plot.y_label ;
  end

% Fontsize for title, axis labels, and tick marks
  if ~isfield( opt.plot, 'fontsize' )
  opt.plt.fnt_sz = 14 ;
  else
  opt.plt.fnt_sz = opt.plot.fontsize ;
  end

% Check consistency for these plots (smae units on both X/Y axis)
  if ~strcmp( opt.plt.x_lbl, opt.plt.y_lbl )
  disp( sprintf( 'Figure: The x label (%s) and y label (%s) are different. Fix opt.plot.x_label, opt.plot.y_label. Stopped.', opt.plt.x_lbl, opt.plt.y_lbl ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% If labels are AU, use the proper units
  if strcmp( lower( opt.plt.x_lbl ), 'au' )
  opt.plt.px_unt = opt.px_cmr_mas * opt.mas2au ;
  end
  if strcmp( lower( opt.plt.x_lbl ), 'mas' )
  opt.plt.px_unt = opt.px_cmr_mas ;
  end
  % Check either case has been used
  if ~isfield( opt.plt, 'px_unt' )
  disp( sprintf( 'Please, define either opt.plot.x_label and opt.plot.y_label to be both either AU or mas. Right now, opt.plot.x_label and opt.plot.y_label are %s', opt.plt.x_lbl ) )
  end

% Whether mas or AU are used for the planet positions
  if strcmp( lower( opt.plt.x_lbl ), 'au' )
  opt.plt.fct_plt_mas = opt.mas2au ;
  else
  opt.plt.fct_plt_mas = 1 ;
  end

% Whether a colorbar will be plotted or not
  if ~isfield( opt.plot, 'colorbar' )
  opt.plt.clrbr = 1 ; % Default, 1, it will be plotted. Set it to 0 in order to not plot it.
  else
  opt.plt.clrbr = opt.plot.colorbar ;
  end

% Units of the colorbar
  if ~isfield( opt.plot, 'colorbar_unit' )
  opt.plt.clrbr_unt = sprintf( 'e (background-3+%1.0f %s_{N} )', opt.plt.n_sgm_ns, '\sigma' ) ;
    if ~( opt.ns.do )
    opt.plt.clrbr_unt = 'e' ;
    end
  else
  opt.plt.clrbr_unt = opt.plot.colorbar_unit ;
  end

% Font size of the colorbar label
  if ~isfield( opt.plot, 'colorbar_fontsize' )
  opt.plt.clrbr_lbl_sz = 12 ;
  else
  opt.plt.clrbr_lbl_sz = opt.plot.colorbar_fontsize ;
  end

% PBaspect
  if ~isfield( opt.plot, 'pbaspect' )
  opt.plt.pbaspect = 1 ;
  else
  opt.plt.pbaspect = opt.plot.pbaspect ;
  end

% Some time stamp
  if isfield( opt.kplr, 'tm_yr_obs' )
  opt.plt.tm_stmp = sprintf( ' YEAR=%4.2f', opt.kplr.tm_yr_obs( i_lp ) ) ;
  else
  opt.plt.tm_stmp = '' ;
  end

% Title of the scene
  if ~isfield( opt.plot, 'title' )
  % The number of frames only 
     if ( opt.ns.do )
       if ( opt.ns.exp_tm_s / 86400 >= 1 )
       n_dy = opt.ns.exp_tm_s / 86400 ;
         if ( n_dy == 1 ), lbl_dy = 'day' ; else, lbl_dy = 'days' ; end
       opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f %s of observation (%i frames, %s_{N}=%3.0f e) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), n_dy, lbl_dy, opt.ns.n_frm, '\sigma', opt.plt.std_ns, opt.plt.tm_stmp ) ;
       elseif ( opt.ns.exp_tm_s / 3600 >= 1 )
       n_hr = opt.ns.exp_tm_s / 3600 ;
         if ( n_hr == 1 ), lbl_hr = 'hour' ; else, lbl_hr = 'hours' ; end
       opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f %s of observation (%i frames, %s_{N}=%3.0f e) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), n_hr, lbl_hr, opt.ns.n_frm, '\sigma', opt.plt.std_ns, opt.plt.tm_stmp ) ;
       else
       n_mn = opt.ns.exp_tm_s / 60 ;
         if ( n_mn == 1 ), lbl_mn = 'min' ; else, lbl_mn = 'mins' ; end
        opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f %s of observation (%i frames, %s_{N}=%3.0f e) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), n_mn, lbl_mn, opt.ns.n_frm, '\sigma', opt.plt.std_ns, opt.plt.tm_stmp ) ;
       end
    else
      if ( opt.ns.exp_tm_s / 86400 >= 1 )
      n_dy = opt.ns.exp_tm_s / 86400 ;
        if ( n_dy == 1 ), lbl_dy = 'day' ; else, lbl_dy = 'days' ; end
      opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f %s of observation (noiseless) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), n_dy, lbl_dy, opt.plt.tm_stmp ) ;
      elseif ( opt.ns.exp_tm_s / 3600 >= 1 )
      n_hr = opt.ns.exp_tm_s / 3600 ;
         if ( n_hr == 1 ), lbl_hr = 'hour' ; else, lbl_hr = 'hours' ; end
      opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f hours of observation (noiseless) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), opt.ns.exp_tm_s / 3600, opt.plt.tm_stmp ) ;
      else
      n_mn = opt.ns.exp_tm_s / 60 ;
         if ( n_mn == 1 ), lbl_mn = 'min' ; else, lbl_mn = 'mins' ; end
      opt.plt.ttl = sprintf( '%4.0f-%4.0f nm, %1.1f min of observation (noiseless) %s', opt.lmbd_img_1_nm( i_bnd ), opt.lmbd_img_2_nm( i_bnd ), opt.ns.exp_tm_s / 60, opt.plt.tm_stmp ) ;
      end
    end
  else
  opt.plt.ttl = opt.plot.title ;
  end

% Super title of the scene
  if ~isfield( opt.plot, 'suptitle' )
  str_tmp = strrep( opt.scn.nm, '_', ' ' ) ;
  str_tmp( 1 ) = upper( str_tmp( 1 ) ) ;
  opt.plt.spttl = sprintf( '%s at %2.1f pc', str_tmp, opt.str.dst_pc ) ;
    if isfield( opt.str, 'nm' )
    opt.plt.spttl = sprintf( '%s, %2.1f pc', opt.str.nm, opt.str.dst_pc ) ;
    end
  else
  opt.plt.spttl = opt.plot.suptitle ;
  end

% Log10 scale
  if ~isfield( opt.plot, 'log_10' )
  opt.plt.log_10 = 0 ;
  else
  opt.plt.log_10 = opt.plot.log_10 ;
  end

% Noflip of x-axis (AU axis)
  if ~isfield( opt.plot, 'noflip' )
  opt.plt.noflp = 0 ; % N (down), E (right) http://www.ugastro.berkeley.edu/infrared10/adaptiveoptics/binary_orbit.pdf
  else
  opt.plt.noflp = opt.plot.noflip ; % N (up), E (left)
  end

% Combinations to be considered
  if ~isfield( opt.plot, 'combination_list' )
  opt.plt.cmb_lst = [ 1, 2, 3, 4 ] ;
  % if it is noiseless
    if ~( opt.ns.do )
    opt.plt.cmb_lst = [ 3 ] ;
    end
  else
  opt.plt.cmb_lst = opt.plot.combination_list ;
    if ~( opt.ns.do )
      if numel( opt.plt.cmb_lst ) == 1
        if opt.plt.cmb_lst == 1, opt.plt.cmb_lst = 3 ; end
        if opt.plt.cmb_lst == 2, opt.plt.cmb_lst = 4 ; end
      end
      if numel( opt.plt.cmb_lst ) > 1
      opt.plt.cmb_lst = [ 3, 4 ] ; % Only these 2 make sense if the simulation is noiseless.
      end
    end
  end

% With/without planet's labels
  if ~isfield( opt.plot, 'add_planet_label' )
    if ( opt.n_pl ) && ( opt.plnt.add.do )
    opt.plt.add_plnt_lbl = 1 ;
    else
    opt.plt.add_plnt_lbl = 0 ;
    end
  else
  opt.plt.add_plnt_lbl = opt.plot.add_planet_label ;
  end
