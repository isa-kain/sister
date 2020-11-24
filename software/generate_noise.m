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

function scn_ns_e = generate_noise( scn_dt_cnv_e, opt )
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

sz_cnv_1 = size( scn_dt_cnv_e, 1 ) ;
sz_cnv_2 = size( scn_dt_cnv_e, 2 ) ;

% Detector's gain affects noise linearly unless pixels get saturated. Here there's no modelling about saturation.

% Shot noise:
tic_ns = tic ;
% (shot noise is ~sqrt(signal)~sqrt(gain), because the signal comes here already multiplied by the gain, so that we need to add a further factor for the gain)
sht_ns_e = round( sqrt( opt.ns.gn ) * poissnormal_rnd_2d( scn_dt_cnv_e, opt.ns.lmbd_nrml ) ) ;

% Dark current: the random part. That is, the shot noise associated with it. The constant part can be removed after being measured. Assuming a Poisson distribution, although in practice each detector has a different distribution.
drk_ns_e = round( opt.ns.gn * poissnormal_rnd_2d( opt.ns.drk_e_s * opt.ns.exp_tm_s * ones( sz_cnv_1, sz_cnv_2 ), opt.ns.lmbd_nrml ) ) ;
% Clock Induced Charge
cic_ns_e = round( opt.ns.gn * opt.ns.cic_e * sqrt( opt.ns.n_frm ) * randn( sz_cnv_1, sz_cnv_2 ) ) ;
% Read noise (does not get affected by the gain)
rd_ns_e = round( opt.ns.rd_e * sqrt( opt.ns.n_frm ) * randn( sz_cnv_1, sz_cnv_2 ) ) ;
% Optional message with timing
  if ( opt.verbose )
  disp( sprintf( 'The noise generaton took %1.3f s', toc( tic_ns ) ) )
  end

% Detector's gain affects noise linearly as a first approximation (shot noise is ~sqrt(signal)~sqrt(gain) )
scn_ns_e =  ( sht_ns_e + drk_ns_e + cic_ns_e + rd_ns_e ) ;
