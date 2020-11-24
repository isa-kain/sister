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
function imagePlane = mft_shift( imagePlaneDiameterInLambdaOverD, numImagePoints, inputField, x_cntr, y_cntr, px_psf_mas )
% Author: Eric Cady. eric.j.cady@jpl.nasa.gov
% Modifed by Sergi Hildebrandt (srh.jpl.caltech@gmail.com) to include a shifted center.

[x dx] = interval_mft( 1, size( inputField, 1));
[y dy] = interval_mft( 1, size( inputField, 2));
[xpp dxpp] = interval_mft( imagePlaneDiameterInLambdaOverD, numImagePoints);
[ypp dypp] = interval_mft( imagePlaneDiameterInLambdaOverD, numImagePoints);

% Take to image plane
[Xout Yout] = meshgrid(xpp, ypp);
xpp = xpp + x_cntr / numImagePoints / px_psf_mas * ( max( xpp( : ) ) - min( xpp( : ) ) ) ;
ypp = ypp + y_cntr / numImagePoints / px_psf_mas * ( max( ypp( : ) ) - min( ypp( : ) ) ) ;
imagePlane = exp(-2*pi*1i*xpp.'*x)*(inputField)*exp(-2*pi*1i*y.'*ypp)*dx*dy;

% Let's find the centroid falls within half a pixel from the center (notice the order y, x due to matlab ordering)
%[ y_msh x_msh ] = meshgrid( 1 : numImagePoints, 1 : numImagePoints ) ;
%int_src = abs( imagePlane ).^2 ;
%x_cntrd = sum( sum( x_msh .* int_src ) ) / sum( sum( int_src ) ) ;
%y_cntrd = sum( sum( y_msh .* int_src ) ) / sum( sum( int_src ) ) ;
%disp( '' )
%disp( sprintf( '%3.3f, %3.3f', 2 * x_cntrd - 1 - numImagePoints, 2 * y_cntrd - 1 - numImagePoints ) ) ;
%disp( '' )
%make_a_stop

function [x, dx] = interval_mft(width, N)
dx = width/N;
x = -width/2 + dx/2:dx:width/2-dx/2;
