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

function efDefectImg = UTot2ImagePlane(  lambdaIn, opt, pupil, UTotL, mssg )
% Function to translate the UTot field from makeStarshadeImage to the image plane
% Author: Eric J. Cady. eric.j.cady@jpl.nasa.gov
% Modified: Sergi R. Hildebrandt. srh.jpl.caltech@gmail.com (mft_shift.m)

tic
units()
Nx_img = opt.nx_img ;
masPerPixel = opt.diam_img_mas / Nx_img ;
impeak = mft_shift( 1, Nx_img, pupil, 0, 0, 1 ) ; % Use 1 L/D square. D=diameter of the primary mirror
peak = max(abs(impeak(:)));
efDefectImg = zeros(Nx_img, Nx_img, length(lambdaIn));
for ll = 1:length(lambdaIn)
    UTot = UTotL(:,:,ll);
    lambda = lambdaIn(ll);
      if opt.verbose   
      disp(['Wavelength: ' num2str(lambda*1e9) 'nm'])
      end

    imagePlaneDiameterInMAS = Nx_img*masPerPixel;
    imagePlaneDiameterInLambdaOverD = imagePlaneDiameterInMAS*mas*opt.dmtr_tlscp_m/lambda;

    imagePlane = mft_shift(imagePlaneDiameterInLambdaOverD, Nx_img, UTot.*pupil, 0, 0, 1 ) ;
    imagePlane = imagePlane/peak; % Normalize to unocculted on-axis peak = 1

    efDefectImg(:,:,ll) = imagePlane;
end
t = toc ;
  if exist( 'mssg', 'var' )
    if ( mssg )
    disp( sprintf( 'efDefectImg took %3.2f seconds', t ) )
    end
  end
