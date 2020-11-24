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
% makePolygon
% Author: Eric J. Cady. eric.j.cady@jpl.nasa.gov
%
% Inputs:
%  - Profile: A row or column vector of points defining the shape of a
% petal as a function of radial distance.
%  - maxR: The desired outer radius of the petals.  Assumed scalar.
%  - numPetals: The number of petals on the polygon.  Assumed scalar.
%  - minFeatureSize: The smallest gap or point on a petal.  Assumed scalar.
%
% Outputs:
%  - polygonX, polygonY: a set of x- and y-values defining the polygon; each
% one is consecutive, and the last is the same as the first, to make a
% closed shape.  These will be column vectors.
%
% Limitations:
%  - This function assumes that an n-point profile is divided into n
% equal-sized regions, and the profile is defined at the midpoint of those
% regions.  This is to be consistent with our occulter optimization code.
%  - This function assumes maxR and minFeatureSize to be in the desired output units.

function [polygonX polygonY] = makePolygon(Profile, r, maxR, numPetals, minFeatureSize, a, useFullOcculter)

numLowGap = 1000;
numHighGap = 1000;
nInterpR = length(r);

% Standardize r and Profile as row vectors.
r = r(:).';
Profile = Profile(:).';

% Define polar coordinates of points on a petal.
holeTheta = 2*pi/numPetals/2*Profile;

% Truncate the petals; the 0.99 is a buffer, as the optimization doesn't
% get exactly on a value.  We include the isAt0or1 term to ensure that we
% don't get the unnecessary points (open space, solid region) even in the
% absence of truncation.
isAt0or1 = (Profile == 0) | (Profile == 1);
closeTo0 = (pi/numPetals*r.*Profile < minFeatureSize/2*.99) & (pi/numPetals*r >= minFeatureSize/2*.99);
closeTo1 = (pi/numPetals*r.*(1-Profile) < minFeatureSize/2*.99);
notClose = ~(closeTo0 | closeTo1 | isAt0or1);
rr = r(notClose);
theta = holeTheta(notClose);

if (sum(notClose) == 0)
    % It's an empty hole.  Or series of them.  Or something very, very
    % close.
    
    % If it's just an entirely open circle, we'll try to make it.  Otherwise punt--you
    % can't have a disconnected set of occulter rings anyway.
    
    switchOver = (find(Profile == 0, 1, 'first') == find(Profile == 1, 1, 'last')+1);    
    if (sum(isAt0or1) == length(isAt0or1)) && (Profile(1) == 1) && (isempty(switchOver) || switchOver == 1)
        if isempty(switchOver)
            radius = maxR;
        else
            radius = r(find(Profile == 0, 1, 'first'));
        end
        ntheta = length(Profile*numPetals);
        thetaVals = linspace(0, 2*pi, ntheta+1);
        polygonX = radius.*cos(thetaVals);
        polygonY = radius.*sin(thetaVals);
    else
       error('Error: occulter structure is too disconnected to model with makePolygon.')
    end
else
    if rr(1) ~= a
        rr = [a rr];
        theta = [theta(1) theta];
    end
    if rr(end) ~= maxR
        rr = [rr maxR];
        theta = [theta theta(end)];
    end

    tmp = linspace(rr(1), rr(end), nInterpR);
    tmp2 = interp1(rr, theta, tmp);
    rr = tmp;
    theta = tmp2;   

    ldt = 2*(pi/numPetals - theta(end))/numLowGap;
    lowgap = (pi/numPetals - theta(end) + ldt/2:ldt:pi/numPetals + theta(end) - ldt/2).';
%     lowgap = (-theta(1) + ldt/2:ldt:theta(1) - ldt/2).';
    hdt = 2*theta(end)/numHighGap;
    highgap = (-theta(end) + hdt/2:hdt:theta(end) - hdt/2).';
    
    % Rotate this petal (and its x-axis mirror) about the origin
    polygonX = [];
    polygonY = [];
    if useFullOcculter == 1  % Full occulter
        for j = 1:numPetals
            % We calculate both sides of the petal, and flip one of them, so that
            % the polygon traces a path out, then back in.
            polygonX = [polygonX; (rr.*cos(2*pi*(j-1)/numPetals-theta)).'; rr(end)*cos(2*pi*(j-1)/numPetals + highgap); ...
                fliplr(rr.*cos(2*pi*(j-1)/numPetals + theta)).'; flipud(rr(1)*cos(2*pi*(j-1)/numPetals + lowgap))];
            polygonY = [polygonY; (rr.*sin(2*pi*(j-1)/numPetals-theta)).'; rr(end)*sin(2*pi*(j-1)/numPetals + highgap); ...
                fliplr(rr.*sin(2*pi*(j-1)/numPetals + theta)).'; flipud(rr(1)*sin(2*pi*(j-1)/numPetals + lowgap))];
        end
        
        % Make the polygon closed
        polygonX = [polygonX; polygonX(1)];
        polygonY = [polygonY; polygonY(1)];       
    elseif useFullOcculter == 2  % One petal
        polygonX = [polygonX; (rr.*cos(-theta)).'; fliplr(rr.*cos(theta)).'];
        polygonY = [polygonY; (rr.*sin(-theta)).'; fliplr(rr.*sin(theta)).'];
    elseif useFullOcculter == 3 % One edge
        polygonX = [polygonX; (rr.*cos(theta)).'];
        polygonY = [polygonY; (rr.*sin(theta)).'];
    else
        error('Bad edge specification.')
    end
end

