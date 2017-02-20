function quasicrystaldemo(N,scale,angles,frames,useLogPolar)
%QUASICRYSTALDEMO  Generate quasicrystal animation via sum of plane waves.
%   QUASICRYSTALDEMO(N,SCALE,ANGLES,FRAMES,USELOGPOLAR) displays a quasicrystal
%   animation. The input N is the dimension of the animated images. N can be a
%   scalar or a two element vector, where width is N(1) and height is N(2). The
%   SCALE input specifies the number of periods per wave across the largest
%   dimension of N and ANGLES specifies the integer number of planes waves that
%   are to be summed. The magnitude of the scalar input FRAMES specifies the
%   integer number of divisions, or phase shifts, between 0 and 2*PI to render.
%   The sign of FRAMES denotes the if the phase shifts are counter-clockwise
%   (positive) or clockwise (negative). The USELOGPOLAR input is a logical value
%   that, when true, specifies output in the Log-Polar coordinate system,
%   instead of Cartesian coordinates.
%
%   Example:
%       quasicrystaldemo([640 480],32,7,30,1);
%
%   See also QUASICRYSTAL
%
%   Andrew D. Horchler, horchler @ gmail . com, Created 11-15-11
%   Revision: 1.1, 3-29-12

%   More details and original inspiration from these sources:
%
%       http://wealoneonearth.blogspot.com/search/label/quasicrystal
%       http://mainisusuallyafunction.blogspot.com/search/label/quasicrystal


% Initialization
if useLogPolar
    x = (0.5*(1-N(1)):0.5*(N(1)-1));
    y = (0.5*(1-N(end)):0.5*(N(end)-1)).';
else
    x = (0.5*(1-N(1)):0.5*(N(1)-1))*2*pi*scale/(max(N)-1);
    y = (0.5*(1-N(end)):0.5*(N(end)-1)).'*2*pi*scale/(max(N)-1);
end
x = x(ones(1,N(end)),:);
y = y(:,ones(1,N(1)));

% Apply log-polar transformation is specified
if useLogPolar
    r = -log(sqrt(x.*x+y.*y))*scale;
    x = atan2(y,x)*scale;
    y = r;
end

% Rotation angles for waves to be summed
rx = cos(0:pi/angles:pi);
ry = sin(0:pi/angles:pi);

% Phase shift counter-clockwise or clockwise
for t=2*pi/abs(frames)*[0:frames-1 0:-1:1+frames]
    qc = cos(x+t);
    for i=2:angles
        qc = qc+cos(x*rx(i)+y*ry(i)+t);
    end
    
    % Scale and output to figure
    imshow(qc*(0.5/angles)+0.5,'Parent',gca);
    drawnow;
end