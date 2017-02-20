function varargout=quasicrystal(varargin)
%QUASICRYSTAL  Generate movie or image of a quasicrystal via sum of plane waves.
%   QUASICRYSTAL(FILENAME,N,SCALE,ANGLES,FRAMES) produces an AVI movie, animated
%   GIF image, or PNG still image depending on the file extension of FILENAME,
%   '.avi', '.gif', or '.png', respectively. The input N is the dimension of the
%   movie or image. N can be a scalar or a two element vector, where width is
%   N(1) and height is N(2). The SCALE input specifies the number of periods per
%   wave across the largest dimension of N and ANGLES specifies the integer
%   number of planes waves that are to be summed. In the context of a movie or
%   animatd GIF, the magnitude of the scalar input FRAMES specifies the integer
%   number of divisions, or phase shifts, between 0 and 2*PI to render. The sign
%   of FRAMES denotes the if the phase shifts are counter-clockwise (positive)
%   or clockwise (negative). For a PNG still image, FRAMES may be any real
%   scalar value specifying a phase shift in radians relative to zero.
%   
%   QUASICRYSTAL(FILENAME,N,SCALE,ANGLES,FRAMES,...,COORDS,...) additionally
%   specifies the coordinates, COORDS, used, which are either 'Log-Polar' or the
%   default 'Cartesian'.
%   
%   QUASICRYSTAL(FILENAME,N,SCALE,ANGLES,FRAMES,...,CNTRST,...) specifies an
%   optional contrast adjustmest, CNTRST, to be applied to the colormap of the
%   output movie or image. CNTRST is a real numeric value in the range [-1, 1].
%   Positive and negative values of CNTRST signify contrast increases and
%   decreases, respectively, relative to the default of zero.
%   
%   QUASICRYSTAL(FILENAME,N,SCALE,ANGLES,FRAMES,...,CMAP,...) specifies a
%   colormap, CMAP, to use instead of the default gray(256) grayscale colormap.
%   CMAP is an Nx3 real numeric matrix, where 2 <= N <= 256. The datatype of
%   CMAP may be logical, uint8, single, or double. The values of CMAP must be in
%   the range [0, N-1] in the case of uint8 and in the range [0, 1] for single
%   and double colormaps.
%   
%   IMG = QUASICRYSTAL(N,SCALE,ANGLES,FRAMES,...) will output image data to the
%   matrix IMG without creating a file. The output IMG is transformed by any
%   applied colormap. The COORDS, CNTRST, and CMAP inputs are optional as above
%   and can be applied in any order.
%   
%   [IMG CMAP] = QUASICRYSTAL(N,SCALE,ANGLES,FRAMES,...) will output image data
%   to the matrix IMG and an assosciated colormap, CMAP, without creating a
%   file. The data in IMG is not transformed by the the colormap before output.
%   
%   Examples:
%   	quasicrystal('out.gif',512,32,7,30);
%       
%   	quasicrystal('out.avi',[640 480],32,5,30,'Log-Polar',jet(256));
%       
%       n = 30;
%    	for i=0:n-1
%       	img = quasicrystal(256,32,7,2*pi*i/n,0.9);
%        	imshow(img);
%         	drawnow;
%     	end
%
%   Andrew D. Horchler, horchler @ gmail . com, Created 11-10-11
%   Revision: 1.8, 4-16-16

%   More details and original inspiration from these sources:
%
%       http://wealoneonearth.blogspot.com/search/label/quasicrystal
%       http://mainisusuallyafunction.blogspot.com/search/label/quasicrystal


% Check inputs
if ischar(varargin{1})
    filename = varargin{1};

    % Complex handling of filename and path
    if ~ischar(filename) || isempty(filename)
        error('quasicrystal:InvalidFileName',...
              'The file name must be a non-empty string.');
    end

    filename = regexprep(filename,'[\/\\]',filesep);
    [pathString,baseFile,extensionProvided] = fileparts(filename);
    
    if isempty(pathString)
        pathString = pwd;
    end
    if exist(pathString,'dir') ~= 7
        error('quasicrystal:NonExistantFolder',...
              'The specified folder, %s, does not exist.',pathString);
    end

    % Extension determines if a move or image is output
    if ~any(strcmpi(extensionProvided,{'.avi','.gif','.png'}))
        error('quasicrystal:InavalidExtension',...
             ['The file name must have a ''.avi'', ''.gif'', or ''.png'' '...
              'extension.']);
    end
    if any(strcmpi(extensionProvided,{'.avi','.gif'}))
        isVideo = true;
        isGIF = strcmpi(extensionProvided,'.gif');
    else
        isVideo = false;
    end

    % Check if the file can be created
    filename = fullfile(pathString,[baseFile extensionProvided]);
    fileExisted = (exist(filename,'file') ~= 0);
    [fid,fidMessage] = fopen(filename,'a');

    if fid ~= -1
        [success,info] = fileattrib(filename);  %#ok<*ASGLU>
        if usejava('jvm')
            jf = java.io.File(info.Name);
            filename = char(jf.getCanonicalPath());
        else
            filename = info.Name;
        end
        fclose(fid);
    end

    % Delete file if it already exists
    if ~fileExisted && (fid ~= -1)
        delete(filename);
    end

    if fid == -1
        error('quasicrystal:FileCreation',...
              'Cannot create file %s:\n\n%s',filename,fidMessage);
    end
    
    narginchk(5,8);
    nargoutchk(0,0);
    fileOutput=1;
else
    narginchk(4,8);
    nargoutchk(0,2);
    isVideo = false;
    fileOutput = 0;
end

% Handle other inputs and shift if a file is being output
N = varargin{1+fileOutput};
if ~isvector(N) || isempty(N) || ~isnumeric(N) || ~isreal(N) || length(N) > 2
    error('quasicrystal:InvalidN',...
          'Dimension input must be one or two element numeric vector.');
end
if N(1) < 1 || mod(N(1),1) || (length(N) > 2 && (N(2) < 1 || mod(N(2),1)))
    error('quasicrystal:InvalidNValue','Dimension input must be vector >= 1.');
end

scale = varargin{2+fileOutput};
if ~isscalar(scale) || isempty(scale) || ~isnumeric(scale) || scale < 0
    error('quasicrystal:InvalidScale',...
          'Scale input must be positive scalar numeric value.');
end

angles = varargin{3+fileOutput};
if ~isscalar(angles) || isempty(angles) || ~isnumeric(angles) || ~isreal(angles)
    error('quasicrystal:InvalidAngles',...
          'Angles input must be real scalar numeric value.');
end
if angles < 1 || mod(angles,1)
    error('quasicrystal:InvalidAnglesValue',...
          'Angles input must be integer >= 1.');
end
angi = 1/angles;

frames = varargin{4+fileOutput};
if ~isscalar(frames) || isempty(frames) || ~isnumeric(frames)
    error('quasicrystal:InvalidFrames',...
          'Frames input must be scalar numeric value.');
end
if isVideo && (mod(frames,1) || frames == 0)
    error('quasicrystal:InvalidFramesValue',...
          'Frames input must be integer >= 1 or <= -1.');
end

% Default values
c = 0;
cmap = [];
nc = 255;
nci = 1;
useLogPolar = false;
isGrayscale = true;
isLinearGrayscale = true;

% Handle coordinates, contrast, and colormap input an any order
for i=5+fileOutput:nargin
    arg = varargin{i};
    if ischar(arg)
        if ~any(strcmpi(arg,{'Cartesian','Rectilinear','Log-Polar','LogPolar'}))
            error('quasicrystal:InvalidCoords',...
                  'Coords input must be ''Cartesian'' or ''Log-Polar''.');
        end
        useLogPolar = any(strcmpi(arg,{'Log-Polar','LogPolar'}));
    elseif isvector(arg)
        c = arg;
        if length(c) ~= 1 && length(c) ~= 3
            error('quasicrystal:InvalidContrastLength',...
                  'Contrast input must be one or three element vector.');
        end
        if isempty(c) || ~isnumeric(c) || ~isreal(c)
            error('quasicrystal:InvalidContrast',...
                  'Contrast input must be real vector numeric value.');
        end
        if any(abs(c)) > 1
            error('quasicrystal:InavlidContrastRange',...
                  'Contrast input must be in the range [-1 1].');
        end
    elseif ~isvector(arg) && ndims(arg) == 2    %#ok<ISMAT>
        if isempty(arg) || ~isreal(arg) || ~(isnumeric(arg) || islogical(arg))
            error('quasicrystal:InvalidColormap',...
                  'Colormap must be non-empty real numeric or logical matrix.');
        end
        if ~isa(arg,'uint8') && ~isfloat(arg) && ~islogical(arg)
            error('quasicrystal:InvalidColormapType',...
                 ['Colormap must be logical, uint8, double, or single ' ...
                  'datatype.']);
        end
        
        nc = size(arg,1)-1;
        if nc == 0 || nc > 255 || size(arg,2) ~= 3
            error('quasicrystal:InvalidColormapSize',...
                  'Colormap input must be Nx3 matrix with 2 <= N <= 256.');
        end

        if isa(arg,'uint8')
            if max(max(arg)) > nc
                error('quasicrystal:InvalidColormapRangeUint8',...
                      'Colormap must be in range [0, N-1] for uint8.');
            end
            cmap = double(arg)*(1/255);
        elseif islogical(arg)
            cmap = double(arg);
        elseif max(max(arg)) > 1 || min(min(arg)) < 0
            error('quasicrystal:InvalidColormapRange',...
                  'Colormap must be in range [0, 1] for %s.',class(arg));
        else
            cmap = arg;
        end
    else 
        error('quasicrystal:UnrecognizedInput',...
              'Invalid or unrecognized input: argument %d.',i);
    end
end

% Adjust contrast and apply to colormap
if any(c ~= 0)
    if isempty(cmap)
        cmap = gray(256);
    end
    if length(c) == 1
        cin = [0 1]+0.5*max(c,0)*[1 -1];
        cout = [0 1]-0.5*min(c,0)*[1 -1];
        cr = (cout(2)-cout(1))/(cin(2)-cin(1));
        cmap = (max(cin(1),min(cin(2),cmap))-cin(1))*cr+cout(1);
    else
        cin1 = 0.5*max(c,0);
        cin2 = 1-cin1;
        cout1 = -0.5*min(c,0);
        cr = (1-2*cout1)./(cin2-cin1);
        cmap(:,1) = (max(cin1(1),min(cin2(1),cmap(:,1)))-cin1(1))*cr(1)+cout1(1);
        cmap(:,2) = (max(cin1(2),min(cin2(2),cmap(:,2)))-cin1(2))*cr(2)+cout1(2);
        cmap(:,3) = (max(cin1(3),min(cin2(3),cmap(:,3)))-cin1(3))*cr(3)+cout1(3);
    end
end

% Detect grayscale and linear grayscale colormaps
if ~isempty(cmap) && ((isVideo && ~isGIF) || fileOutput == 0)
    if all(cmap(:,1) == cmap(:,2)) && all(cmap(:,1) == cmap(:,3))
        isGrayscale = true;
        if all(abs(diff(cmap(:,1))+cmap(1,1)-cmap(2,1))<eps)
            isLinearGrayscale = true;
            nci = 255/nc;
        else
            isLinearGrayscale = false;
        end
    else
        isGrayscale = false;
    end
end
nc = 0.5*nc;

% Initialization
if length(N) == 2 && N(1) ~= N(2)
    if useLogPolar
        x = (0.5*(1-N(1)):0.5*(N(1)-1));
        y = (0.5*(1-N(2)):0.5*(N(2)-1)).';
    else
        sc = 2*pi*scale/(max(N)-1);
        x = (0.5*(1-N(1)):0.5*(N(1)-1))*sc;
        y = (0.5*(1-N(2)):0.5*(N(2)-1)).'*sc;
    end
    x = x(ones(1,N(2)),:);
    y = y(:,ones(1,N(1)));
else
    if useLogPolar
        x = (0.5*(1-N(1)):0.5*(N(1)-1));
    else
        x = (0.5*(1-N(1)):0.5*(N(1)-1))*(2*pi*scale/(N(1)-1));
    end
    x = x(ones(1,N(1)),:);
    y = x.';
end

% Apply log-polar transformation is specified
if useLogPolar
    r = -log(sqrt(x.*x+y.*y))*scale;
    x = atan2(y,x)*scale;
    y = r;
end

% Rotation angles for waves to be summed
rx = cos(0:pi*angi:pi);
ry = sin(0:pi*angi:pi);

% Output video or image
if isVideo
    % NTSC standard frame rate
    fps = 29.97;
    if isGIF
        fmt = imformats('gif');
        wgifa = @(img,varargin)feval(fmt.write,img,cmap,filename,varargin{:});
        
        % Phase shift counter-clockwise or clockwise
        for t = 2*pi/abs(frames)*[0:frames-1 0:-1:1+frames]
            qc = cos(x+t);
            for i=2:angles
                qc = qc+cos(x*rx(i)+y*ry(i)+t);
            end
            
            % Scale and output animated GIF frame
            qc = uint8(qc*(nc*angi)+nc);
            if t == 0
                feval(wgifa,qc,'delaytime',1/fps,'loopcount',Inf);
            else
                feval(wgifa,qc,'delaytime',1/fps,'writemode','append');
            end
        end
    else
        % Handle older versions of Matlab that don't have VideoWriter class
        if verLessThan('matlab','7.11') || exist('VideoWriter','class') ~= 8
            vid = avifile(filename);    %#ok<DAVIFL>
            vid.Colormap = cmap;
            vid.Compression = 'None';
            vid.FPS = fps;

            % Phase shift counter-clockwise or clockwise
            for t=2*pi/abs(frames)*[0:frames-1 0:-1:1+frames]
                qc = cos(x+t);
                for i=2:angles
                    qc = qc+cos(x*rx(i)+y*ry(i)+t);
                end
                
                % Scale and output movie frame
                vid = addframe(vid,uint8(qc*(nc*angi)+nc));
            end
            vid = close(vid); %#ok<*NASGU>
        else
            vid = VideoWriter(filename,'Motion JPEG AVI'); 
            vid.FrameRate = fps;
            vid.Quality = 100;

            open(vid);

            if ~isGrayscale
                rgbout = zeros(N(end),N(1),3);
                rgb = zeros(N(end),N(1));
            end

            % Phase shift counter-clockwise or clockwise
            for t=2*pi/abs(frames)*[0:frames-1 0:-1:1+frames]
                qc = cos(x+t);
                for i=2:angles
                    qc = qc+cos(x*rx(i)+y*ry(i)+t);
                end
                
                % Scale and output movie frame
                if isGrayscale
                    if isLinearGrayscale
                        writeVideo(vid,nci*uint8(qc*(nc*angi)+nc));
                    else
                        qc(:) = cmap(uint8(qc*(nc*angi)+nc)+1,1);
                        writeVideo(vid,qc);
                    end
                else
                    % Apply colormap
                    qc = uint8(qc*(nc*angi)+nc)+1;
                    rgb(:) = cmap(qc,1);
                    rgbout(:,:,1) = rgb;
                    rgb(:) = cmap(qc,2);
                    rgbout(:,:,2) = rgb;
                    rgb(:) = cmap(qc,3);
                    rgbout(:,:,3) = rgb;

                    writeVideo(vid,rgbout);
                end
            end
            close(vid);
        end
    end
else
    % Just one image/frame/phase-shift
    qc = cos(x+frames);
    for i=2:angles
        qc = qc+cos(x*rx(i)+y*ry(i)+frames);
    end
    
    % Scale and save to file or output matrix
    if fileOutput
        fmt = imformats('png');
        fmt.write(uint8(qc*(nc*angi)+nc),cmap,filename);
    else
        if nargout == 1
            if isGrayscale
                if isLinearGrayscale
                    varargout{1} = nci*uint8(qc*(nc*angi)+nc);
                else
                    qc(:) = cmap(uint8(qc*(nc*angi)+nc)+1,1);
                    varargout{1} = qc;
                end
            else
                rgbout = zeros(N(end),N(1),3);
                rgb = zeros(N(end),N(1));
                
                % Apply colormap
                qc = uint8(qc*(nc*angi)+nc)+1;
                rgb(:) = cmap(qc,1);
                rgbout(:,:,1) = rgb;
                rgb(:) = cmap(qc,2);
                rgbout(:,:,2) = rgb;
                rgb(:) = cmap(qc,3);
                rgbout(:,:,3) = rgb;
                
                varargout{1} = rgbout;
            end
        else
            varargout{1} = nci*(uint8(qc*(nc*angi)+nc));
            varargout{2} = cmap;
        end
    end
end