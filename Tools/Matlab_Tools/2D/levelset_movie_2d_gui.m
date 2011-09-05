function varargout = levelset_movie_2d_gui(varargin)

if nargin == 0  % LAUNCH GUI
	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);

	handles.directory = '.';
	handles.number_of_frames = 0;
	while (fopen(sprintf('%s/levelset.%d',handles.directory,handles.number_of_frames + 1))) >= 0,
		handles.number_of_frames = handles.number_of_frames + 1;
	end
	[handles.m, handles.n] = get_dimensions(handles);

	guidata(fig, handles);

	set(handles.slider, 'SliderStep',[1/handles.number_of_frames 1/handles.number_of_frames]);
	set(handles.slider, 'Min', 0);
	set(handles.slider, 'Max', handles.number_of_frames);
	set(handles.text_frame, 'String', 'frame');

%	update_plot(handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


% --------------------------------------------------------------------
function [m,n] = get_dimensions(handles)

fid = fopen(sprintf('%s/header.0',handles.directory),'rb','l');
m = fread(fid,1,'int');
n = fread(fid,1,'int');

% --------------------------------------------------------------------
function update_plot(handles)

frame = round(get(handles.slider, 'Value'));
set(handles.text_frame, 'String', sprintf('%d', frame));
fid = fopen(sprintf('%s/levelset.%d',handles.directory,frame),'rb','l');
phi = fread(fid,[handles.n,handles.m],'double');
axes(handles.axes);
orient tall;
contour(phi,0,'b');
axis equal;
axis([1 handles.m 1 handles.n]);

% --------------------------------------------------------------------
function varargout = slider_Callback(h, eventdata, handles, varargin)

update_plot(handles);
