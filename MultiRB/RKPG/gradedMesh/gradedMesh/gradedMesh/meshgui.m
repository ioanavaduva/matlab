function varargout = meshgui(varargin)
% MESHGUI A dynamic mesh function generator.
%   MeshGui allows user to choose a curvature and weight paramater
%   dynamically. When the user presses the Accept button the GUI will
%   return a function for designing meshes with the [c, w] parameters 
%   baked into the function. 

%   See meshfunc for more details.
%
%   Ben Postlethwaite 2012
%   benpostlethwaite.ca


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @meshgui_OpeningFcn, ...
                   'gui_OutputFcn',  @meshgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function showmesh(handles)
% Show the mesh interactively
c = handles.c;
w = handles.w;
n = handles.n;

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Create new data
[x, dx] = handles.gmesh( linspace(0, 1 , n) );
xc = centre(x);
% Transform dx into percentages
dx = 1./(max(dx)) .* dx; 

% Plot
plot(xc, dx, 'r*')
plotmesh(x, x ./max(x) ) % plot in percents on Y axis
title(sprintf('Graded mesh with curve = %i, weight = %1.2f, ncells = %i',...
    c, w, n), 'fontSize', 14)
ylabel('1D Cell size %')
legend('Calculated dx','mesh line \rho','Location','NorthEastOutside')


% --- Executes just before meshgui is made visible.
function meshgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to meshgui (see VARARGIN)

handles.c = 4;
handles.w = 0.1;
handles.n = 30;
handles.gmesh = meshfunc(handles.c, handles.w);
showmesh(handles)

% Choose default command line output for meshgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes meshgui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = meshgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.gmesh;
delete(handles.figure1)


% --- Executes on slider movement.
function weightSlider_Callback(hObject, eventdata, handles)
% hObject    handle to weightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
w = get(hObject,'Value');
handles.w = w;
handles.gmesh = meshfunc(handles.c, handles.w);
showmesh(handles)

% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function weightSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function curveSlider_Callback(hObject, eventdata, handles)
% hObject    handle to curveSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c = ceil((get(hObject,'Value')));
if mod(c,2)
    c = c + 1;
end
handles.c = c;
handles.gmesh = meshfunc(handles.c, handles.w);
showmesh(handles)
set(handles.curveSlider,'Value',c)

% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function curveSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curveSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in OkayButton.
function OkayButton_Callback(hObject, eventdata, handles)
% hObject    handle to OkayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
uiresume(handles.figure1);


% --- Executes on slider movement.
function nslider_Callback(hObject, eventdata, handles)
% hObject    handle to nslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = round((get(hObject,'Value')));
handles.n = n;
showmesh(handles)
set(handles.nslider,'Value', n)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function nslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
