% Adam Benabbou
% Fall 2019 - Summer 2020
% -------------------------------------------------------------------------

function varargout = DataVisualizer(varargin)
% DATAVISUALIZER MATLAB code for DataVisualizer.fig
%      DATAVISUALIZER, by itself, creates a new DATAVISUALIZER or raises the existing
%      singleton*.
%
%      H = DATAVISUALIZER returns the handle to a new DATAVISUALIZER or the handle to
%      the existing singleton*.
%
%      DATAVISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAVISUALIZER.M with the given input arguments.
%
%      DATAVISUALIZER('Property','Value',...) creates a new DATAVISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataVisualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataVisualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataVisualizer

% Last Modified by GUIDE v2.5 03-Jun-2020 05:14:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataVisualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @DataVisualizer_OutputFcn, ...
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


% --- Executes just before DataVisualizer is made visible.
function DataVisualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataVisualizer (see VARARGIN)

storedStructure = load('A_6-02_999Samples_Max Speed_V3.mat');
lb = storedStructure.lb;
ub = storedStructure.ub;
handles.net_mass = storedStructure.net_mass;
handles.net_disp = storedStructure.net_disp;
handles.net_stress = storedStructure.net_stress;
handles.net_CL = storedStructure.net_CL;
handles.net_CD = storedStructure.net_CD;
handles.net_CL_approx = storedStructure.net_CL_approx;

handles.FlightCondition = storedStructure.MainFlightCondition;
handles.LW_ratio = storedStructure.LW_ratio;


handles.costfunction = storedStructure.f;
handles.constraints = storedStructure.g;

%handles.CL_CD = @(x) -handles.net_CL(x)./handles.net_CD(x);

handles.f_surr_min_set = storedStructure.f_min;
handles.x_surr_min_set = storedStructure.x_min;

[handles.f_surr_min, handles.best_index] = min(handles.f_surr_min_set);

handles.x_min = handles.x_surr_min_set(:,handles.best_index); %shared by surrogate and real

handles.f_real_min = storedStructure.f_min_result;

handles.real_total_mass_opt = storedStructure.total_mass_opt;
handles.real_max_disp_opt = storedStructure.max_disp_opt;
handles.real_KS_stress_opt = storedStructure.KS_stress_opt;
handles.real_CD_opt = storedStructure.CD_opt(2,:);
handles.real_CL_opt = storedStructure.CL_opt(2,:);
handles.real_CL_approx_opt = storedStructure.CL_approx_opt;

handles.M = 100;
M=handles.M;

handles.xn1 = linspace(lb(1),ub(1),M);
handles.xn2 = linspace(lb(2),ub(2),M);
handles.xn3 = linspace(lb(3),ub(3),M);
handles.xn4 = linspace(lb(4),ub(4),M);
handles.xn5 = linspace(lb(5),ub(5),M);
handles.xn6 = linspace(lb(6),ub(6),M);

handles.x_axis = handles.xn1;
handles.y_axis = handles.xn4;


handles.plotstring = [];
handles.plotstring1= [1 0 0 0 0 0];
handles.plotstring2 = [1 0 0 0 0 0];
handles.view = '2D';
handles.title = 'Total Mass Plot';
handles.x_axis_title = 'Span (m)';
handles.y_axis_title = 'Span (m)';
handles.mode = 0;

handles.point_x1 = 1;
handles.point_x2 = 1;
%handles.real_response = -handles.real_CL_opt/handles.real_CD_opt;
handles.real_response = handles.real_total_mass_opt;

handles.fig_h = [];

handles.contour_number = 20;

handles.net = handles.net_mass;
% Choose default command line output for DataVisualizer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataVisualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataVisualizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_x_axis.
function popupmenu_x_axis_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_x_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_x_axis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_x_axis
str = get(hObject, 'String');
val = get(hObject, 'Value');

switch str{val}
    case 'Span'
        handles.x_axis = handles.xn1;
        handles.plotstring1 = [1 0 0 0 0 0];
        handles.x_axis_title = 'Span (m)';
    case 'Root Chord'
        handles.x_axis = handles.xn2;
        handles.plotstring1 = [0 1 0 0 0 0];
        handles.x_axis_title = 'Root Chord Length (m)';
    case 'Tip Chord'
        handles.x_axis = handles.xn3;
        handles.plotstring1 = [0 0 1 0 0 0];
        handles.x_axis_title = 'Tip Chord Length (m)';
    case 'Panel Thickness'
        handles.x_axis = handles.xn4;
        handles.plotstring1 = [0 0 0 1 0 0];
        handles.x_axis_title = 'Panel Thickness (m)';
    case 'Rib Thickness'
        handles.x_axis = handles.xn5;     
        handles.plotstring1 = [0 0 0 0 1 0];
        handles.x_axis_title = 'Rib Thickness(m)';
    case 'Spar Thickness'
        handles.x_axis = handles.xn6;
        handles.plotstring1 = [0 0 0 0 0 1];
        handles.x_axis_title = 'Spar Thickness (m)';
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_x_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_x_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_y_axis.
function popupmenu_y_axis_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_y_axis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_y_axis
str = get(hObject, 'String');
val = get(hObject, 'Value');

%I know i have a distinct style of coding. Don't judge. It works. 
switch str{val}
    case 'Span'
        handles.y_axis = handles.xn1;
        handles.plotstring2 = [2 0 0 0 0 0];
        handles.y_axis_title = 'Span (m)';
    case 'Root Chord'
        handles.y_axis = handles.xn2;
        handles.plotstring2 = [0 2 0 0 0 0];
        handles.y_axis_title = 'Root Chord Length (m)';        
    case 'Tip Chord'
        handles.y_axis = handles.xn3;
        handles.plotstring2 = [0 0 2 0 0 0];
        handles.y_axis_title = 'Tip Chord Length (m)';
    case 'Panel Thickness'
        handles.y_axis = handles.xn4;
        handles.plotstring2 = [0 0 0 2 0 0];
        handles.y_axis_title = 'Panel Thickness (m)';
    case 'Rib Thickness'
        handles.y_axis = handles.xn5;     
        handles.plotstring2 = [0 0 0 0 2 0];
        handles.y_axis_title = 'Rib Thickness (m)';
    case 'Spar Thickness'
        handles.y_axis = handles.xn6;
        handles.plotstring2 = [0 0 0 0 0 2];
        handles.y_axis_title = 'Spar Thickness (m)';
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_y_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_y_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_surrogate_model.
function popupmenu_surrogate_model_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_surrogate_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(hObject, 'String');
val = get(hObject, 'Value');

handles.mode = 0;
switch str{val}
    case 'Total Mass'
        handles.net = handles.net_mass;
        handles.title = 'Total Mass';
        handles.real_response = handles.real_total_mass_opt;

    case 'Max Displacement'
        handles.net = handles.net_disp;
        handles.title = 'Max Displacement';
        handles.real_response = handles.real_max_disp_opt;
        
    case 'KS Stress'
        handles.net = handles.net_stress;
        handles.title = 'KS Stress';
        handles.real_response = handles.real_KS_stress_opt;

    case 'Coefficient of Lift'
        handles.net = handles.net_CL_approx;
        handles.title = 'Lift Approximation';
        handles.real_response = handles.real_CL_approx_opt;
        
    case 'CL/CD Constraint' %DONT USE!
        handles.mode = 5;
        handles.optfunc = handles.CL_CD;
        handles.title = ' - CL to CD';
        handles.real_response = -handles.real_CL_opt/handles.real_CD_opt; 

    case 'Max Displacement Constraint'
        handles.mode = 1;
        handles.optfunc = handles.constraints;
        handles.title = 'Max Displacement Constraint';
        handles.real_response = handles.real_max_disp_opt - handles.x_min(1)/10;

    case 'KS Stress Constraint'
        handles.mode = 2;
        handles.optfunc = handles.constraints;
        handles.title = 'KS Stress Constraint';
        handles.real_response = handles.real_KS_stress_opt - 1;
        
    case 'L/W Constraint'
        handles.mode = 3;
        handles.optfunc = handles.constraints;
        handles.title = 'Lift-to-Weight Constraint';
        
        rho = handles.FlightCondition.rho;
        v = handles.FlightCondition.velocity;
        S = handles.x_min(1)*0.5*(handles.x_min(2)+handles.x_min(3));
        mass = handles.real_total_mass_opt;
        handles.real_response = -(handles.real_CL_approx_opt*(1/2)*rho*v*v*S - handles.LW_ratio*((1/0.15)*mass + mass));
end

%save handles structure
guidata(hObject,handles)


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_surrogate_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_surrogate_model


% --- Executes during object creation, after setting all properties.
function popupmenu_surrogate_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_surrogate_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Surf.
function pushbutton_Surf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Surf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold on
h0 = surf(handles.x_axis, handles.y_axis, handles.currentData);
title(handles.title);
colorbar
grid on
h1 = plot3(handles.x_surr_min_set(handles.point_x1,:), handles.x_surr_min_set(handles.point_x2,:), handles.opt_responses, 'o','MarkerEdgeColor','c','MarkerFaceColor','m','MarkerSize',6);
h2 = plot3(handles.x_min(handles.point_x1), handles.x_min(handles.point_x2), handles.best_response, 's','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',12);
h3 = plot3(handles.x_min(handles.point_x1), handles.x_min(handles.point_x2), handles.real_response, 'p','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',12);
h4 = plot3([handles.x_min(handles.point_x1) handles.x_min(handles.point_x1)], [handles.x_min(handles.point_x2) handles.x_min(handles.point_x2)], [handles.best_response, handles.real_response], '-.r', 'LineWidth', 2);

set(gca,'FontSize',14,'FontName','Minion Pro');

lgd = legend ('Response Surface','Optimum Points','Best Optimum Condition','Real Model at Optimum');
lgd.Location = 'NorthWest';

if handles.view == '2D'
    view(90,0);
    ylabel(handles.y_axis_title)
else
    view(-37.5,30);
    
    %view(0,0);
    xlabel(handles.x_axis_title)
    ylabel(handles.y_axis_title)
end
hold off
% picture_name = sprintf('%s - %s vs %s (3D)', handles.title, handles.x_axis_title, handles.y_axis_title);
 pause(3);
% saveas(hObject, picture_name, 'png')
delete(h0);delete(h1);delete(h2);delete(h3); delete(h4)


% --- Executes on button press in pushbutton_Contour.
function pushbutton_Contour_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contour(handles.x_axis, handles.y_axis, handles.currentData, handles.contour_number);
title(handles.title);
colorbar
xlabel(handles.x_axis_title)
ylabel(handles.y_axis_title)


% --- Executes on button press in PlotUpdater.
function PlotUpdater_Callback(hObject, eventdata, handles)
% hObject    handle to PlotUpdater (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

M = handles.M;

[the_x, the_y] = meshgrid(handles.x_axis,handles.y_axis);

xn = [];
handles.plotstring = handles.plotstring1 + handles.plotstring2;
for i = 1:length(handles.plotstring)
    switch handles.plotstring(i)
        case 1 
            xn = [xn reshape(the_x, M*M,1)];
            handles.view = '3D';
            handles.point_x1 = i; %setting which variable to graph
            
        case 2
            xn = [xn reshape(the_y, M*M,1)];
            handles.point_x2 = i;
            
        case 3
            xn = [xn reshape(the_y, M*M,1)];
            handles.view = '2D';
            handles.point_x1 = i;
            handles.point_x2 = i;
        otherwise
            temp = handles.x_min(i)*ones(M*M,1);
            xn = [xn temp];
    end
end
handles.opt_responses = [];
handles.best_response = [];

switch handles.mode
    case 0
        f_sim = sim(handles.net,xn');
        handles.currentData = reshape(f_sim,M,M);
        handles.opt_responses = sim(handles.net,handles.x_surr_min_set);
        handles.best_response = sim(handles.net,handles.x_min);

    case 1 % max disp constraint
        f_sim = feval(handles.optfunc,xn');
        handles.currentData = reshape(f_sim(1,:),M,M);    
        dataset = feval(handles.optfunc,handles.x_surr_min_set);
        handles.opt_responses = dataset(1,:);
        dataset2 = feval(handles.optfunc,handles.x_min);
        handles.best_response = dataset2(1); %dataset2(1,handles.best_index); Just one.
    case 2 %ks stress constraint
        f_sim = feval(handles.optfunc,xn');
        handles.currentData = reshape(f_sim(2,:),M,M);
        dataset = feval(handles.optfunc,handles.x_surr_min_set);
        handles.opt_responses = dataset(2,:);
        dataset2 = feval(handles.optfunc,handles.x_min);
        handles.best_response = dataset2(2); %dataset2(2,handles.best_index);
        
    case 3 %L/W constraint
        f_sim = feval(handles.optfunc,xn');
        handles.currentData = reshape(f_sim(3,:),M,M);
        dataset = feval(handles.optfunc,handles.x_surr_min_set);
        handles.opt_responses = dataset(3,:);
        dataset2 = feval(handles.optfunc,handles.x_min);
        handles.best_response = dataset2(3);%dataset2(2,handles.best_index);

    case 5 %CL/CD - MASS IS SAME AS JUST SURROGATE. 
        f_sim = feval(handles.optfunc,xn');
        handles.currentData = reshape(f_sim,M,M);
        handles.opt_responses = feval(handles.optfunc,handles.x_surr_min_set);
        handles.best_response = feval(handles.optfunc,handles.x_min);
    otherwise
        error('nonexistant mode');
end
guidata(hObject,handles)



function ContourNumberUpdater_Callback(hObject, eventdata, handles)
% hObject    handle to ContourNumberUpdater (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ContourNumberUpdater as text
%        str2double(get(hObject,'String')) returns contents of ContourNumberUpdater as a double
handles.contour_number = str2double(get(hObject, 'String'));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function ContourNumberUpdater_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ContourNumberUpdater (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Save.
function pushbutton_Save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

picture_name = sprintf('%s - %s vs %s (3D)', handles.title, handles.x_axis_title, handles.y_axis_title);

saveas(hObject, picture_name, 'png')
