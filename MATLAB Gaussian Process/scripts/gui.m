function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 26-Aug-2018 22:29:07


% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_OpeningFcn, ...
    'gui_OutputFcn',  @gui_OutputFcn, ...
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

% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global testset trainsetall trainsetdiff trainsetsim
load('data/traindatasim.mat');
trainsetsim = Xtrain;
load('data/traindatadiffclass17.mat');
trainsetdiff = Xtrain;
load('data/traindataallclass1.mat');
trainsetall = Xtrain;
load('data/testdataclass1.mat');
testset = Xpred;
global iter;
if isempty(iter)
    iter = 1;
end
global xselected vselected;
xselected = testset{1}(1,:);
vselected = testset{1}(4,:);
global stepstopredict
stepstopredict = 40;
makeplot(handles);
% This sets up the initial plot - only do when we are invisible
% so window can get raised using gui.

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% increase
axes(handles.axes1);
cla;

global iter
iter = iter - 5;
if iter < 1
    iter = 1;
end

makeplot(handles);

function pushbutton3_Callback(hObject, eventdata, handles)
% decrease
global iter xselected
iter = iter + 5;
if iter >= size(xselected,2)
    iter = size(xselected,2) - 1;
end
makeplot(handles);


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

3
% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
    
end
% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

makeplot(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Trained on different class', 'Trained on really similiar', 'Trained on same class, almost all trajectories', 'Trained on all without time as a feature','Trained on similiar, parameters optimized'});


% --- Executes on button press in pushbutton3.


function  makeplot(handles)
popup_sel_index = get(handles.popupmenu1, 'Value');
global xselected vselected iter stepstopredict;
%1 - Trained on different class
%2 - Trained on really similiar
%3 - Trained on same class, almost all trajectories
%4 - Trained on all without time as a feature
global trainsetall trainsetdiff trainsetsim
plot(1:size(xselected,2), xselected, 'k-','LineWidth',1.5)
hold on;
h3 = plot(iter, xselected(iter), 'ko');
switch popup_sel_index
    case 1
        for i=1:size(trainsetdiff,2)
            h5 = plot(trainsetdiff{i}(1,:),'Color',[222 222 222]/255);
        end
        [xpred xsd] = differentclass(xselected,vselected,iter,stepstopredict);
    case 2
        for i=1:size(trainsetsim,2)
            h5 = plot(trainsetsim{i}(1,:),'Color',[222 222 222]/255);
        end
        [xpred xsd] = similiartrajectories(xselected,vselected,iter,stepstopredict);
    case 3
        for i=1:size(trainsetall,2)
            h5 = plot(trainsetall{i}(1,:),'Color',[222 222 222]/255);
        end
        [xpred xsd] = sameclassalltrajectories(xselected,vselected,iter,stepstopredict);
    case 4
        for i=1:size(trainsetall,2)
            h5 = plot(trainsetall{i}(1,:),'Color',[222 222 222]/255);
        end
        [xpred xsd] = withouttime(xselected,vselected,iter,stepstopredict);
    case 5
        for i=1:size(trainsetsim,2)
            h5 = plot(trainsetsim{i}(1,:),'Color',[222 222 222]/255);
        end
        tic
        [xpred xsd] = simoptimized(xselected,vselected,iter,stepstopredict);
        temptotal = toc
end
h0 = plot(1:size(xselected,2), xselected, 'k-','LineWidth',1.5);


grid on;
tplot = iter+1:size(xpred,2)+iter;
h1 = plot(tplot, xpred, 'b-','LineWidth',1.5);
h2 = plot(tplot, xpred + 2*xsd, 'g-','LineWidth',1.5);
plot(tplot, xpred - 2*xsd, 'g-','LineWidth',1.5);
legend([h0 h1 h2 h3 h5],'True trajectory','Predicted trajectory','2 standard deviations above and below'...
    ,'Current point','Training trajectories');
if popup_sel_index ~=1
    ylim([min(xselected)-0.07 max(xselected)+0.07]);
end
ylabel('Position(m)')
xlabel('Timestep')
hold off
%     plot(tmarked, xmarked, 'ko');



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iter
iter = [];


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
popup_sel_index = get(handles.popupmenu2, 'Value');
global xselected vselected testset;
xselected = testset{popup_sel_index}(1,:);
vselected = testset{popup_sel_index}(4,:);
global iter
if iter >= size(xselected,2)
    iter = size(xselected,2) - 1;
end
makeplot(handles);

function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global testset;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
names = {};
for i = 1:size(testset,2)
    names = [names; num2str(i)];
end
set(hObject, 'String', names)
hObject



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global stepstopredict
stepstopredict = str2num(get(handles.edit1, 'String'));
makeplot(handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
