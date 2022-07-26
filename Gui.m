function varargout = Gui(varargin)
% GUI MATLAB code for Gui.fig
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
%      applied to the GUI before Gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gui

% Last Modified by GUIDE v2.5 20-Mar-2021 17:18:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_OutputFcn, ...
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


% --- Executes just before Gui is made visible.
function Gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gui (see VARARGIN)

% Choose default command line output for Gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global N detener WF X1 X2 Y ;


% UIWAIT makes Gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Entrenar.
function Entrenar_Callback(hObject, eventdata, handles)

global N detener WF range1 range2 minX1 minX2 X1 X2 Y; 
detener=false;

axes(handles.axes1); % Make averSpec the current axes.
cla reset; % Do a complete and total reset of the axes.

axes(handles.axes2); % Make averSpec the current axes.
cla reset; % Do a complete and total reset of the axes.

if get(handles.radiobutton3,'Value')==1
datos1=str2double(get(handles.uitable1,'Data'));

X1=datos1(:,1);
X2=datos1(:,2);
Y=datos1(:,3);
end

range1 = max(X1(:)) - min(X1(:));
Xn01 = (X1 - min(X1(:))) / range1;
X1n = 2 * Xn01 - 1

minX1=min(X1(:));
minX2=min(X2(:));

range2 = max(X2(:)) - min(X2(:));
Xn02 = (X2 - min(X2(:))) / range2;
X2n = 2 * Xn02 - 1

Xc=[];
Xx=[];

for i=1:N 
    if Y(i)==1
   Xc=vertcat(Xc,[X1n(i),X2n(i)])
    else
   Xx=vertcat(Xx,[X1n(i),X2n(i)])
    end
end

 plot(handles.axes1, Xc(:,1), Xc(:,2), 'o') 
 hold(handles.axes1, 'on' )
 plot(handles.axes1, Xx(:,1), Xx(:,2), 'x')

 xlim(handles.axes1,[-2 2])
 ylim(handles.axes1,[-2 2])
 
 Wk=[(-1+2*rand);(-1+2*rand);(-1+2*rand)]
 %coeff=rand/10;
 coeff=0.1;
 error=1;
 ep=[];
 errorT=[];
 kepT=[];
 kep=0;
 WkT=[];

 if get(handles.radiobutton1,'Value')==1
   while abs(error) > 0.01
       for i=1:N         
           Xk=[1;X1n(i);X2n(i)];
           S=dot(Xk,Wk);
           Phi=tanh(S);
           ek=Y(i)-Phi;
           Phid=1-(Phi^2);
           Wk1=Wk+coeff*ek*Phid*Xk;
           Wk=Wk1;
           ep = ep + ek^2;
           
           axes(handles.axes1); % Make averSpec the current axes.
           cla reset; % Do a complete and total reset of the axes.

           axes(handles.axes2); % Make averSpec the current axes.
           cla reset; % Do a complete and total reset of the axes.

           Xg2=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3));
           plot(handles.axes1,X1n,Xg2);
           
           hold(handles.axes1, 'on' );
           Xg22=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3))-(0.9/Wk(3));
           plot(handles.axes1,X1n,Xg22);
           
           hold(handles.axes1, 'on' );
           Xg23=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3))+(0.9/Wk(3));
           plot(handles.axes1,X1n,Xg23);
           
           hold(handles.axes1, 'on' );
           plot(handles.axes1, Xc(:,1), Xc(:,2), 'o') ;
           hold(handles.axes1, 'on' );
           plot(handles.axes1, Xx(:,1), Xx(:,2), 'x');
           
           xlim(handles.axes1,[-2 2]);
           ylim(handles.axes1,[-2 2]);
           
       end
       kep=kep+1;
       
       error=0.5*ep;
       errorT=horzcat(errorT,error);
       kepT=horzcat(kepT,kep);
       ep=[];
       
       plot(handles.axes2,kepT,errorT);
       xlim(handles.axes2,[0 kep]); 
       
   end
   WF=Wk;
 end

 
  if get(handles.radiobutton2,'Value')==1
  ep=0;
   WkT=zeros(3,1)
   while abs(error) > 0.01
       for i=1:N 
           Xk=[1;X1n(i);X2n(i)];
           S=dot(Xk,Wk)
           Phi=tanh(S);
           ek=Y(i)-Phi
           Phid=1-(Phi^2)
           Wk1=coeff*ek*Phid*Xk
           WkT=WkT+Wk1         
           ep=ep+ek^2
       end
       Wk=Wk+WkT;
       WkF=Wk;
       
           axes(handles.axes1); % Make averSpec the current axes.
           cla reset; % Do a complete and total reset of the axes.

           axes(handles.axes2); % Make averSpec the current axes.
           cla reset; % Do a complete and total reset of the axes.

           Xg2=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3));
           plot(handles.axes1,X1n,Xg2);
           
           hold(handles.axes1, 'on' );
           Xg22=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3))-(0.9/Wk(3));
           plot(handles.axes1,X1n,Xg22);
           
           hold(handles.axes1, 'on' );
           Xg23=-(Wk(2)/Wk(3))*X1n-(Wk(1)/Wk(3))+(0.9/Wk(3));
           plot(handles.axes1,X1n,Xg23);
           
           hold(handles.axes1, 'on' );
           plot(handles.axes1, Xc(:,1), Xc(:,2), 'o') ;
           hold(handles.axes1, 'on' );
           plot(handles.axes1, Xx(:,1), Xx(:,2), 'x');
           
           xlim(handles.axes1,[-2 2]);
           ylim(handles.axes1,[-2 2]);
        
       kep=kep+1;
       
       error=0.5*ep;
       errorT=horzcat(errorT,error);
       kepT=horzcat(kepT,kep);
       
       
       plot(handles.axes2,kepT,errorT);
       xlim(handles.axes2,[0 kep]); 
       
       WkT=[];
       ep=0;
   end
 WF=WkF
end
 

 

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

global detener;

detener=true



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
    set(handles.radiobutton2,'Value',0)




% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
    set(handles.radiobutton1,'Value',0)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global WF range1 range2 minX1 minX2;

WF

X1p=str2double(get(handles.edit2, 'String'));
X2p=str2double(get(handles.edit3, 'String'));

Xn01 = (X1p - minX1) / range1;
X1n = 2 * Xn01 - 1

Xn02 = (X2p - minX2) / range2;
X2n = 2 * Xn02 - 1

Xk=[1;X1n;X2n]

S=dot(Xk,WF);
Phi=tanh(S);

if Phi>0.8
set(handles.edit4,'string','clase 1'); 
elseif Phi<-0.8
    set(handles.edit4,'string','clase 0');    
else
  set(handles.edit4,'string',num2str(Phi));     
end


    
    



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Asignar.
function Asignar_Callback(hObject, eventdata, handles)

global N X1 X2 Y; 

if get(handles.radiobutton3,'Value')==1 
    N=str2double(get(handles.edit5, 'String'));
    num_elem=cell(N,3);
    num_elem(:,:)={''};
    set(handles.uitable1,'Data',num_elem)
    set(handles.uitable1,'ColumnEditable',true(1,3))
end


if get(handles.radiobutton4,'Value')==1
    
Nx=str2double(get(handles.edit6, 'String'));
Nc=str2double(get(handles.edit7, 'String'));


% define maximum number of allowable mouse clicks to avoid infinite loop
% and also to preallocate array to hold click data
maxPoints = Nc;

% define plot limits so it doesn't keep rescaling which might be confusing
xLimits = [-5 5];
yLimits = [-5 5];

% open new blank figure with defined limits
figure('Name','Marque su primera clase');
xlim(xLimits)
ylim(yLimits)
hold on

%

% instruct user on how to enter points and how to terminate
disp('Click the mouse wherever in the figure; press ENTER when finished.');

% preallocate array to hold mouse click coordinates
mousePointCoords = zeros(maxPoints,2);

% set up loop to collect and display mouse click points
count = 0;
for k = 1:maxPoints
    % get the mouse click point, or terminate if user presses enter
    %  in which case the coordinates will be returned empty
    coords = ginput(1);
    if isempty(coords)
        break
    end

    count = count + 1;
    mousePointCoords(count,:) = coords;
    plot(coords(:,1),coords(:,2),'o','MarkerSize',8);
    
    
end
% clean up
hold off
mousePointCoords = mousePointCoords(1:count,:); % trim off unused array

% display results
disp(mousePointCoords)

X1=mousePointCoords(:,1);
X2=mousePointCoords(:,2);
Y=ones(Nc,1);

% define maximum number of allowable mouse clicks to avoid infinite loop
% and also to preallocate array to hold click data
close 

maxPoints = Nx;

% define plot limits so it doesn't keep rescaling which might be confusing
xLimits = [-5 5];
yLimits = [-5 5];

% open new blank figure with defined limits
figure('Name','Marque su segunda clase');
xlim(xLimits)
ylim(yLimits)
hold on

%

% instruct user on how to enter points and how to terminate
disp('Click the mouse wherever in the figure; press ENTER when finished.');

% preallocate array to hold mouse click coordinates
mousePointCoords = zeros(maxPoints,2);

% set up loop to collect and display mouse click points
count = 0;
for k = 1:maxPoints
    % get the mouse click point, or terminate if user presses enter
    %  in which case the coordinates will be returned empty
    coords = ginput(1);
    if isempty(coords)
        break
    end

    count = count + 1;
    mousePointCoords(count,:) = coords;
    plot(coords(:,1),coords(:,2),'x','MarkerSize',8);
    
    
end
% clean up
hold off
mousePointCoords = mousePointCoords(1:count,:); % trim off unused array

% display results
disp(mousePointCoords)

X1=vertcat(X1,mousePointCoords(:,1));
X2=vertcat(X2,mousePointCoords(:,2));
Y=vertcat(Y,ones(Nx,1).*-1);
close 
N=Nx+Nc;
end






% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton4,'Value',0)

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
set(handles.radiobutton3,'Value',0)



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
