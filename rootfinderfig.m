function varargout = rootfinderfig(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rootfinderfig_OpeningFcn, ...
                   'gui_OutputFcn',  @rootfinderfig_OutputFcn, ...
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
end

function rootfinderfig_OpeningFcn(hObject, eventdata, handles, varargin)
    set(handles.xi,'visible','off');
    set(handles.xii,'visible','off');
    set(handles.X_initial_1,'visible','off');
    set(handles.xinit_1,'visible','off');
    set(handles.xu,'visible','off');
    set(handles.xuu,'visible','off');
    set(handles.xl,'visible','off');
    set(handles.xll,'visible','off');
    set(handles.mtext1,'visible','off');
    set(handles.mtext2,'visible','off');
    set(handles.func,'visible','off');
    set(handles.ftext,'visible','off');
    set(handles.fpath,'visible','off');
    set(handles.maxtext,'visible','off');
    set(handles.errtext,'visible','off');
    set(handles.maxi,'visible','off');
    set(handles.es,'visible','off');
    set(handles.magictext,'visible','off');
    set(handles.mgcminp,'visible','off');
    set(handles.mgcfinp,'visible','off');
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.gx,'visible','off');
    set(handles.funccalc,'visible','off');

    global stepiterations;
    stepiterations=1;


    % Choose default command line output for rootfinderfig
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes rootfinderfig wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

function varargout = rootfinderfig_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
end

function number_of_iterations_Callback(hObject, eventdata, handles)

end

function number_of_iterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function execution_time_Callback(hObject, eventdata, handles)
end

function execution_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function Approximate_root_Callback(hObject, eventdata, handles)
end

function Approximate_root_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function precision_Callback(hObject, eventdata, handles)
end

function precision_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function hp_Callback(hObject, eventdata, handles)
end

function hp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function coff_Callback(hObject, eventdata, handles)
end

function coff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function bisection_Callback(hObject, eventdata, handles)
if get(handles.bisection,'value')==1
    set(handles.newton_raphson,'value',0);
    set(handles.fixed,'value',0);
    set(handles.secant,'value',0);
    set(handles.false,'value',0);
    set(handles.xu,'visible','on');
    set(handles.xuu,'visible','on');
    set(handles.xl,'visible','on');
    set(handles.xll,'visible','on');
    set(handles.xi,'visible','off');
    set(handles.xii,'visible','off');
    set(handles.X_initial_1,'visible','off');
    set(handles.xinit_1,'visible','off');
    set(handles.maxtext,'visible','on');
    set(handles.errtext,'visible','on');
    set(handles.maxi,'visible','on');
    set(handles.es,'visible','on');,
    set(handles.magictext,'visible','off');
    set(handles.mgcminp,'visible','off');
    set(handles.mgcfinp,'visible','off');
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.gx,'visible','off');
    set(handles.funccalc,'visible','off');
end 
end

function false_Callback(hObject, eventdata, handles)
if get(handles.false,'value')==1
    set(handles.newton_raphson,'value',0);
    set(handles.fixed,'value',0);
    set(handles.secant,'value',0);
    set(handles.bisection,'value',0);
    set(handles.xu,'visible','on');
    set(handles.xuu,'visible','on');
    set(handles.xl,'visible','on');
    set(handles.xll,'visible','on');
    set(handles.xi,'visible','off');
    set(handles.xii,'visible','off');
    set(handles.X_initial_1,'visible','off');
    set(handles.xinit_1,'visible','off');
    set(handles.maxtext,'visible','on');
    set(handles.errtext,'visible','on');
    set(handles.maxi,'visible','on');
    set(handles.es,'visible','on');
    set(handles.magictext,'visible','off');
    set(handles.mgcminp,'visible','off');
    set(handles.mgcfinp,'visible','off');
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.gx,'visible','off');
    set(handles.funccalc,'visible','off');
end 
end

function fixed_Callback(hObject, eventdata, handles)
if get(handles.fixed,'value')==1
    set(handles.newton_raphson,'value',0);
    set(handles.false,'value',0);
    set(handles.secant,'value',0);
    set(handles.bisection,'value',0);
    set(handles.xu,'visible','off');
    set(handles.xuu,'visible','off');
    set(handles.xl,'visible','off');
    set(handles.xll,'visible','off');
    set(handles.xi,'visible','on');
    set(handles.xii,'visible','on');
    set(handles.X_initial_1,'visible','off');
    set(handles.xinit_1,'visible','off');
    set(handles.maxtext,'visible','on');
    set(handles.errtext,'visible','on');
    set(handles.maxi,'visible','on');
    set(handles.es,'visible','on');
    set(handles.magictext,'visible','on');
    set(handles.mgcminp,'visible','on');
    set(handles.mgcfinp,'visible','on');
    set(handles.funccalc,'visible','on');
end 
end

function newton_raphson_Callback(hObject, eventdata, handles)
if get(handles.newton_raphson,'value')==1
    set(handles.false,'value',0);
    set(handles.fixed,'value',0);
    set(handles.secant,'value',0);
    set(handles.bisection,'value',0);
    set(handles.xu,'visible','off');
    set(handles.xuu,'visible','off');
    set(handles.xl,'visible','off');
    set(handles.xll,'visible','off');
    set(handles.xi,'visible','on');
    set(handles.xii,'visible','on');
    set(handles.X_initial_1,'visible','off');
    set(handles.xinit_1,'visible','off');
    set(handles.maxtext,'visible','on');
    set(handles.errtext,'visible','on');
    set(handles.maxi,'visible','on');
    set(handles.es,'visible','on');
    set(handles.magictext,'visible','off');
    set(handles.mgcminp,'visible','off');
    set(handles.mgcfinp,'visible','off');
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.gx,'visible','off');
    set(handles.funccalc,'visible','off');
end 
end

function secant_Callback(hObject, eventdata, handles)
if get(handles.secant,'value')==1
    set(handles.newton_raphson,'value',0);
    set(handles.fixed,'value',0);
    set(handles.false,'value',0);
    set(handles.bisection,'value',0);
    set(handles.xu,'visible','off');
    set(handles.xuu,'visible','off');
    set(handles.xl,'visible','off');
    set(handles.xll,'visible','off');
    set(handles.xi,'visible','on');
    set(handles.xii,'visible','on');
    set(handles.X_initial_1,'visible','on');
    set(handles.xinit_1,'visible','on');
    set(handles.maxtext,'visible','on');
    set(handles.errtext,'visible','on');
    set(handles.maxi,'visible','on');
    set(handles.es,'visible','on');
    set(handles.magictext,'visible','off');
    set(handles.mgcminp,'visible','off');
    set(handles.mgcfinp,'visible','off');
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.gx,'visible','off');
    set(handles.funccalc,'visible','off');
end 
end

function true_value_Callback(hObject, eventdata, handles)
end

function true_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function edit15_Callback(hObject, eventdata, handles)
end

function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function xu_Callback(hObject, eventdata, handles)
end

function xu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function xi_Callback(hObject, eventdata, handles)
end

function xi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function runbut_Callback(hObject, eventdata, handles)
set(handles.hello,'visible','on');
syms x;
if isempty(get(handles.maxi,'string'))
    imax=50;
    set(handles.maxi,'string','50');
else
    imax=str2num(get(handles.maxi,'string'));
end
if isempty(get(handles.es,'string'))
    errors=0.00001;
    set(handles.es,'string','0.00001');
else
    errors=str2num(get(handles.es,'string'));
end
es=errors;
%%%
filepath=get(handles.fpath,'String');
if get(handles.finp,'value')==1  
    if isfile(filepath)
        fileID = fopen(filepath,'r');
        s = fscanf(fileID,'%s');
    else
        warningMessage = sprintf('This file does not exists');
        uiwait(msgbox(warningMessage));
    return;
    end
end   
%%%
if get(handles.minp,'value')==1
    s=get(handles.func,'String');
end
iterations = 0;
s_fun = str2func(['@(x)' s]);
s_funv = str2func(['@(x)' vectorize(s)]);% Vectorized
true_roots = solve(s_funv,x);
x = linspace(-10, 10, 25);
axes(handles.graph);
plot(x, s_funv(x));
grid
time1 = clock;
%Bisection 
if get(handles.bisection,'value')==1
    ea=inf;
    xl=str2num(get(handles.xl,'string'));
    xu= str2num(get(handles.xu,'string'));
    es=errors;
    vector = [];
    check=(s_fun(xl))*(s_fun(xu));
    if(check>0)
        warningMessage = sprintf('No Bracket');
        uiwait(msgbox(warningMessage));
        return;
    end
    xr=xl;
    for i=1:1:imax
        iterations = iterations + 1;
        xrold=xr;
        xr=(xu+xl)/2;
        if(i>1)
            ea = abs((xr-xrold)/xr);
        end
        vector = [vector; xl xu xr s_fun(xl) s_fun(xu) s_fun(xr) ea];

        test=(s_fun(xl))*(s_fun(xr));
        if(test>0)
            xl=xr;
        elseif(test<0)
            xu=xr;
        else
            ea=0;
        end
        if (ea < es)
            break;
        end
    end
  axes(handles.graph);
      es=errors;

  hold on
  plot(xr, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
  grid on;
  set(handles.hello,'ColumnName',{'xl'; 'xu'; 'xr'; 'f(xl)'; 'f(xu)';'f(xr)'; 'Ea'},'Data',vector);
  approximate_root = xr;
end 

%False Position
if get(handles.false,'value')==1
    ea=inf;
    xl=str2num(get(handles.xl,'string'));
    xu= str2num(get(handles.xu,'string'));
    vector = [];
    check=(s_fun(xl))*(s_fun(xu));
    if(check>0)
        warningMessage = sprintf('No Bracket');
        uiwait(msgbox(warningMessage));
        return;
     end
    iterations = 0;
    for i=1:1:imax
        iterations = iterations + 1;
        fl=s_fun(xl);
        fu=s_fun(xu);
        xr = ((xl * fu) - (xu * fl)) / (fu - fl);
        if(i>1)
            ea = abs((xr-xrold)/xr);
        end
        xrold=xr;
        vector = [vector; xl xu xr s_fun(xl) s_fun(xu) s_fun(xr) ea];
        test=(s_fun(xl))*(s_fun(xr));
        if(test>0)
            xl=xr;
        elseif(test<0)
            xu=xr;
        else
            ea=0;
        end
        if (ea < es)
            break;
        end
    end
  axes(handles.graph);
  hold on
  plot(xr, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
  grid on;
  set(handles.hello,'ColumnName',{'xl'; 'xu'; 'xr'; 'f(xl)'; 'f(xu)';'f(xr)'; 'Ea'},'Data',vector);
  approximate_root = xr;
end 

%Fixed Point
if get(handles.fixed,'value')==1
    ea=inf;
    iterations=0;
    xi=str2num(get(handles.xi,'string'));
    vector=[];
    if get(handles.mgcfinp,'value')==1
    filepath2=get(handles.mgcpath,'String');
    if isfile(filepath2)
        fileID = fopen(filepath2,'r');
        g = fscanf(fileID,'%s');
    else
        warningMessage = sprintf('This file does not exists');
        uiwait(msgbox(warningMessage));
        return;
    end
    elseif get(handles.mgcminp,'value')==1
    g=get(handles.mgcfunc,'String');
    elseif get(handles.funccalc,'value')==1
        syms x;
        g_dash = diff(s_fun, x);
        g_dash=matlabFunction(g_dash);
        if((abs(g_dash(xi))+1)<1)
             warningMessage = sprintf('Convergence is guaranteed');
             uiwait(msgbox(warningMessage));
        else
             warningMessage = sprintf('Convergence is not guaranteed');
             uiwait(msgbox(warningMessage));
        end
        for i = 1:imax
            iterations=iterations+1;
            xnew = s_fun(xi)+xi;
                ea = abs((xnew - xi) / xnew);
            vector = [vector; xi xnew ea];
            xi = xnew;
            test = s_fun(xnew)+xnew;
            
            if (i > 1)
                if (test == 0)
                    ea = 0;
                end
            
                if (ea < es)
                    break
                end
            end
        end
        set(handles.hello,'ColumnName',{'xi'; 'xi+1'; 'Ea'},'Data',vector);
        approximate_root = xnew;
        es=errors;
        axes(handles.graph);
        hold on
        grid on;
        plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
        time2 = clock;
differ = time2-time1;
elapsed_time=differ(6);
roots = transpose(true_roots);
n = length(roots);
closestRoot=roots(1);
minDiff = abs(roots(1)-approximate_root);
for k=2:n
    if abs(roots(k)-approximate_root) < minDiff 
        closestRoot = roots(k);
        minDiff =  abs(roots(k)-approximate_root);
    end
end
set(handles.precision, 'string', ea);
set(handles.Approximate_root, 'string', approximate_root);
set(handles.true_value, 'string', num2str(double(closestRoot)));
Et = abs((closestRoot - approximate_root) / closestRoot);
set(handles.true_error, 'string', num2str(double(Et)));
set(handles.number_of_iterations, 'string', iterations);
set(handles.execution_time, 'string', elapsed_time);
        return;
    end
    syms x;
    g_fun = str2func(['@(x)' g]);
    g_dash = diff(g_fun, x);
    g_dash=matlabFunction(g_dash);
    if(abs(g_dash(xi))<1)
        warningMessage = sprintf('Convergence is guaranteed');
        uiwait(msgbox(warningMessage));
    else
        warningMessage = sprintf('Convergence is not guaranteed');
        uiwait(msgbox(warningMessage));
    end
    for i=1:1:imax
        iterations=iterations+1;
        xnew=g_fun(xi);
        ea = abs((xnew-xi)/xnew);
        vector = [vector; xi xnew ea];
        test=s_fun(xnew);
        if(test==0)
            ea=0;
        end
        if (ea < es)
            break;
        end
        xi=xnew;
    end
    set(handles.hello,'ColumnName',{'xi'; 'xi+1';'Ea'},'Data',vector);
    approximate_root = xnew;
    axes(handles.graph);
    hold on
    grid on;
    plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
    
end 

%Newton Raphson
if get(handles.newton_raphson,'value')==1
xold=str2double(get(handles.xi,'string'));
syms x;
B =  diff(s_fun, x);
if B==0
     disp('Derivative equals 0');
     return;
end
C=matlabFunction(B);
vector = [];

for i = 1:imax
    fx = s_fun(xold);
    D=diff(B, x);
    if D==0
        fdash = double(B);
    else
        fdash = C(xold);
    end
    
    xnew = xold - (fx / fdash);
    iterations=iterations+1;
    Ea=inf;
    if (i > 1)
        Ea = abs((xnew - xold) / xnew);
    end
    vector = [vector; xold xnew Ea];
    xold = xnew;
    test = s_fun(xnew);
   
    if (i > 1)
        if (test == 0)
            Ea = 0;
        end
    
        if (Ea < es)
            break
        end
    end
end
ea=Ea;
set(handles.hello,'ColumnName',{'xi'; 'xi+1'; 'Ea'},'Data',vector);
approximate_root = xnew;
es=errors;
axes(handles.graph);
hold on
grid on;
plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
end

%Secant Method
if get(handles.secant,'value')==1
    xold=str2double(get(handles.xinit_1,'string'));
    xi=str2double(get(handles.xi,'string'));
    vector=[];
    for i=1:1:imax
        fxi=s_fun(xi);
        fxold=s_fun(xold);
        xnew=xi-(fxi*(xold-xi))/(fxold-fxi);
        ea = abs((xnew-xi)/xnew);
        vector = [vector; xold xi xnew ea];
        xold=xi;
        xi=xnew; 
        test=s_fun(xnew);
        if(test==0)
            ea=0;
            break;
        end
        if (ea < es)
            break;
        end
    end
    set(handles.hello,'ColumnName',{'xi-1'; 'xi';'xi+1'; 'Ea'},'Data',vector);
    approximate_root = xnew;

    es=errors;
    axes(handles.graph);
    hold on
    grid on;
    plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
end 
time2 = clock;
differ = time2-time1;
elapsed_time=differ(6);
roots = transpose(true_roots);
n = length(roots);
f=0;
if(~any(imag(roots(1))))
closestRoot=roots(1);
minDiff = abs(roots(1)-approximate_root);
f=1;
end
for k=2:n
    if(~any(imag(roots(k))))
         if abs(roots(k)-approximate_root) < minDiff 
             closestRoot = roots(k);
             minDiff =  abs(roots(k)-approximate_root);
             f=1;
         end
    end
end
if(f==0)
     warningMessage = sprintf('Function has no real roots');
     handles.hello.Data = {};
     uiwait(msgbox(warningMessage));
     return;
end    
set(handles.precision, 'string', ea);
set(handles.Approximate_root, 'string', approximate_root);
set(handles.true_value, 'string', num2str(double(closestRoot)));
Et = abs((closestRoot - approximate_root) / closestRoot);
set(handles.true_error, 'string', num2str(double(Et)));

set(handles.number_of_iterations, 'string', iterations);
set(handles.execution_time, 'string', elapsed_time);
if isnan (ea)
     warningMessage = sprintf('Divergence occurs');
     uiwait(msgbox(warningMessage));
end        
end

function func_Callback(hObject, eventdata, handles)
end

function func_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function xl_Callback(hObject, eventdata, handles)
end

function xl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function xinit_1_Callback(hObject, eventdata, handles)
end

function xinit_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function maxi_Callback(hObject, eventdata, handles)
end

function maxi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function es_Callback(hObject, eventdata, handles)
end

function es_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function fpath_Callback(hObject, eventdata, handles)
end

function fpath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function minp_Callback(hObject, eventdata, handles)
    if get(handles.minp,'value')==1
    set(handles.finp,'value',0);
    set(handles.mtext1,'visible','on');
    set(handles.mtext2,'visible','on');
    set(handles.func,'visible','on');
    set(handles.ftext,'visible','off');
    set(handles.fpath,'visible','off');
    end 
end

function finp_Callback(hObject, eventdata, handles)
    if get(handles.finp,'value')==1
    set(handles.minp,'value',0);
    set(handles.mtext1,'visible','off');
    set(handles.mtext2,'visible','off');
    set(handles.func,'visible','off');
    set(handles.ftext,'visible','on');
    set(handles.fpath,'visible','on');
    end 
end

function resetbut_Callback(hObject, eventdata, handles)
global stepiterations;
    set(handles.stepbut,'Enable','on');
stepiterations=1;
set(handles.minp,'value',0);
set(handles.finp,'value',0);
set(handles.newton_raphson,'value',0);
set(handles.fixed,'value',0);
set(handles.false,'value',0);
set(handles.bisection,'value',0);
set(handles.func,'String','');
set(handles.fpath,'String','');
set(handles.maxi,'String','');
set(handles.es,'String','');
set(handles.xi,'String','');
set(handles.xl,'String','');
set(handles.xinit_1,'String','');
set(handles.xu,'String','');
set(handles.xi,'visible','off');
set(handles.xii,'visible','off');
set(handles.X_initial_1,'visible','off');
set(handles.xinit_1,'visible','off');
set(handles.xu,'visible','off');
set(handles.xuu,'visible','off');
set(handles.xl,'visible','off');
set(handles.xll,'visible','off');
set(handles.mtext1,'visible','off');
set(handles.mtext2,'visible','off');
set(handles.func,'visible','off');
set(handles.ftext,'visible','off');
set(handles.fpath,'visible','off');
set(handles.maxtext,'visible','off');
set(handles.errtext,'visible','off');
set(handles.maxi,'visible','off');
set(handles.es,'visible','off');
set(handles.magictext,'visible','off');
set(handles.mgcminp,'visible','off');
set(handles.mgcfinp,'visible','off');
set(handles.mgcminp,'value',0);
set(handles.mgcfinp,'value',0);
set(handles.mgcfpathtext,'visible','off');
set(handles.mgcminptext,'visible','off');
set(handles.mgcfunc,'visible','off');
set(handles.mgcfunc,'String','');
set(handles.mgcpath,'visible','off');
set(handles.mgcpath,'String','');
set(handles.gx,'visible','off');
set(handles.funccalc,'visible','off');
set(handles.funccalc,'value',0);
set(handles.Approximate_root, 'string', '');
set(handles.precision, 'string', '');
set(handles.true_value, 'string', '');
set(handles.number_of_iterations, 'string', '');
set(handles.true_error, 'string', '');
set(handles.execution_time, 'string', '');
set(handles.hello, 'visible', 'off');
handles.hello.Data = {};
cla(handles.graph);
end

function stepbut_Callback(hObject, eventdata, handles)
cla(handles.graph);
set(handles.hello,'visible','on');
syms x;
global stepiterations;
iterations=0;
if isempty(get(handles.maxi,'string'))
    imax=50;
    set(handles.maxi,'string','50');
else
    imax=str2num(get(handles.maxi,'string'));
end
if isempty(get(handles.es,'string'))
    errors=0.00001;
    set(handles.es,'string','0.00001');
else
    errors=str2num(get(handles.es,'string'));
end
filepath=get(handles.fpath,'String');
if get(handles.finp,'value')==1  
    if isfile(filepath)
        fileID = fopen(filepath,'r');
        s = fscanf(fileID,'%s');
    else
        warningMessage = sprintf('This file does not exists');
        uiwait(msgbox(warningMessage));
    return;
    end
end   
if get(handles.minp,'value')==1
    s=get(handles.func,'String');
end
s_fun = str2func(['@(x)' s]);
s_funv = str2func(['@(x)' vectorize(s)]);% Vectorized
true_roots = solve(s_funv,x);
x = linspace(-10, 10, 25);
axes(handles.graph);
plot(x, s_funv(x));
grid on;
es=errors;
time1=clock;
if get(handles.bisection,'value')==1
    ea=inf;
    xl=str2num(get(handles.xl,'string'));
    xu= str2num(get(handles.xu,'string'));
    es=errors;
    vector = [];
    check=(s_fun(xl))*(s_fun(xu));
    if(check>0)
        disp('No Bracket');
        return;
    end
    xr=xl;
    limit=min(stepiterations,imax);
    for i=1:1:limit
        iterations = iterations + 1;
        xrold=xr;
        xr=(xu+xl)/2;
        if(i>1)
            ea = abs((xr-xrold)/xr);
        end
                vector = [vector; xl xu xr s_fun(xl) s_fun(xu) s_fun(xr) ea];

        test=(s_fun(xl))*(s_fun(xr));
        if(test>0)
            xl=xr;
        elseif(test<0)
            xu=xr;
        else
            ea=0;
        end
        if (ea < es)
             set(handles.stepbut,'Enable','off');
            break;
        end
    end
  axes(handles.graph);
  hold on
  plot(xr, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
  approximate_root = xr;
  grid on
  set(handles.hello,'ColumnName',{'xl'; 'xu'; 'xr'; 'f(xl)'; 'f(xu)';'f(xr)'; 'Ea'},'Data',vector);

end 

%False Position
if get(handles.false,'value')==1
    ea=inf;
    xl=str2num(get(handles.xl,'string'));
    xu= str2num(get(handles.xu,'string'));
    es=errors;
    vector = [];
    check=(s_fun(xl))*(s_fun(xu));
    if(check>0)
        disp('No Bracket');
            set(handles.stepbut,'Enable','off');
    end
    limit=min(stepiterations,imax);
    for i=1:1:limit
        iterations = iterations + 1;
        fl=s_fun(xl);
        fu=s_fun(xu);
        xr = ((xl * fu) - (xu * fl)) / (fu - fl);
        if(i>1)
            ea = abs((xr-xrold)/xr);
        end
        xrold=xr;
        vector = [vector; xl xu xr s_fun(xl) s_fun(xu) s_fun(xr) ea];
        test=(s_fun(xl))*(s_fun(xr));
        if(test>0)
            xl=xr;
        elseif(test<0)
            xu=xr;
        else
            ea=0;
                set(handles.stepbut,'Enable','off');
        end
        if (ea < es)
                set(handles.stepbut,'Enable','off');
            break;
        end
    end
  axes(handles.graph);
  approximate_root = xr;

  hold on
  plot(xr, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
  grid on
  set(handles.hello,'ColumnName',{'xl'; 'xu'; 'xr'; 'f(xl)'; 'f(xu)';'f(xr)'; 'Ea'},'Data',vector);
end 
%Fixed Point
if get(handles.fixed,'value')==1
    ea=inf;
    iterations=0;
    xi=str2num(get(handles.xi,'string'));
    vector=[];
    if get(handles.mgcfinp,'value')==1
    filepath2=get(handles.mgcpath,'String');
     if isfile(filepath2)
        fileID = fopen(filepath2,'r');
        g = fscanf(fileID,'%s');
    else
        warningMessage = sprintf('This file does not exists');
        uiwait(msgbox(warningMessage));
        return;
    end
    elseif get(handles.mgcminp,'value')==1
    g=get(handles.mgcfunc,'String');
    elseif get(handles.funccalc,'value')==1
     limit=min(stepiterations,imax);
        for i=1:1:limit
                iterations=iterations+1;
            xnew = s_fun(xi)+xi;
                ea = abs((xnew - xi) / xnew);
 
            vector = [vector; xi xnew ea];
            xi = xnew;
            test = s_fun(xnew)+xnew;
            
            if (i > 1)
                if (test == 0)
                    ea = 0;
                       set(handles.stepbut,'Enable','off');
                end
            
                if (ea < es)
                       set(handles.stepbut,'Enable','off');
                    break
                end
            end
        end
set(handles.hello,'ColumnName',{'xi'; 'xi+1'; 'Ea'},'Data',vector);
%set(handles.Approximate_root, 'string', xnew);
approximate_root = xnew;
es=errors;
axes(handles.graph);
hold on
plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
stepiterations=stepiterations+1;
time2 = clock;
differ = time2-time1;
elapsed_time=differ(6)
roots = transpose(true_roots);
n = length(roots);
closestRoot=roots(1);
minDiff = abs(roots(1)-approximate_root);
for k=2:n
    if abs(roots(k)-approximate_root) < minDiff 
        closestRoot = roots(k);
        minDiff =  abs(roots(k)-approximate_root);
    end
end
set(handles.Approximate_root, 'string', approximate_root);
set(handles.precision, 'string', ea);
set(handles.true_value, 'string', num2str(double(closestRoot)));
set(handles.number_of_iterations, 'string', iterations);
Et = abs((closestRoot - approximate_root) / closestRoot);
set(handles.true_error, 'string', num2str(double(Et)));
set(handles.execution_time, 'string', elapsed_time);
        return;
    end
    g_fun = str2func(['@(x)' g]);
       limit=min(stepiterations,imax);
        for i=1:1:limit
        iterations=iterations+1;
        xnew=g_fun(xi);
        ea = abs((xnew-xi)/xnew);
        vector = [vector; xi xnew ea];
        test=s_fun(xnew);
        if(test==0)
            set(handles.stepbut,'Enable','off');
            ea=0;
        end
        if (ea < es)
            set(handles.stepbut,'Enable','off');
            break;
        end
        xi=xnew;
    end
    set(handles.hello,'ColumnName',{'xi'; 'xi+1';'Ea'},'Data',vector);
    approximate_root = xnew;
    axes(handles.graph);
    hold on
    plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
end 


%Newton Raphson
if get(handles.newton_raphson,'value')==1
xold=str2double(get(handles.xi,'string'))
syms x;
B =  diff(s_fun, x);
if B==0
     disp('Derivative equals 0');
     return;
end
C=matlabFunction(B);
vector = [];
limit=min(stepiterations,imax);
    for i=1:1:limit
    fx = s_fun(xold);
    D=diff(B, x);
    if D==0
        fdash = double(B);
    else
        fdash = C(xold);
    end
    
    xnew = xold - (fx / fdash);
    iterations=iterations+1;
    Ea=inf;
    if (i > 1)
        Ea = abs((xnew - xold) / xnew);
    end
    vector = [vector; xold xnew Ea];
    xold = xnew;
    %xr=xnew;
    test = s_fun(xnew);
   
    if (i > 1)
        if (test == 0)
            Ea = 0;
                set(handles.stepbut,'Enable','off');
        end
    
        if (Ea < es)
                set(handles.stepbut,'Enable','off');
            break
        end
    end
end
set(handles.hello,'ColumnName',{'xi'; 'xi+1'; 'Ea'},'Data',vector);
approximate_root = xnew;
es=errors;
ea=Ea;
axes(handles.graph);
hold on
plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
end
 
%Secant Method
if get(handles.secant,'value')==1
    xold=str2double(get(handles.xinit_1,'string'));
    xi=str2double(get(handles.xi,'string'));
    vector=[];
    limit=min(stepiterations,imax);
    for i=1:1:limit
        fxi=s_fun(xi);
        fxold=s_fun(xold);
        xnew=xi-(fxi*(xold-xi))/(fxold-fxi);
        ea = abs((xnew-xi)/xnew);
        vector = [vector; xold xi xnew ea];
        xold=xi;
        xi=xnew; 
        test=s_fun(xnew);
        if(test==0)
            ea=0;
            set(handles.stepbut,'Enable','off');
            break;
        end
        if (ea < es)
            set(handles.stepbut,'Enable','off');
            break;
        end
    end
    set(handles.hello,'ColumnName',{'xi-1'; 'xi';'xi+1'; 'Ea'},'Data',vector);
    
    approximate_root = xnew;
    axes(handles.graph);
    hold on
    plot(xnew, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25);
end 
stepiterations=stepiterations+1;
time2 = clock;
differ = time2-time1;
elapsed_time=differ(6)
roots = transpose(true_roots);
n = length(roots);
closestRoot=roots(1);
minDiff = abs(roots(1)-approximate_root);
for k=2:n
    if abs(roots(k)-approximate_root) < minDiff 
        closestRoot = roots(k);
        minDiff =  abs(roots(k)-approximate_root);
    end
end
set(handles.Approximate_root, 'string', approximate_root);
set(handles.precision, 'string', ea);
set(handles.true_value, 'string', num2str(double(closestRoot)));
set(handles.number_of_iterations, 'string', iterations);
Et = abs((closestRoot - approximate_root) / closestRoot);
set(handles.true_error, 'string', num2str(double(Et)));
set(handles.execution_time, 'string', elapsed_time);
if isnan (ea)
     warningMessage = sprintf('Divergence occurs');
     uiwait(msgbox(warningMessage));
end        
end

function true_error_Callback(hObject, eventdata, handles)
end

function true_error_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function mgcminp_Callback(hObject, eventdata, handles)
if get(handles.mgcminp,'value')==1
    set(handles.mgcfinp,'value',0);
    set(handles.funccalc,'value',0);
    set(handles.mgcminptext,'visible','on');
    set(handles.mgcfunc,'visible','on');
    set(handles.gx,'visible','on');
end 
end

function mgcfinp_Callback(hObject, eventdata, handles)
if get(handles.mgcfinp,'value')==1
    set(handles.mgcminp,'value',0);
    set(handles.funccalc,'value',0);
    set(handles.mgcfpathtext,'visible','on');
    set(handles.mgcpath,'visible','on');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.gx,'visible','off');
end 
end

function funccalc_Callback(hObject, eventdata, handles)
if get(handles.funccalc,'value')==1
    set(handles.mgcminp,'value',0);
    set(handles.mgcfinp,'value',0);
    set(handles.mgcfpathtext,'visible','off');
    set(handles.mgcpath,'visible','off');
    set(handles.mgcminptext,'visible','off');
    set(handles.mgcfunc,'visible','off');
    set(handles.gx,'visible','off');
end 
end

function mgcfunc_Callback(hObject, eventdata, handles)
end

function mgcfunc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function mgcpath_Callback(hObject, eventdata, handles)
end

function mgcpath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit4_Callback(hObject, eventdata, handles)
end

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
