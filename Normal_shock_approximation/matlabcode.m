%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables having 'd' at the beginning represents the dimentional quantities. 
% And those without it represents non-dimentional quantities. 


clear;
clc;


%M1 = 10.0;
%n  = 100;                    % no of grid points
%alpha = 0.14;    
%CFL   = 0.2;

M1 = str2num(argv(){1});
n  = str2num(argv(){2});                    % no of grid points
alpha = str2num(argv(){3});    
CFL   = str2num(argv(){4});

%% VARIABLE INITIALIZATION
U = zeros(3,n);
F = zeros(3,n);
u   = zeros(1,n);
rho = zeros(1,n);
p   = zeros(1,n);
T   = zeros(1,n);
rhoe= zeros(1,n);
x   = zeros(1,n);

%% CONSTANTS
gamma = 1.4;
Rsp   = 287;                            % sp. gas constant R
Rspnd = 1/gamma;                        % non-dimensional sp. gas constant


dp1   = 1.01325e5;                      % Dimensional pressure
drho1 = 1.225;                          % Dimensional density
dT1   = 288; a1  = sqrt(gamma*Rsp*dT1); % Dimensional temperature and freestream speed of sound 

u1   = M1;                              % velocity non-dimensionalized by a1
rho1 = drho1/drho1;                     % density  non-dimensionalized by its freestream value 1.225
T1   = dT1/dT1;                         % temperature non-dimensionalized by its freestream value 288
p1   = dp1/(drho1*a1*a1);               % pressure non-dimensionalized by rho*a^2

%% RANKINE-HUGONIOT CONDITIONS
M2 = sqrt((M1^2*(gamma-1)+2)/(2*gamma*M1^2-(gamma-1)));

p2 = p1*((1+gamma*M1^2)/(1+gamma*M2^2));
T2 = T1*((1+(gamma-1)/2*M1^2)*((2*gamma/(gamma-1)*M1^2-1)))/(M1^2*(2*gamma/(gamma-1)+(gamma-1)/2));
dT2 = T2*dT1; a2  = sqrt(gamma*Rsp*dT2);
u2 = M2*(a2/a1);
rho2 = u1/u2;


%% DOMAIN DETAILS
L     = 35;                 % Domain = [-5,30] shock being at x = 0
dx    = 35/(n-1);
dt    = CFL*dx/(abs(u1+1)); % Computing the time-step based on largest eigen value
nt    = 10000;              % Number of iterations

%% INITIAL CONDITIONS
for i = 1:n    
     x(i)= -5 +(i-1)*dx;     
     if(x(i)<0)           % free-stream conditions before shock   
        rho(i) = rho1;
        u(i)   = u1;
        p(i)   = p1;
        T(i)   = T1;
     else                 % conditions after shock   
        rho(i) = rho2;
        u(i)   = u2;
        p(i)   = p2;
        T(i)   = T2;
     end     
        rhoe(i)   = p(i)/(gamma-1)+0.5*rho(i)*u(i)^2 ;  
end

step = 0;
t=0;
 

%% SOLVING
while(t<nt*dt)    
    
    % Defining Matrices
    U = [rho; rho.*u; rhoe]; 
    F = [rho.*u; rho.*u.^2+p; u.*(rhoe+p)];  
    
    % Solving and updating the matrix variables
    f(1:3,1:n-1) =  0.5*(F(1:3,1:n-1)+F(1:3,2:n)) - 0.5*alpha*(u1+1)*(U(1:3,2:n)-U(1:3,1:n-1)); 
    U(1:3,2:n-1) =  U(1:3,2:n-1)  - dt/dx*( (f(1:3,2:n-1)-f(1:3,1:n-2)));                    
         
    
    % Variables update
    rho=U(1,:);          
    u=U(2,:)./rho;
    rhoe=U(3,:);           
    p=(gamma-1)*(rhoe-0.5*rho.*u.^2 );      
    T  = p./(rho*Rspnd); 
      
    % Time-step update
    t = t + dt;         
    step = step + 1;        
    fprintf('Count %d \n', step);  % Printing iteration count to the screen
       
end


%% PLOTTING
%%% For in-built plotting
%figure(1);  
%subplot(221); plot(x,rho,'-','DisplayName',num2str(alpha),'LineWidth',1.5); xlabel('x'); ylabel('\rho'); xlim('auto'); ylim('auto'); title('Density profile');     legend('-DynamicLegend'); hold all;
%subplot(222); plot(x,u  ,'-','DisplayName',num2str(alpha),'LineWidth',1.5); xlabel('x'); ylabel('u');    xlim('auto'); ylim('auto'); title('Velocity profile');    legend('-DynamicLegend'); hold all;
%subplot(223); plot(x,p  ,'-','DisplayName',num2str(alpha),'LineWidth',1.5); xlabel('x'); ylabel('p');    xlim('auto'); ylim('auto'); title('Pressure profile');    legend('-DynamicLegend'); hold all;
%subplot(224); plot(x,T  ,'-','DisplayName',num2str(alpha),'LineWidth',1.5); xlabel('x'); ylabel('T');    xlim('auto'); ylim('auto'); title('Temperature profile'); legend('-DynamicLegend'); hold all;

%%% For plotting in tecplot
Data   = [x; rho; u; p; T];
fileID = fopen(['dat_files/data_',num2str(M1),"_",num2str(n),"_",num2str(alpha),"_",num2str(CFL),"_",'.dat'],'w');
fprintf(fileID,'ZONE alpha = %.1f" \n',alpha);
fprintf(fileID,'%6.2f  %8.4f  %8.4f  %8.4f  %8.4f   \n',Data);
fclose(fileID);













