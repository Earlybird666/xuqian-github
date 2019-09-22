% script mps_array_model_int
% This script solves for the normalized velocity wave field of 
% an array of 1-D elements radiating waves through a fluid/solid
% interface using the MATLAB function ps_3Dint. Both time delay 
% and apodization laws can be specified for the array to steer 
% it and focus it in the solid.


clear
% ------------------input parameters ----------------
tic          % time the calculations
lx = 0.15;   % element length in x-direction (mm)
ly = 0.15;   % element length in y-direction (mm)
gx=0.05;     % gap length in x-direction (mm)
gy = 0.05;   % gap length in y-direction (mm)
f= 5;        % frequency (MHz)
d1=1.0;      % density, medium one (arbitrary units)
cp1 = 1480;  % compressional wave speed, medium one (m/sec)
d2=7.9;      % density, medium two  (same arbitrary units)
cp2 =5900;   % compressional wave speed, medium two (m/sec)
cs2=3200;    % shear wave speed, medium two (m/sec)
type='p';    % wave type, medium two
mat=[d1,cp1,d2,cp2,cs2,type];   % form material vector

L1 =11;        % number of elements in x-direction
L2 =11;        % number of elements in y-direction
angt =10.217;  % angle of the array (deg)
Dt0=50.8;      % height of array center from interface (mm)

theta2 =20;   % steering angle in theta direction (deg)
phi =0;      %steering angle in phi direction (deg)
DF = inf;    % focal distance (mm)
%weighting choices are 'rect','cos', 'Han', 'Ham', 'Blk', 'tri' 
ampx_type ='rect';   % weighting coeffcients in x-direction
ampy_type ='rect';   % weighting coefficients in y-direction

% field points to evaluate
xs= linspace(-5,20, 100);
zs= linspace(1, 20, 100);
y=0;
[x,z]=meshgrid(xs,zs);

% ------------- end input parameters -----------------------

c1=cp1;
if strcmp(type,'p')
    c2=cp2;
elseif strcmp(type, 's')
    c2=cs2;
else
    error('type incorrect')
end



% calculate array pitches
sx = lx+gx;
sy = ly+gy;

% compute centroid locations for the elements
Nx = 1:L1;
Ny = 1:L2;
ex =(2*Nx -1-L1)*(sx/2);
ey =(2*Ny -1 -L2)*(sy/2);

% generate time delays, put in exponential 
% and calculate amplitude weights
td =delay_laws3Dint(L1,L2,sx,sy,angt,phi,theta2,Dt0,DF,c1,c2,'n');
delay = exp(1i.*2.*pi.*f.*td);
Cx = discrete_windows(L1,ampx_type);
Cy = discrete_windows(L2,ampy_type);


% calculate normalized velocities
vx=0;
vy=0;
vz=0;
for nn=1:L1
    for ll=1:L2
        [vxe,vye,vze]= ps_3Dint(lx,ly,f,mat,ex(nn),ey(ll),angt,Dt0,x,y,z,1,1);
        vx = vx + Cx(nn)*Cy(ll)*delay(nn,ll)*vxe;
        vy = vy + Cx(nn)*Cy(ll)*delay(nn,ll)*vye;
        vz = vz + Cx(nn)*Cy(ll)*delay(nn,ll)*vze;
    end
end
  
% ------------------ outputs --------------------------
%plot results
vmag=sqrt(abs(vx).^2 +abs(vy).^2 +abs(vz).^2);
imagesc(xs,zs,vmag)       
  toc    % end of time calculations  