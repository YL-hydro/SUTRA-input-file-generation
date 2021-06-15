% WRITE INPUT FILE
% NOTE THAT TO RUN THIS FILE, slsetpath.m from the SUTRALAB_INPwrite should
% be run first

clear all
close all
c=ConstantObj();  % create Constant Object for all the constants, we suggest to use c.g rather than 9.81 in the script to enhance readerbility;

%% create SUTRA.fil

fil=filObj('read_from_file','no');

% the default value for 'fid' is -1, meaning, this will not write to SUTRA.fil if export
fil.terms.('inp').('fname')  = 'SAND.inp' ; fil.terms.('inp').('fid')   = 50;
fil.terms.('ics').('fname')  = 'SAND.ics' ; fil.terms.('ics').('fid')   = 55;
fil.terms.('lst').('fname')  = 'SAND.lst' ; fil.terms.('lst').('fid')   = 60;
fil.terms.('rst').('fname')  = 'SAND.rst' ; fil.terms.('rst').('fid')   = 66;
fil.terms.('nod').('fname')  = 'SAND.nod' ; fil.terms.('nod').('fid')   = 30;
fil.terms.('ele').('fname')  = 'SAND.ele' ; fil.terms.('ele').('fid')   = 40;
fil.terms.('obs').('fname')  = 'SAND.obs' ; fil.terms.('obs').('fid')   = 70;
fil.terms.('smy').('fname')  = 'SAND.smy' ; fil.terms.('smy').('fid')   = 80;
fil.terms.('bcof').('fname') = 'SAND.bcof'; fil.terms.('bcof').('fid')  = 91;
fil.terms.('bcos').('fname') = 'SAND.bcos'; fil.terms.('bcos').('fid')  = 93;
fil.terms.('bcop').('fname') = 'SAND.bcop'; fil.terms.('bcop').('fid')  = 92;
fil.terms.('bcou').('fname') = 'SAND.bcou'; fil.terms.('bcou').('fid')  = 94;


fil.export_to_file();

%% Constants and values
c_saltwater_kgPkg          = 0.035;
% c_freshwater_kgPkg         = 0.001;
initial_head_aquifer_m     = 0;
% permeability_silt_m2       = 5.38e-15;
permeability_sand_m2       = 8.00e-12;



%% CALCULATE THE NODES AND GENERATE NODES X,Y INFO
% control point input for block one along y axis
x1 = 0;
x2 = 250;
nex = 500; %Number of segments along x

y(1,1)=-5;
y(1,2)=-5;

ney_section(1)=30;

y(2,1)=-2;
y(2,2)=-2;

ney_section(2)=60;

y(3,1)=1;
y(3,2)=0.4;

ney_section(3)=50;

y(4,1)=1.5;
y(4,2)=0.9;

% generate coordinates of nodes in BLOCK ONE ALONG X for dateset14
ney     =sum(ney_section);
nx      =nex+1;
ny      =ney+1;

x_array = x1:(x2-x1)/nex:x2;
y_keyline = size (y,1); %control lines in height

for i=1:y_keyline

	y_key_interval=(y(i,2)-y(i,1))/nex;

	if y_key_interval==0
		y_array(i,1:nx) = y(i,1);
	else
		y_array(i,1:nx) = y(i,1):y_key_interval:y(i,2);
	end

end

for i=1:nx

        x_nod_array ( (i-1)*ny+1: i*ny )= x_array(i);

	for j=1:y_keyline-1
        
        y_interval (j,i) = (y_array(j+1,i)-y_array(j,i))./ney_section(j);      
        y_location (1)   = 0;
        y_location (j+1) = sum(ney_section(1:j));
                
        y_nod_array ( ney*(i-1)+i + y_location(j):ney*(i-1)+i + y_location (j+1) )= y_array(j,i): y_interval(j,i) : y_array(j+1,i) ; 
    end
    
	
end
% BLOCK TWO ALONG X
x1 = 250;
x2 = 260;
nex = 40; %Number of segments along x


y(1,1)=-5;
y(1,2)=-5;

ney_section(1)=30; % NEED TO MAKE SURE THIS IS THE SAME ALONG ALL BLOCKS

y(2,1)=-2;
y(2,2)=-2;

ney_section(2)=60;

y(3,1)=0.4;
y(3,2)=-1.3;

ney_section(3)=50;

y(4,1)=0.9;
y(4,2)=-0.8;

dz=0.01;

%% process coordinate of nodes for dateset14
p2 = length (x_nod_array);
ney     =sum(ney_section);
nx      =nex+1;
ny      =ney+1;

x_array = x1:(x2-x1)/nex:x2;
y_keyline = size (y,1); %control lines in height

for i=1:y_keyline

	y_key_interval=(y(i,2)-y(i,1))/nex;

	if y_key_interval==0
		y_array(i,1:nx) = y(i,1);
	else
		y_array(i,1:nx) = y(i,1):y_key_interval:y(i,2);
	end

end

for i=1:nx

        x_nod_array ( (i-1)*ny+1+p2-(ney+1): i*ny+p2-(ney+1) )= x_array(i);

	for j=1:y_keyline-1
        
        y_interval (j,i) = (y_array(j+1,i)-y_array(j,i))./ney_section(j);      
        y_location (1)   = 0;
        y_location (j+1) = sum(ney_section(1:j));
                
        y_nod_array ( ney*(i-1)+i + y_location(j)+p2-(ney+1):ney*(i-1)+i + y_location (j+1)+p2-(ney+1) )= y_array(j,i): y_interval(j,i) : y_array(j+1,i) ; 
    end
    
	
end
% BLOCK THREE ALONG X
x1 = 260;
x2 = 280;
nex = 40; %Number of segments along x


y(1,1)=-5;
y(1,2)=-5;

ney_section(1)=30;

y(2,1)=-2;
y(2,2)=-2;

ney_section(2)=60;

y(3,1)=-1.3;
y(3,2)=-1.3;

ney_section(3)=50;

y(4,1)=-0.8;
y(4,2)=-0.8;

dz=0.01;

p2 = length (x_nod_array);
ney     =sum(ney_section);
nx      =nex+1;
ny      =ney+1;

x_array = x1:(x2-x1)/nex:x2;
y_keyline = size (y,1); %control lines in height

for i=1:y_keyline

	y_key_interval=(y(i,2)-y(i,1))/nex;

	if y_key_interval==0
		y_array(i,1:nx) = y(i,1);
	else
		y_array(i,1:nx) = y(i,1):y_key_interval:y(i,2);
	end

end

for i=1:nx

        x_nod_array ( (i-1)*ny+1+p2-(ney+1): i*ny+p2-(ney+1) )= x_array(i);

	for j=1:y_keyline-1
        
        y_interval (j,i) = (y_array(j+1,i)-y_array(j,i))./ney_section(j);      
        y_location (1)   = 0;
        y_location (j+1) = sum(ney_section(1:j));
                
        y_nod_array ( ney*(i-1)+i + y_location(j)+p2-(ney+1):ney*(i-1)+i + y_location (j+1)+p2-(ney+1) )= y_array(j,i): y_interval(j,i) : y_array(j+1,i) ; 
    end
    
	
end

%% CHECK THE NODES
figure
scatter (x_nod_array,y_nod_array,'.');

%%
nx = 500+40+40+1; % note here that nx needs to be changed to the total nx
nex = 580; % change the nex to 250
ny = ny; % ny keeps the same 

nn=nx*ny; % number of nodes
ne=nex*ney; % number of elements

ii = (1:nn)';

nod_idx_mtx = reshape(ii,ny,nx);% note that nod_idx_mtx requires to be flipped later
x_nod_mtx = reshape(x_nod_array,ny,nx);
y_nod_mtx = reshape(y_nod_array,ny,nx); % note that y_nod_mtx requires to be flipped later

%calculate dx, dy
keynodes            = zeros(size(x_nod_mtx,1),size(x_nod_mtx,2)+1);
keynodes(:,2:end-1) = (x_nod_mtx(:,1:end-1)+x_nod_mtx(:,2:end))/2;
keynodes(:,1)       = x_nod_mtx(:,1);
keynodes(:,end)     = x_nod_mtx(:,end);
dx_cell_mtx         = diff(keynodes,1,2); 

keynodes            = zeros(size(y_nod_mtx,1)+1,size(y_nod_mtx,2));
keynodes(2:end-1,:) = (y_nod_mtx(1:end-1,:)+y_nod_mtx(2:end,:))/2;
keynodes(1,:)       = y_nod_mtx(1,:);
keynodes(end,:)     = y_nod_mtx(end,:);
dy_cell_mtx         = diff(keynodes,1,1); 

nod_idx_mtx_g = flip(nod_idx_mtx); % flip the mtx
x_nod_mtx_g = flip(x_nod_mtx);
y_nod_mtx_g = flip(y_nod_mtx);
dx_cell_mtx_g = flip(dx_cell_mtx);
dy_cell_mtx_g = flip(dy_cell_mtx);

% to generate the ele arrary - the central point of each ele
idx = 1; % node count

iin1 = zeros(ne,1);
iin2 = zeros(ne,1);
iin3 = zeros(ne,1);
iin4 = zeros(ne,1);
x_ele_array=zeros(ne,1);
y_ele_array=zeros(ne,1);

for j  = 1:nex
    for i = 1:ney
        iin1(idx) = nod_idx_mtx(i,j);
        iin2(idx) = nod_idx_mtx(i,j+1);
        iin3(idx) = nod_idx_mtx(i+1,j+1);
        iin4(idx) = nod_idx_mtx(i+1,j);
        x_ele_array(idx)= mean([x_nod_array(iin1(idx)),x_nod_array(iin2(idx)),x_nod_array(iin3(idx)),x_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        y_ele_array(idx)= mean([y_nod_array(iin1(idx)),y_nod_array(iin2(idx)),y_nod_array(iin3(idx)),y_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        idx       = idx+1;

    end
end

x_ele_mtx_m = reshape(x_ele_array, ney,nex); % x,y coordiante mtx of the central point of each element
y_ele_mtx_m = reshape(y_ele_array, ney,nex);


x_ele_mtx_m_g = flip(x_ele_mtx_m); % flip
y_ele_mtx_m_g = flip(y_ele_mtx_m);

%% NOW START TO PRINT OUT .INP FILE AND THE CONTENT

inp = inpObj('SAND','read_from_file','no');   % setup a empty inpObj

% dataset 1
inp.title1 = 'SAND ET RAIN';
inp.title2 = 'INP by SutraLab';

% dataset 2a
%'SUTRA VERSION 2.2 SOLUTE TRANSPORT'
%'2D REGULAR MESH' 105 4001

inp.vermin = '2.2';  % note this should be a string not a number;
inp.simula = 'SOLUTE';


% dataset 2b
 %        2d mesh          ==>   ktype(1) = 2
 %        3d mesh          ==>   ktype(1) = 3
 %        irregular mesh   ==>   ktype(2) = 0
 %        layered mesh     ==>   ktype(2) = 1
 %        regular mesh     ==>   ktype(2) = 2
 %        blockwise mesh   ==>   ktype(2) = 3

inp.ktype(1)  = 2;  % 2D mesh
inp.mshtyp{1} = '2D';
inp.mshtyp{2} = 'REGULAR';

inp.nn1 = ny;
inp.nn2 = nx;


% ##  DATASET 3:  Simulation Control Numbers
inp.nn   = nn;
inp.ne   = ne;
inp.npbc = 0;   %length(ipbc);  revised after dataset 19
inp.nubc = 0;
inp.nsop = 0;       % dataset 17
inp.nsou = 0;
inp.nobs = 0;


%%##  DATASET 4:  Simulation Mode Options

inp.cunsat = 'UNSATURATED';
inp.cssflo = 'TRANSIENT FLOW';
inp.csstra = 'TRANSIENT TRANSPORT';
inp.cread  = 'COLD' ;
inp.istore = 9999;


%%##  DATASET 5:  Numerical Control Parameters
inp.up   = 0;
inp.gnup = 0.01;
inp.gnuu = 0.01;


%  DATASET 6:  Temporal Control and Solution Cycling Data
%  
inp.nsch  = 1;
inp.npcyc = 1;
inp.nucyc = 1;

%DATASET 6:  Temporal Control and Solution Cycling Data

inp.schnam = 'TIME_STEPS';
inp.schtyp = 'TIME CYCLE';
inp.creft  = 'ELAPSED';
%inp.scalt  = 6000;   %reduce to 3000
%inp.scalt  = 3000;   %reduce to 3000
inp.scalt  = 600;   %reduce to 3000
%inp.scalt  = 150;   %reduce to 3000
%inp.scalt  = 6000;   %reduce to 3000
%inp.scalt  = 0.6;   %reduce to 3000
inp.ntmax  = 2160;
inp.timei  = 0;
inp.timel  = 1.e99;
inp.timec  = 1.;
inp.ntcyc  = 9999;
inp.tcmult = 1;
inp.tcmin  = 1.e-20;
inp.tcmax  = 1.e99;


%##  DATASET 7:  ITERATION AND MATRIX SOLVER CONTROLS
%##  [ITRMAX]        [RPMAX]        [RUMAX]
inp.itrmax = 1000;
inp.rpmax  = 1e+5;
%inp.rpmax = 5e-2;  TO200317 too strigent for the first step
inp.rumax  = 1.0e-1;
% ##  [CSOLVP]  [ITRMXP]         [TOLP]
inp.csolvp = 'ORTHOMIN' ;
inp.itrmxp = 3000;
inp.tolp   = 1e-11;

%##  [CSOLVU]  [ITRMXU]         [TOLU]
inp.csolvu = 'ORTHOMIN';
inp.itrmxu = 3000;
inp.tolu   = 1e-11;


%##  DATASET 8:  Output Controls and Options
%## [NPRINT]  [CNODAL]  [CELMNT]  [CINCID]  [CPANDS]  [CVEL]  [CCORT] [CBUDG]   [CSCRN]  [CPAUSE]
%   2920        'N'        'N'        'N'        'Y'     'Y'        'Y'    'Y'      'Y' 'Y' 'Data Set 8A'

inp.nprint = 10;
inp.cnodal = 'N';
inp.celmnt = 'N';
inp.cincid = 'N';
inp.cpands = 'N';
inp.cvel   = 'N';
inp.ccort  = 'N';
inp.cbudg  = 'N';
inp.cscrn  = 'N';   % screen output
inp.cpause = 'N';


%## [NCOLPR]    [NCOL]
%     -1000  'N'  'X'  'Y'  'P'  'U'  'S'  '-' 
%## [LCOLPR]    [LCOL]
%     1000 'E'  'X'  'Y'  'VX' 'VY' '-' 
%## [NOBCYC]    [INOB]
%
%##  [NBCFPR]  [NBCSPR]  [NBCPPR]  [NBCUPR]  [CINACT]
%     1000         1000       1000      1000       'Y'
%##
%##  DATASET 9:  Fluid Properties
%##     [COMPFL]           [CW]       [SIGMAW]        [RHOW0]       [URHOW0]        [DRWDU]        [VISC0]
%         0.0                1.         3.890D-10         1.0E+3             0.     7.0224E+02         1.0E-3
%##
%##  DATASET 10:  Solid Matrix Properties
%##     [COMPMA]           [CS]       [SIGMAS]         [RHOS]
%         0.0             0.             0.             1.
%##
%##  DATASET 11:  Adsorption Parameters


inp.ncolpr = -10;
%inp.ncol  = 'N'  'X'  'Y'  'P'  'U'  'S'  '-';
inp.ncol   = {['N'],['X' ],[ 'Y'  ],['P' ],[ 'U' ],[ 'S' ],[ '-']};

inp.lcolpr = 10;
inp.lcol   = {[ 'E' ],[ 'X' ],[ 'Y'  ],['VX' ],['VY' ],['-']};


inp.nbcfpr = 10;
inp.nbcspr = 10;
inp.nbcppr = 10;
inp.nbcupr = 10;
inp.cinact = 'Y';



%##    DATASET 9:  FLUID PROPERTIES
inp.compfl = 1.e-11;
inp.cw     = 0;
inp.sigmaw = 1.e-9;
inp.rhow0  = 1000;
inp.urhow0 = 0;
inp.drwdu  = 700;
inp.visc0  = 0.001;


%##      DATASET 10:  SOLID MATRIX PROPERTIES
inp.compma = 1.e-7;
inp.cs     = 0;
inp.sigmas = 0;
inp.rhos   = 2600;   %solid density of sodium chloride

%##  DATASET 11:  Adsorption Parameters
%##     [ADSMOD]         [CHI1]         [CHI2]
%#'NONE'
%'FREUNDLICH' 1.D-47 0.05  #less rigid

%inp.adsmod = 'FREUNDLICH';
inp.adsmod = 'FREUNDLICH';
inp.chi1   = 1.D-46;
inp.chi2   = 0.05 ;



%##
%##  DATASET 12:  Production of Energy or Solute Mass
%##     [PRODF0]       [PRODS0]       [PRODF1]       [PRODS1]
%0.             0.             0.             0.
inp.prodf0 = 0;
inp.prods0 = 0;
inp.prodf1 = 0;
inp.prods1 = 0;


%##
%##  DATASET 13:  Orientation of Coordinates to Gravity
%##      [GRAVX]        [GRAVY]        [GRAVZ]
%0.           -9.81          0.
inp.gravx = 0;
inp.gravy = -9.81;
inp.gravz = 0;

%##  DATASET 14:  NODEWISE DATA
%%##                              [SCALX] [SCALY] [SCALZ] [PORFAC]
inp.scalx  = 1.;
inp.scaly  = 1.;
inp.scalz  = 1.;
inp.porfac = 1.;

%##      [II]    [NRE    G(II)]  [X(II)] [Y(II)] [Z(II)] [POR(II)]
inp.ii   = (1:nn)';
inp.nreg = zeros(nn,1)+1;
inp.x    = x_nod_array';
inp.y    = y_nod_array';
inp.z    = zeros(nn,1)+1; % FOR 2D mesh, z does not matter, here set as 1
% inp.z    = 2*pi*inp.x+0.1;
inp.por  = zeros(nn,1)+0.43;

% dataset 15: ELEMENTWISE DATA
%##                              [PMAXFA]        [PMINFA]        [ANG1FA]        [ALMAXF]        [ALMINF]        [ATMAXF]        [ATMINF]
%'ELEMENT'               1.0000000D+00   1.0000000D+00   1.0000000D+00   2 2 2 2
%##     [L]      [LREG(L)]       [PMAX(L)]       [PMIN(L)]       [ANGLEX(L)]     [ALMAX(L)]      [ALMIN(L)]      [ATMAX(L)]      [ATMIN(L)]
inp.pmaxfa = 1.;
inp.pminfa = 1.;
inp.angfac = 0.;
inp.almaxf = 1.;
inp.alminf = 1.;
inp.atmaxf = 1.;
inp.atminf = 1.;
inp.l      = (1:ne)';
inp.lreg   = zeros(ne,1)+1;

pmax_mtx = zeros(size(x_ele_mtx_m_g))+permeability_sand_m2;
% consider the following codes for layered aquifer
% mask_ele_mtx_silt_layer_gravity_compensated = y_ele_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for element matrix 
% 
% pmax_mtx_gravity_compensated_m2 (mask_ele_mtx_silt_layer_gravity_compensated) =  permeability_silt_m2;   % silt layer permeability
% 
% pmax_mtx_m2=flip(pmax_mtx_gravity_compensated_m2);
pmax_array = pmax_mtx(:);

inp.pmax   = pmax_array;
inp.pmin   = pmax_array;
inp.anglex = zeros(ne,1);
inp.almax  = zeros(ne,1)+0.5e-0;
inp.almin  = zeros(ne,1)+0.5e-0;
inp.atmax  = zeros(ne,1)+0.5e-1;
inp.atmin  = zeros(ne,1)+0.5e-1;

% dataset 17: FLUID SOURCES AND SINKS, used for ET in SUTRASET
% note: for SUTRASET when node number is negative, the second input is
% surface area of the node (which is dx while z  = 1);
% and the third input is the thickness (dy)of the cell.

inp.iqcp  = -nod_idx_mtx_g(1,:)';
inp.qinc  = dx_cell_mtx_g(1,:)'*1; % dx*dz
inp.uinc  = dy_cell_mtx_g(1,:)'; %dy
inp.nsop = length(inp.iqcp); % here to change the initial NSOP 0 value to actual size. this has to be changed or this dataset cannot be printed

% dataset 18: Energy or Solute sources and sinks

% dataset 19: Specified pressure
inp.ipbc = -nod_idx_mtx_g(1,:)'; % for tide pressure
inp.pbc  = zeros(size(inp.ipbc)) - 7.37325e+03; % note that even though inp.ipbc is designed as negative, inp.pbc, .ubc, .npbc has to be specified in order to print this dataset
inp.ubc  = zeros(size(inp.ipbc)) + c_saltwater_kgPkg;
inp.npbc = length(inp.ipbc); % here to change the initial given 0 size to actual size

%##  DATASET     22:  Ele        ment Incid      ence Data
%##    [LL]      [IIN(1)]        [IIN(2)]        [IIN(3)]        [IIN(4)]

%ne_mtx=reshape(l,ney,nex);


% DATASET 22, user need to check if the sequence of the node is in clockwise manner as suggested by the manual

inp.iin1= iin1;
inp.iin2= iin2;
inp.iin3= iin3;
inp.iin4= iin4;

%% WRITE .INP file
inp.export_to_file();

%% ICS FILE

pm1_mtx_pa = 1024.99*9.84*(initial_head_aquifer_m-y_nod_mtx);
% um1_mtx_kgPkg = zeros(size(nod_idx_mtx_g))+c_saltwater_kgPkg;

ics       = icsObj('SAND','read_from_file','no');
ics.tics  = 0;
%ics.cpuni = 'UNIFORM';
%ics.pm1   = -30;
ics.cpuni = 'NONUNIFORM';
ics.pm1   = pm1_mtx_pa(:);
ics.cuuni = 'UNIFORM';
ics.um1   = c_saltwater_kgPkg ;

ics.export_to_file();