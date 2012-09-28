%FSEM3D_SNAPSHOT reads and plots fault fields at a given time
%
% d = FSEM3D_snapshot(isnap, [fig, dir, fault])
%
% INPUTS	isnap	snapshot index, as in Snapshot*.bin file names
%		fig	[1] plot horizontal slip rate or not
%		        *** this option has been temporarily disabled, the value is ignored ***
%		dir	["."] directory containing the SPECFEM3D output data Snapshot*.bin
%		fault	[1] fault id
%
% OUTPUTS	d	structure containing fault fields:
%			X,Y,Z	  fault node coordinates (km)
%			Dx,Dz	  slip (m)
%			Vx,Vz	  slip rate (m/s)
%			Tx,Ty,Tz  stress change (MPa)
%			S	  slip "state" variable in friction law (m)
%			Sg	  strength relative to initial stress (MPa)
%			Trup	  rupture time (s)
%			Tpz	  process zone time, when slip=Dc (s)
%
% Jean-Paul Ampuero	ampuero@erdw.ethz.ch modified by 
% Percy Galvez  percy.galvez@sed.ethz.ch 19/01/2011.
%
% WARNING : This script only do snapshots of single precision files so far.

function d = FSEM3D_snapshot(isnap,fig,DATA_DIR,fault)

NDAT = 14; 
VM = 3;	% saturate the color scale of the slip rate figure at VM m/s

if nargin<2, fig=1; end
if nargin<3, DATA_DIR = '.'; end
if nargin<4, fault = 1; end

BinFile = sprintf('%s/Snapshot%u_F%u.bin',DATA_DIR,isnap,fault);

if ~exist(BinFile,'file'), error(sprintf('File %s does not exist',BinFile)), end
fid=fopen(BinFile);
BinRead = fread(fid,[1,inf],'single')' ;
fclose(fid);

BinRead = reshape( BinRead(:),length(BinRead)/(NDAT),NDAT);
BinRead = BinRead(2:end-1,:);

% Reorder fault nodes (lexicographic z,x)
%[LOC,IND] = sortrows( BinRead(:,[1 3]) );
%BinRead = BinRead(IND,:);

d.X  = BinRead(:,1)/1e3; % in km
d.Y  = BinRead(:,2)/1e3; % in km
d.Z  = BinRead(:,3)/1e3; % in km
d.Dx = BinRead(:,4);
d.Dz = BinRead(:,5);
d.Vx = BinRead(:,6);
d.Vz = BinRead(:,7);
d.Tx = BinRead(:,8); % in MPa
d.Ty = BinRead(:,9);
d.Tz = BinRead(:,10); % in MPa
d.S  = BinRead(:,11);
d.Sg = BinRead(:,12); % in MPa
d.Trup = BinRead(:,13);
d.Tpz = BinRead(:,14);

clear BinRead


return
