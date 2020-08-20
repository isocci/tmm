% Set toplevel path to GCMs configuration
base_path='~/Desktop/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicForcing=1
periodicMatrix=1

dt=43200; % time step to use

rearrangeProfiles=1
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0
writePCFiles=0

READ_SWRAD=0
useCarbon=1
useAtmModel=1
pCO2atm_ini=280.0
useVirtualFlux=0
empScaleFactor=1.0

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

explicitMatrixFileBase=fullfile(base_path,explicitMatrixFileBase);
implicitMatrixFileBase=fullfile(base_path,implicitMatrixFileBase);

explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);

preconditionerMatrixFile=fullfile(base_path,preconditionerMatrixFile);

gcmDataPath=fullfile(base_path,'GCM');
bgcDataPath=fullfile(base_path,'BiogeochemData');
freshWaterForcingFile=fullfile(gcmDataPath,'FreshWaterForcing_gcm');
empFixFile=fullfile(gcmDataPath,empFixFile);

iceFile=fullfile(bgcDataPath,'ice_fraction');
windFile=fullfile(bgcDataPath,'wind_speed');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','dz','da','x','y','z','deltaT','gridType')

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
nbb=length(Ib);

if rearrangeProfiles || bigMat
  load(profilesFile,'Ip_pre','Ir_pre','Ip_post','Ir_post','Irr')
  Ip=Ip_pre;
  Ir=Ir_pre;
end

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

if periodicForcing
  nm=12;
else
  nm=1;
end

% Use steady state T/S from GCM. Note we always load seasonal data here.
load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);

if ~useCarbon
  useAtmModel=0;
  useVirtualFlux=0;
end

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);

% Surface forcing
% Compute total E-P for virtual flux:
% Fv = TRg*(E-P), E-P in m/s
% d[TR]/dt = ... + Fv/dz
load(freshWaterForcingFile,'EmPgcm','Srelaxgcm','saltRelaxTimegcm')
zeroNetEmP=1;
EmP=get_surface_emp_for_virtual_flux(gridFile,boxFile,EmPgcm,Srelaxgcm,saltRelaxTimegcm,Salt,zeroNetEmP);
if fixEmP
  load(empFixFile,'empFixX','empFixY')
  nEmPFix=length(empFixX);
% Points to fix: each SET is individually fixed.
  for k=1:nEmPFix
	if length(empFixX{k})>1 % polygon  
	  Ifix{k}=find(inpolygon(Xboxnom(Ib),Yboxnom(Ib),empFixX{k},empFixY{k})); % indexed to Ib
	else % single point
	  Ifix{k}=find(Xboxnom(Ib)==empFixX{k} & Yboxnom(Ib)==empFixY{k}); % referenced to Ib
	end
  end
% Rest of ocean (not fixed)
  Infix=find(~ismember([1:nbb]',cat(1,Ifix{:}))); % Everything else indexed to Ib	
% Now fix E-P so that annual and spatial integral is 0.
  for k=1:length(Ifix) % loop over each SET of problematic points
	if useAreaWeighting
	  areabfix=dab_surf(Ifix{k});
	  EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(areabfix'*EmP(Ifix{k},:))/sum(areabfix);
	else
	  EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(mean(EmP(Ifix{k},:)));
	end
  end
  if useAreaWeighting
	areabnfix=dab_surf(Infix);
	EmP(Infix,:) = EmP(Infix,:) - mean(areabnfix'*EmP(Infix,:))/sum(areabnfix);
  else
	EmP(Infix,:) = EmP(Infix,:) - mean(mean(EmP(Infix,:)));
  end
end
EmP=EmP*empScaleFactor;

if useVirtualFlux
  volFracSurf=zeros(nb,1);
  volFracSurf(Ib)=volb(Ib)/sum(volb(Ib));  
end

% Grid variables
dzb=gridToMatrix(dz,[],boxFile,gridFile);

% Surface forcing data
load(iceFile,'Fice')
load(windFile,'windspeed')
Ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
windb=gridToMatrix(windspeed,Ib,boxFile,gridFile,1);

atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
  EmP=mean(EmP,2);
  Ficeb=mean(Ficeb,2);
  windb=mean(windb,2);  
  atmospb=mean(atmospb,2);
end

if ~periodicMatrix
  if rescaleForcing
	Rfs=mean(Rfs,2);    
  end  
end

if READ_SWRAD
% Read SW radiation from file
  swradb=repmat(1000,[nbb 1]);
else
% Interpret par as latitude
  latb=Yboxnom(Ib);
end

% Some constants for carbonate chemistry
rho0=1024.5; % kg/m^3 nominal density from OCMIP-2

if rearrangeProfiles
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  Theta=Theta(Ir,:);
  Salt=Salt(Ir,:);
  if useVirtualFlux
    volFracSurf=volFracSurf(Ir);
  end  
  if rescaleForcing
    Rfs=Rfs(Ir,:);
  end  
  dzb=dzb(Ir);
  Ib=find(izBox==1);
%
  Ip=Ip_post;
  Ir=Ir_post;
end  

% Initial condition
PO4=repmat(2.17,[nb 1]); % [mmol P/m3]
DOP=repmat(1e-4,[nb 1]); % [mmol P/m3]
OXY=repmat(300.0,[nb 1]); % [mmol O2/m3]
PHY=repmat(1e-4,[nb 1]); % [mmol P/m3]
ZOO=repmat(1e-4,[nb 1]); % [mmol P/m3]
DET=repmat(1e-4,[nb 1]); % [mmol P/m3]
NO3=repmat(31.23,[nb 1]); % [mmol N/m3]
if useCarbon
  DIC=repmat(2251,[nb 1]); % [mmol C/m3]
  ALK=repmat(2360,[nb 1]); % [mmol eq/m3]
end

if useCoarseGrainedMatrix
% Coarse grain initial conditions
end

if writeFiles
  calc_periodic_times_for_tmm('monthly-365-day year','periodic_times_365d.bin');
  calc_periodic_times_for_tmm('monthly-360-day year','periodic_times_360d.bin');  
% Transport matrices
  if writeTMs
%   Explicit transport matrix
	I=speye(nb,nb);
	if ~periodicMatrix
      disp('loading annual mean explicit TM')	
      load(explicitAnnualMeanMatrixFile,'Aexpms')	
	  if rearrangeProfiles
		Aexpms=Aexpms(Ir_pre,Ir_pre); % rearrange
	  end      
	  % make discrete
	  Aexpms=dt*Aexpms;
	  Aexpms=I+Aexpms;
	  if useCoarseGrainedMatrix
		Aexpms=Beta*Aexpms*M; % coarse-grained explicit transport matrix
	  end        
	  writePetscBin('Ae.petsc',Aexpms,[],1)
	else
      % load each month from separate file
      disp('loading monthly mean explicit TMs')	      
	  for im=1:12 
		fn=[explicitMatrixFileBase '_' sprintf('%02d',im)];
		load(fn,'Aexp')
		if rearrangeProfiles
		  Aexp=Aexp(Ir_pre,Ir_pre); % rearrange
		end
		% make discrete
		Aexp=dt*Aexp;
		Aexp=I+Aexp;
		if useCoarseGrainedMatrix
%         Not sure if this is really kosher!		  
		  Aexp=Beta*Aexp*M; % coarse-grained explicit transport matrix
		end
		writePetscBin(['Ae_' sprintf('%02d',im-1)],Aexp,[],1)
		clear Aexp
	  end
	end
%   Implicit transport matrix
	if ~periodicMatrix
      disp('loading annual mean implicit TM')		
      load(implicitAnnualMeanMatrixFile,'Aimpms')
      if dtMultiple~=1
		if bigMat % big matrix. do it a block at a time.
		  for is=1:nbb % change time step multiple
			Aimpms(Ip_pre{is},Ip_pre{is})=Aimpms(Ip_pre{is},Ip_pre{is})^dtMultiple;
		  end
		else
		  Aimpms=Aimpms^dtMultiple;
		end  
	  end	
	  if rearrangeProfiles
		Aimpms=Aimpms(Ir_pre,Ir_pre); % rearrange
	  end
	  if useCoarseGrainedMatrix
		Aimpms=Beta*Aimpms*M; % coarse-grained implicit transport matrix
	  end
	  writePetscBin('Ai.petsc',Aimpms,[],1)
	else
	  % load each month from separate file
      disp('loading monthly mean implicit TMs')	      	  
	  for im=1:12
		fn=[implicitMatrixFileBase '_' sprintf('%02d',im)];		
		load(fn,'Aimp')
		if dtMultiple~=1
		  if bigMat % big matrix. do it a block at a time.		
			for is=1:nbb % change time step multiple
			  Aimp(Ip_pre{is},Ip_pre{is})=Aimp(Ip_pre{is},Ip_pre{is})^dtMultiple;
			end
		  else
			Aimp=Aimp^dtMultiple;		
		  end
		end  
		if rearrangeProfiles
		  Aimp=Aimp(Ir_pre,Ir_pre); % rearrange
		end
		if useCoarseGrainedMatrix
		  Aimp=Beta*Aimp*M; % coarse-grained implicit transport matrix		
		end
		writePetscBin(['Ai_' sprintf('%02d',im-1)],Aimp,[],1)
		clear Aimp
	  end
	end
  end	  	  
% Initial conditions  
  writePetscBin('po4ini.petsc',PO4)
  writePetscBin('dopini.petsc',DOP)
  writePetscBin('oxyini.petsc',OXY)
  writePetscBin('phyini.petsc',PHY)
  writePetscBin('zooini.petsc',ZOO)  
  writePetscBin('detini.petsc',DET)
  writePetscBin('no3ini.petsc',NO3)  
  if useCarbon
	writePetscBin('dicini.petsc',DIC)
	writePetscBin('alkini.petsc',ALK)
  end    
  if ~periodicForcing
	write_binary('fice.bin',Ficeb,'real*8')
	write_binary('wind.bin',windb,'real*8')	
	write_binary('atmosp.bin',atmospb,'real*8')	
  else
    for im=1:nm
	  write_binary(['fice_' sprintf('%02d',im-1)],Ficeb(:,im),'real*8')
	  write_binary(['wind_' sprintf('%02d',im-1)],windb(:,im),'real*8')	  
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')	  
	end
  end
  if READ_SWRAD
	write_binary('swrad.bin',swradb,'real*8')
  else   
    write_binary('latitude.bin',latb,'real*8')
  end
  if ~periodicForcing
	writePetscBin('Ts.petsc',Theta)
	writePetscBin('Ss.petsc',Salt)
  else
    for im=1:nm
	  writePetscBin(['Ts_' sprintf('%02d',im-1)],Theta(:,im))
	  writePetscBin(['Ss_' sprintf('%02d',im-1)],Salt(:,im))
    end    
  end
  if ~periodicForcing
	write_binary('EmP.bin',EmP,'real*8')
  else
	for im=1:nm
	  write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
	end
  end        
  if useVirtualFlux
    writePetscBin('surface_volume_fraction.petsc',volFracSurf)
  end  
  if rescaleForcing
	if ~periodicMatrix
	  writePetscBin('Rfs.petsc',Rfs)
	else
	  for im=1:nm
		writePetscBin(['Rfs_' sprintf('%02d',im-1)],Rfs(:,im))
	  end    
	end    
  end  
% Grid data
  write_binary('drF.bin',nz,'int')  
  write_binary('drF.bin',dznom,'real*8',1)
  writePetscBin('dz.petsc',dzb)
  if useAtmModel
    write_binary('pCO2atm_ini.bin',pCO2atm_ini,'real*8')
    write_binary('dA.bin',dab_surf,'real*8')    
  end  
% Profile data  
  if rearrangeProfiles
    if ~useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ip);
	  gEndIndices=cellfun(@(x)x(end),Ip);
    else % useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ipcg);
	  gEndIndices=cellfun(@(x)x(end),Ipcg);
    end  
    write_binary('gStartIndices.bin',[length(gStartIndices);gStartIndices],'int')
    write_binary('gEndIndices.bin',[length(gEndIndices);gEndIndices],'int')
  end
end

if useCoarseGrainedMatrix
  numProfiles=nbbcg;
else  
  numProfiles=nbb;
end
disp(['Number of Profiles in this Configuration: ' int2str(numProfiles)])

if writePCFiles
  pc=load(preconditionerMatrixFile,'Aexpms');
  if rearrangeProfiles
    A=pc.Aexpms(Ir_pre,Ir_pre);
  else
    A=pc.Aexpms;
  end
  clear pc  
  if useCoarseGrainedMatrix
    A=Beta*A*M;
    save pc_cg_data A nbbcg CGgrid CG Ipcg Ibcg dt
  else
    save pc_data A nbb nz nb Ip Ib dt
  end
end
