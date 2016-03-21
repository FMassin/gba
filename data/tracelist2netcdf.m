clear all;
 
% Load the data
load zTrainingHonshu08.mat;
load hTrainingHonshu08.mat;

% Define the time steps
tmin = 0.5;
tmax = 20;
dt = 0.5;
t = tmin:dt:tmax;

% Define the filter coordinate axes as the upper corner frequency of 
% each bandpass filter
f = [0.1875,0.375,0.75,1.5,3,6,12,24,48];

% Determine the array shape
ntraces = length(zTraining.amax);
[nfilter,ntimes] = size(zTraining.amax{1});
ntracesh = length(hTraining.amax);
[nfilterh,ntimesh] = size(hTraining.amax{1});

if nfilter ~= nfilterh || ntimes ~= ntimesh || ntraces ~= ntracesh
    fprintf(1,'Dimensions in Z and H training set not identical')
    exit
end

if length(t) ~= ntimes
    fprintf(1,'Time axes inconsistent with training data')
    exit
end

ztrain = zeros(ntimes,ntraces,nfilter);
htrain = zeros(ntimes,ntraces,nfilter);

% Convert cell arrays into 3 dimensional matrices
for k=1:ntraces
    ztrain(:,k,:) = log10(zTraining.amax{k})';
    htrain(:,k,:) = log10(hTraining.amax{k})';
end

% Write everything to a netcdf file
ncid = netcdf.create('GbA_training.nc','NC_CLOBBER');
dimid1 = netcdf.defDim(ncid,'time',ntimes);
dimid2 = netcdf.defDim(ncid,'traces',ntraces);
dimid3 = netcdf.defDim(ncid,'filter',nfilter);

% Describe the coordinate axes
varid_t = netcdf.defVar(ncid,'time','double',[dimid1]);
netcdf.putAtt(ncid,varid_t,'units','s');
varid_f = netcdf.defVar(ncid,'filter','double',[dimid3]);
netcdf.putAtt(ncid,varid_f,'units','Hz');
varid_m = netcdf.defVar(ncid,'magnitude','double',[dimid2]);
varid_ed = netcdf.defVar(ncid,'epicdist','double',[dimid2]);
netcdf.putAtt(ncid,varid_ed,'units','km');
varid_hd = netcdf.defVar(ncid,'hypodist','double',[dimid2]);
netcdf.putAtt(ncid,varid_hd,'units','km');

% Matlab stores matrices in column major format, meaning the 
% 1st index varies the fastest and the last the slowest. So to have the 
% the same results for the same indices in C++, Python, and Matlab, the 
% dimensions have to be reversed before writing.
varidz = netcdf.defVar(ncid,'z','double',[dimid3,dimid2,dimid1]);
varidh = netcdf.defVar(ncid,'h','double',[dimid3,dimid2,dimid1]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varidz,permute(ztrain,[3,2,1]));
netcdf.putVar(ncid,varidh,permute(htrain,[3,2,1]));

netcdf.putVar(ncid,varid_t,t);
netcdf.putVar(ncid,varid_f,f);
netcdf.putVar(ncid,varid_m,zTraining.m);
netcdf.putVar(ncid,varid_ed,zTraining.epiDist);
netcdf.putVar(ncid,varid_hd,zTraining.hypDist);

netcdf.close(ncid);
