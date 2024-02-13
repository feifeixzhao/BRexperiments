% caesar_discharge_file.m: Creates text file for discharge timeseries to
% feed to CAESAR model. 

%%% Edit these parameters
filename = 'run2_scan_015_discharge.txt';
linear_Q_ramp = false;
Qinit = 0.12e-3; % initial discharge, m^3/s 1e-9 is a good place to start
Qfinal = 0.12e-3;
% Qfinal = 0.247e-3; % final discharge, m^3/s  % 11e-5 looks like the one for full inundation.
Qtimestep_minutes = 60; % use this because it's roughly the equilibration time and makes run time in hours have a decimal of finite length (needed for caesar)
timesteps_initialize = 2; % number of timesteps at the start of the run to keep the discharge fixed at the initial value. Should be at least 2 to leave time to fill reservoir.
N_Q_levels = 10; % number of discharge levels between Qinit and Qfinal, inclusive.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if linear_Q_ramp
    Q_all = [repmat(Qinit,timesteps_initialize,1);linspace(Qinit,Qfinal,N_Q_levels)'];
else
    Q_all = [repmat(Qinit,timesteps_initialize,1);logspace(log10(Qinit),log10(Qfinal),N_Q_levels)'];
end
timestep = ((1:numel(Q_all))-1)'; % i.e., starting from 0
% Arrange in a column as timestep-Q; timestep-Q; ...
tQ = [timestep,Q_all];
tQ = tQ';
tQ = tQ(:);
fid = fopen(filename,'wt'); % for text files, use 'wt', it works better than just 'w'
formatspec = '%1.0f %5.12f 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n';
fprintf(fid,formatspec,tQ);
fclose(fid);
type(filename)

sprintf('Run time = %3.2f hours',max(timestep)*Qtimestep_minutes/60)
