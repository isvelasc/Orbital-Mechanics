clc; clear; close all;

addpath('..\Functions\');
addpath('..\Functions\SGP4\');

%% Part 1 Info
latitude  = 40; % deg
longitude = -60; % deg
altitude  = -0.0269114; % km

%% Part 1 
tleinfo  = readTLE('ISS TLE.txt'); % acquired from https://celestrak.com/NORAD/elements/stations.txt
tleinfo.epoch(1) = tleinfo.epoch(1) + 2000; % fix date component
viewtime = datevec('26 Jan 2022 02:03:00.000');
time     = etime(viewtime,tleinfo.epoch);
TLEfname = TLEpropagator('ISS TLE.txt',time,'f',0);
tleinfo  = readTLE(TLEfname);

coes     = coesfromTLE(tleinfo);
[r, v]   = coes2sv(coes.h, coes.ecc, coes.inc, coes.raan, coes.w, coes.TA, 398600);
tleinfo.epoch(1) = tleinfo.epoch(1) + 2000; % fix date component

% sidereal time
lst      = LST(viewtime, longitude); % deg

% site location
Rsite    = RS_ECI(deg2rad(latitude),altitude,deg2rad(lst)); % km

% Az/el
[Az,el,rho] = sv2razel(r',Rsite',lst,latitude);

% right ascension and declination
[RA,DEC] = RA_DEC_from_Az_el(latitude,Az,el,lst,'d');
 


%% Part 2 Info
latitude  = 35.25402; % deg
longitude = -120.69235; % deg
altitude  = 42/1000; % km

%% Part 2
tleinfo   = readTLE('CZ-4 RB.txt'); % acquired from www.heavens-above.com
tleinfo.epoch(1) = tleinfo.epoch(1) + 2000;

% set start time prior to S/C viewing time
viewtime  = datetime('6-Feb-2022 18:04:31','TimeZone','America/Los_Angeles');
viewtime.TimeZone = 'UTC';
vtstart = datevec(viewtime);

% do first propagation to prior
time     = etime(vtstart,tleinfo.epoch);
TLEfname = TLEpropagator('CZ-4 RB.txt',time,'f',0);

% set end time post to S/C viewing time
viewtime  = datetime('6-Feb-2022 20:04:31','TimeZone','America/Los_Angeles');
viewtime.TimeZone = 'UTC';
vtend     = datevec(viewtime);

% do second propagation from prior to post
time           = etime(vtend,vtstart);
[TLEfname,t,s] = TLEpropagator(TLEfname,time,'f',1);

% make time vector for viewtime window
dtview = linspace(datetime(vtstart),datetime(vtend),length(t));
for tstep = 2:length(t)-1
    dtview(tstep) = dtview(1) + seconds(t(tstep)); % fix times to correspond to each ode45 iteration
end
dtview = datevec(dtview);

criteria = zeros(length(t),1);
Az       = zeros(length(t),1);
el       = zeros(length(t),1);
RA       = zeros(length(t),1);
DEC      = zeros(length(t),1);
% compute criteria for each iteration
for tstep = 1:length(t)

    % preliminary data
    lst    = LST(dtview(tstep,:),longitude);
    Rsite  = RS_ECI(deg2rad(latitude),altitude,deg2rad(lst));
    [Az(tstep),el(tstep),rho,rhosez] = sv2razel(s(tstep,1:3),Rsite',lst,latitude);
    rhosez = rhosez/rho;
    [RA(tstep),DEC(tstep)] = RA_DEC_from_Az_el(latitude,Az(tstep),el(tstep),lst,'d');
    
    % convert position vector to ECEF
    tempvec = rotz(lst)\s(tstep,1:3)';
    recef   = roty(90-latitude)\tempvec;

    % get solar vector and convert to ECEF
    suneci  = solar_position(juliandate(dtview(tstep,:)));
    tempvec = rotz(lst)\suneci;
    sunecef = roty(90-latitude)\tempvec;

    % get distance for viewability
    phi  = asin(cross(sunecef,recef)/(norm(sunecef)*norm(recef)));
    dist = norm(norm(recef)*cos(phi - pi/2));
    
    % calculate three criteria: 1) above horizon 2) site in darkness 3) S/C visible
    if rhosez(3) > 0 && el(tstep) > 10 && dot(suneci,Rsite) < 0 && dist > 6378
        criteria(tstep) = 1;
    else
        criteria(tstep) = 0;
    end

end

% find times where all criteria was met
idx = find(criteria == 1);
viewtimes = datetime(dtview(idx,:),'TimeZone','UTC');
viewtimes.TimeZone = 'America/Los_Angeles';
vttable = table(viewtimes,Az(idx),el(idx),RA(idx),DEC(idx),'VariableNames',{'Viewing Window','Azimuth','Elevation','Right Ascension','Declination'});
disp(vttable)



%% Part 3 Info
latitude        = 35.30; % deg
longitude       = -120.66; % deg
altitude        = 105.8/1000; % km
single_pass     = readtable('Single Pass Observation.txt');
double_pass     = readtable('Double Pass Observation.txt');
single_passTLE1 = readTLE('Single Unknown 1.txt');
single_passTLE2 = readTLE('Single Unknown 2.txt');
single_passTLE3 = readTLE('Single Unknown 3.txt');
double_passTLE1 = readTLE('Double Unknown 1.txt');
double_passTLE2 = readTLE('Double Unknown 2.txt');
double_passTLE3 = readTLE('Double Unknown 3.txt');
% update TLE epoch
single_passTLE1.epoch(1) = single_passTLE1.epoch(1) + 2000;
single_passTLE2.epoch(1) = single_passTLE2.epoch(1) + 2000;
single_passTLE3.epoch(1) = single_passTLE3.epoch(1) + 2000;
double_passTLE1.epoch(1) = double_passTLE1.epoch(1) + 2000;
double_passTLE2.epoch(1) = double_passTLE2.epoch(1) + 2000;
double_passTLE3.epoch(1) = double_passTLE3.epoch(1) + 2000;

%% Part 3 Single Pass
sp_alpha = single_pass.RA_degs; % deg
sp_delta = single_pass.DEC_degs; % deg
sp_date  = single_pass.Date_UT{1};
sp_time  = [datevec(strcat(sp_date," ",string(single_pass.Time_UT(1))));...
            datevec(strcat(sp_date," ",string(single_pass.Time_UT(2))));...
            datevec(strcat(sp_date," ",string(single_pass.Time_UT(3))))];

% sidereal time
sp_lst   = [LST(sp_time(1,:),longitude);...
            LST(sp_time(2,:),longitude);...
            LST(sp_time(3,:),longitude)]; % deg

% site location
sp_Rsite = [RS_ECI(deg2rad(latitude),altitude,deg2rad(sp_lst(1)))...
            RS_ECI(deg2rad(latitude),altitude,deg2rad(sp_lst(2)))...
            RS_ECI(deg2rad(latitude),altitude,deg2rad(sp_lst(3)))];

% distance to object
sp_rhoh  = [rhohat(sp_alpha(1),sp_delta(1),'d')...
            rhohat(sp_alpha(2),sp_delta(2),'d')...
            rhohat(sp_alpha(3),sp_delta(3),'d')];

% time differences
sp_tau1  = etime(sp_time(1,:),sp_time(2,:)); % s
sp_tau3  = etime(sp_time(3,:),sp_time(2,:)); % s

% Dr. Kira Abercromby Code
% Double R Method
[sp_rvec2, sp_vvec2] = AERO557doubleR(sp_rhoh(:,1),sp_rhoh(:,2),sp_rhoh(:,3),...
                    sp_Rsite(:,1),sp_Rsite(:,2),sp_Rsite(:,3),sp_tau1,sp_tau3); % [km km/s]


% back poropagate to to t1
options  = odeset('RelTol',1e-8,'AbsTol',1e-8);
sp_dtb   = [0 sp_tau1];
[~,sp1]  = ode45(@orbitprop,sp_dtb,[sp_rvec2' sp_vvec2'],options,398600);
sp_coes1 = sv2coes(sp1(end,1:3),sp1(end,4:6));

% back propagate TLE's to t1
time1     = etime(sp_time(1,:),single_passTLE1.epoch);
time2     = etime(sp_time(1,:),single_passTLE2.epoch);
time3     = etime(sp_time(1,:),single_passTLE3.epoch);
TLEfname1 = TLEpropagator('Single Unknown 1.txt',time1,'f',0);
TLEfname2 = TLEpropagator('Single Unknown 2.txt',time2,'f',0);
TLEfname3 = TLEpropagator('Single Unknown 3.txt',time3,'f',0);
tleinfo1  = readTLE(TLEfname1);
tleinfo2  = readTLE(TLEfname2);
tleinfo3  = readTLE(TLEfname3);
coes1     = coesfromTLE(tleinfo1);
coes2     = coesfromTLE(tleinfo2);
coes3     = coesfromTLE(tleinfo3);


% compare COEs from TLE specs and observations
obs1    = [rad2deg(sp_coes1.inc);rad2deg(sp_coes1.raan);sp_coes1.ecc;rad2deg(sp_coes1.omega);rad2deg(sp_coes1.Me);sqrt(398600/sp_coes1.a^3)*24*3600/(2*pi);sp_coes1.a];
TLE1    = [rad2deg(coes1.inc);rad2deg(coes1.raan);coes1.ecc;rad2deg(coes1.w);rad2deg(coes1.Me);coes1.N*24*3600/(2*pi);coes1.a];
TLE2    = [rad2deg(coes2.inc);rad2deg(coes2.raan);coes2.ecc;rad2deg(coes2.w);rad2deg(coes2.Me);coes2.N*24*3600/(2*pi);coes2.a];
TLE3    = [rad2deg(coes3.inc);rad2deg(coes3.raan);coes3.ecc;rad2deg(coes3.w);rad2deg(coes3.Me);coes3.N*24*3600/(2*pi);coes3.a];
col_n   = {'Obs','TLE 1','TLE 2','TLE 3'};
row_n   = {'inc  (deg)','raan (deg)','ecc','argp (deg)','me   (deg)','n    (rev/day)','a    (km)'};
sp_comp = table(obs1,TLE1,TLE2,TLE3,'VariableNames',col_n,'RowNames',row_n);
disp(sp_comp)


%% Part 3 Double Pass
dp1_alpha = double_pass.RA_degs(1:4); % deg
dp2_alpha = double_pass.RA_degs(5:8); % deg
dp1_delta = double_pass.DEC_degs(1:4); % deg
dp2_delta = double_pass.DEC_degs(5:8); % deg
dp_date   = double_pass.Date_UT{1};
dp1_time  = [datevec(strcat(dp_date," ",string(double_pass.Time_UT(1))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(2))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(3))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(4))))];

dp2_time  = [datevec(strcat(dp_date," ",string(double_pass.Time_UT(5))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(6))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(7))));...
             datevec(strcat(dp_date," ",string(double_pass.Time_UT(8))))];

% sidereal time
dp1_lst   = [LST(dp1_time(1,:),longitude);...
             LST(dp1_time(2,:),longitude);...
             LST(dp1_time(3,:),longitude);...
             LST(dp1_time(4,:),longitude)]; % deg
dp2_lst   = [LST(dp2_time(1,:),longitude);...
             LST(dp2_time(2,:),longitude);...
             LST(dp2_time(3,:),longitude);...
             LST(dp2_time(3,:),longitude)]; % deg

% site location
dp1_Rsite = [RS_ECI(deg2rad(latitude),altitude,deg2rad(dp1_lst(1)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp1_lst(2)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp1_lst(3)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp1_lst(4)))];
dp2_Rsite = [RS_ECI(deg2rad(latitude),altitude,deg2rad(dp2_lst(1)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp2_lst(2)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp2_lst(3)))...
             RS_ECI(deg2rad(latitude),altitude,deg2rad(dp2_lst(4)))];

% distance to object
dp1_rhoh  = [rhohat(dp1_alpha(1),dp1_delta(1),'d')...
             rhohat(dp1_alpha(2),dp1_delta(2),'d')...
             rhohat(dp1_alpha(3),dp1_delta(3),'d')...
             rhohat(dp1_alpha(4),dp1_delta(4),'d')];
dp2_rhoh  = [rhohat(dp2_alpha(1),dp2_delta(1),'d')...
             rhohat(dp2_alpha(2),dp2_delta(2),'d')...
             rhohat(dp2_alpha(3),dp2_delta(3),'d')...
             rhohat(dp2_alpha(4),dp2_delta(4),'d')];

% time differences for middle two passes 
dp1a_tau1  = etime(dp1_time(1,:),dp1_time(2,:)); % s \ for r2,v2
dp1a_tau3  = etime(dp1_time(3,:),dp1_time(2,:)); % s /
dp2a_tau1  = etime(dp2_time(1,:),dp2_time(2,:)); % s \ for r6,v6
dp2a_tau3  = etime(dp2_time(3,:),dp2_time(2,:)); % s /

% Dr. Kira Abercromby Code
% Double R Method
[dp1a_rvec, dp1a_vvec] = AERO557doubleR(dp1_rhoh(:,1),dp1_rhoh(:,2),dp1_rhoh(:,3),...
                    dp1_Rsite(:,1),dp1_Rsite(:,2),dp1_Rsite(:,3),dp1a_tau1,dp1a_tau3); % [km km/s]
[dp2a_rvec, dp2a_vvec] = AERO557doubleR(dp2_rhoh(:,1),dp2_rhoh(:,2),dp2_rhoh(:,3),...
                    dp2_Rsite(:,1),dp2_Rsite(:,2),dp2_Rsite(:,3),dp2a_tau1,dp2a_tau3); % [km km/s]

% use Izzo Gooding to find correct velocity at middle vectors
time = etime(dp2_time(2,:),dp1_time(2,:))/24/3600;
[dp_v1_rev0, dp_v2_rev0] = lambert(dp1a_rvec', dp2a_rvec', time, 0, 398600);
[dp_v1_rev1, dp_v2_rev1] = lambert(dp1a_rvec', dp2a_rvec', time, 1, 398600);

% backward propagate rev options to t1
options  = odeset('RelTol',1e-8,'AbsTol',1e-8);
dp_dtb   = [0 etime(dp1_time(1,:),dp2_time(2,:))];
[~,dp_rev0] = ode45(@orbitprop,dp_dtb,[dp2a_rvec' dp_v2_rev0],options,398600);
[~,dp_rev1] = ode45(@orbitprop,dp_dtb,[dp2a_rvec' dp_v2_rev1],options,398600);
dp_coes0 = sv2coes(dp_rev0(end,1:3),dp_rev0(end,4:6));
dp_coes1 = sv2coes(dp_rev1(end,1:3),dp_rev1(end,4:6));

% back propagate TLE's to t1
time1     = etime(dp1_time(1,:),double_passTLE1.epoch);
time2     = etime(dp1_time(1,:),double_passTLE2.epoch);
time3     = etime(dp1_time(1,:),double_passTLE3.epoch);
TLEfname1 = TLEpropagator('Double Unknown 1.txt',time1,'f',0);
TLEfname2 = TLEpropagator('Double Unknown 2.txt',time2,'f',0);
TLEfname3 = TLEpropagator('Double Unknown 3.txt',time3,'f',0);
tleinfo1  = readTLE(TLEfname1);
tleinfo2  = readTLE(TLEfname2);
tleinfo3  = readTLE(TLEfname3);
coes1     = coesfromTLE(tleinfo1);
coes2     = coesfromTLE(tleinfo2);
coes3     = coesfromTLE(tleinfo3);

% compare COEs from TLE specs and observations
rev0    = [rad2deg(dp_coes0.inc);rad2deg(dp_coes0.raan);dp_coes0.ecc;rad2deg(dp_coes0.omega);rad2deg(dp_coes0.Me);sqrt(398600/dp_coes0.a^3)*24*3600/(2*pi);dp_coes0.a];
rev1    = [rad2deg(dp_coes1.inc);rad2deg(dp_coes1.raan);dp_coes1.ecc;rad2deg(dp_coes1.omega);rad2deg(dp_coes1.Me);sqrt(398600/dp_coes1.a^3)*24*3600/(2*pi);dp_coes1.a];
TLE1    = [rad2deg(coes1.inc);rad2deg(coes1.raan);coes1.ecc;rad2deg(coes1.w);rad2deg(coes1.Me);coes1.N*24*3600/(2*pi);coes1.a];
TLE2    = [rad2deg(coes2.inc);rad2deg(coes2.raan);coes2.ecc;rad2deg(coes2.w);rad2deg(coes2.Me);coes2.N*24*3600/(2*pi);coes2.a];
TLE3    = [rad2deg(coes3.inc);rad2deg(coes3.raan);coes3.ecc;rad2deg(coes3.w);rad2deg(coes3.Me);coes3.N*24*3600/(2*pi);coes3.a];
col_n   = {'Obs w/ 0 Rev','Obs w/ 1 Rev','TLE 1','TLE 2','TLE 3'};
row_n   = {'inc  (deg)','raan (deg)','ecc','argp (deg)','me   (deg)','n    (rev/day)','a    (km)'};
dp_comp = table(rev0,rev1,TLE1,TLE2,TLE3,'VariableNames',col_n,'RowNames',row_n);
disp(dp_comp)

%% Part 4 Info
filename = TLEfname1;
dp_rvec  = dp_rev1(end,1:3);
dp_vvec  = dp_rev1(end,4:6);
period   = dp_coes1.T;

%% Part 4
% propagate r and v vectors
options     = odeset('RelTol',1e-8,'AbsTol',1e-8);
dp_dt       = [0 period];
mu          = 398600; % km^3/s^2
dp_state    = [dp_rvec dp_vvec];
[dp_t,dp_s] = ode45(@orbitprop,dp_dt,dp_state,options,mu);

% update vectors & get new coes
dp_rvec     = dp_s(end,1:3)';
dp_vvec     = dp_s(end,4:6)';
dp_coes     = sv2coes(dp_rvec,dp_vvec);

% propagate TLE
[dp_TLEfname, dp_tlet, dp_tles] = TLEpropagator(filename,period,'f',1);
dp_tleinfo  = readTLE(dp_TLEfname);
dp_tlecoes  = coesfromTLE(dp_tleinfo);
[r, v]      = coes2sv(dp_tlecoes.h, dp_tlecoes.ecc, dp_tlecoes.inc, dp_tlecoes.raan, dp_tlecoes.w, dp_tlecoes.TA, 398600);

figure
plot3(dp_tles(:,1),dp_tles(:,2),dp_tles(:,3))
hold on
plot3(dp_s(:,1),dp_s(:,2),dp_s(:,3))
plot3(dp_tles(end,1),dp_tles(end,2),dp_tles(end,3),'o','MarkerSize',8,'MarkerEdgeColor','#7E2F8E','MarkerFaceColor','#7E2F8E')
plot3(dp_s(end,1),dp_s(end,2),dp_s(end,3),'o','MarkerSize',8,'MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F')
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('True Orbit','Guess Orbit','True Position','Guess Position')
hold off

% calc RA/DEC off r vec
dp_rvec = dp_rvec/norm(dp_rvec);
lobs    = dp_rvec(1);
mobs    = dp_rvec(2);
nobs    = dp_rvec(3);
DECobs  = asind(nobs);
if mobs > 0
    RAobs = acosd(lobs/cosd(DECobs));
else
    RAobs = 360 - ascod(lobs/cosd(DECobs));
end

r       = r/norm(r);
ltle    = r(1);
mtle    = r(2);
ntle    = r(3);
DECtle  = asind(ntle);
if mtle > 0
    RAtle = acosd(ltle/cosd(DECtle));
else
    RAtle = 360 - ascod(ltle/cosd(DECtle));
end

figure
scatter(RAtle,DECtle,'filled')
hold on
scatter(RAobs,DECobs,'filled')
xlim([RAtle-1.25 RAtle+1.25])
ylim([DECtle-1.25 DECtle+1.25])
ylabel('Declination')
xlabel('Right Ascension')
legend('TLE','Obs')


