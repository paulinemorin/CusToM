function [time_g,time_d] = TimingFootFlat(Contact_detection)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
delay = 35;
%delay = 40;
[~,time_g]=max(Contact_detection(1,delay:end-delay));
[~,time_d]=max(Contact_detection(2,delay:end-delay));

time_g=time_g+delay;
time_d=time_d+delay;
end

