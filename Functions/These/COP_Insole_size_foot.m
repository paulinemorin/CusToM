function [CoP_size_feet] = COP_Insole_size_foot(filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
data_full=importdata([filename '.txt']);
liste=[22,23,44,45];

CoP_insole_row=[data_full.data(:,liste)]';

[Foot_length,Foot_large] = Foot_size(filename);
[Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large);
taille=[Foot_length_insole,Foot_large_insole];

%% Transformation cop data semelles

for ii=[1,2]
    CoP_size_feet((ii-1)*2+1,:)=[taille(1).*CoP_insole_row((ii-1)*2+1,:)];
    CoP_size_feet((ii-1)*2+2,:)=[taille(2).*CoP_insole_row((ii-1)*2+2,:)];
end
 

end
