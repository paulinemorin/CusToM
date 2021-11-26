function [Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if Foot_length > 0.2400 && Foot_length < 0.2786
    Foot_length_insole = 0.2486 ;
    Foot_large_insole = 0.0902;
elseif Foot_length > 0.2786 && Foot_length < 0.3042
    Foot_length_insole = 0.2742 ;
    Foot_large_insole = 0.0975 ;
else
    warning('Wrong foot or insole size')
end

end

