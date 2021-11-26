function Pressure_200 = double_freq(Pressure_data)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

NbPointsPrediction = size(Pressure_data, 1);
nbframe_pressure = size(Pressure_data, 2);

Pressure_200=zeros(NbPointsPrediction, 2*nbframe_pressure);

for i = 1 : nbframe_pressure-1
    Pressure_200(:,2*i-1) = Pressure_data(:,i);
    Pressure_200(:,2*i) = (Pressure_data(:,i)+Pressure_data(:,i+1))/2;
end
Pressure_200(:,2*nbframe_pressure)=Pressure_data(:,nbframe_pressure);

end

