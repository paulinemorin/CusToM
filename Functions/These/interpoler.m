function Pressure_200 = interpoler(Pressure_data)
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


CoP_200=zeros(4, 2*nbframe_pressure);

for i = 1 : nbframe_pressure-1
    CoP_200(:,2*i-1) = CoP(:,i);
    CoP_200(:,2*i) = (CoP(:,i)+CoP(:,i+1))/2;

    if  or (CoP(1:2,i)==[0;0] , CoP(1:2,i+1)==[0;0])
        CoP_200(1:2,2*i)=[0;0];
    end
    if  or (CoP(3:4,i)==[0;0] , CoP(3:4,i+1)==[0;0])
        CoP_200(3:4,2*i)=[0;0];
    end
end
CoP_200(:,2*nbframe_pressure)=CoP(:,nbframe_pressure);

for i = 1 : nbframe_pressure*2
    if CoP_200(1,i)==0 && CoP_200(2,i)==0
        Pressure_200(1:16,i)=zeros(16,1);
    end
    if CoP_200(3,i)==0 && CoP_200(4,i)==0
        Pressure_200(17:32,i)=zeros(16,1);
    end
end

end

