% from pg. 7 Huff and Feller, "Relation of circulating red cell volume
% to body density and obesity", Journal of Clinical Investigation, 35, 1:7 (1956)
%
function [BV,bfp,bmi] = bloodVolume(height,weight,age,gender,race)

if (nargin < 5) race = 'caucasian'; end;

[bfp,bmi] = bodyFatPercentage(height,weight,age,gender,race);

lean_body_weight = (1-bfp/100)*weight;
fat_weight = (bfp/100)*weight;

if (strcmp(gender,'male'))
    BV = 0.042*lean_body_weight + 0.025*fat_weight + 1.92;
else
    BV = 0.045*lean_body_weight + 0.034*fat_weight + 1.02;
end;

return;