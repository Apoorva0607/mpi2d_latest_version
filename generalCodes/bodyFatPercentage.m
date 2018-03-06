% computation based on Eq. 4 in Gallagher, Heymsfield, Heo, Jebb,
% Murgatroyd, and Sakamoto, Am. J. Clin. Nutr. 72 694:701 (2000).
%
function [bfp,bmi] = bodyFatPercentage(height,weight,age,gender,race)

if (nargin < 5) race = 'caucasian'; end;

if (strcmp(gender,'male')) sex = 1; else sex = 0; end;
if (strcmp(race,'asian')) asian = 1; else asian = 0; end;
if (strcmp(race,'african-american')) african_american = 1; else african_american = 0; end;

bmi = weight/height^2;

bfp = 63.7 - 864/bmi - 12.1*sex + 0.12*age + 129*asian/bmi - 0.091*asian*age - 0.030*african_american*age;

return;

