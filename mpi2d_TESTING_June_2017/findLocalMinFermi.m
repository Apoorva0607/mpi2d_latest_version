%
%
%  Brian Martin - 3/18/13
%
%  This code will find the spot to cut off the data for the Fermi model in
%  mpi2d.  Give it the data (usually deltaSICurve(1,:)) and it should
%  find the minimum right after the first initial circulation peak or allow
%  you to find it manually.

function index = findLocalMinFermi(data)
frames = 1:length(data);

% find where the initial circ peak begins (ish)
for i = 1:length(data)
    if( data(i+5) - data(i) > 150)
        index = i+5;
        break;
    end
    
end

% find the minimum after the initial circ peak
start = index;
good_indexes_found(1) = start;
figure(5833923);
plot(data);
hold on;
plot(frames(index),data(index),'dm');
legend('red x''s are each NEW min found','magenta diamond is where search for min begins');

spot = 2;
for ii = start:length(data) - 4
    
    % Change index to the new min and plot each new minimum found
    if (data(ii) < data(index))
        index = ii;
        good_indexes_found(spot) = index;
        spot = spot + 1;
        plot(frames(index),data(index),'xr');
        pause(.1)
    end
    
end
hold off;

% This will return a boolean so that a manual mode for picking the
% correct cut off point will be allowed if something strange is happening
% in the AIF, like a sudden dip in signal intensity at the end or something
% similar.

% The try/catch block is mainly in case the first index found is the
% minimum, which shouldn't be the case; but if it is, then will allow the
% user to pick the minimum manually/visually.

% Can probably make this cleaner, should do at a later time.
try
if(good_indexes_found(length(good_indexes_found)) - good_indexes_found(length(good_indexes_found) - 2) > 10)
    index = false;
    
    fprintf(2,'\nIt looks like something weird is happening, please pick the cut off manually\n');
    fprintf('A screen will appear with the AIF and some red x''s signaling the minimums \n(not the plot you''re looking at right now)\n');
    disp('Pick the x that is at a minimum after the first circ peak');
    disp('Example: If it is the second x, type in 2 and press enter (return) ');
    disp('Type return to have that next plot appear');
    keyboard;
    
end
catch
    index = false;
    % the '2' is to make the text of that line red, that is all it does.
    fprintf(2,'\nIt looks like something weird is happening, please pick the cut off manually\n');
    fprintf('A screen will appear with the AIF and some red x''s signaling the minimums \n(not the plot you''re looking at right now)\n');
    disp('Pick the x that is at a minimum after the first circ peak');
    disp('Example: If it is the second x, type in 2 and press enter (return) ');
    disp('Type return to have that next plot appear');
    keyboard;
end

end