function mat2txt(temp,mat,framesToSelect,ranget)

%   counter=1;
%  if ~isempty(framesToSelect)
%           for j=framesToSelect    % EVRD 8/17/11, to select ranget. 
%             if (j <=max(ranget)) && (j >1)
%                  actualFrameNumbers(counter)=j;
%                  counter=counter+1;
%             end
%           end
% else
%              for j=ranget
%                 actualFrameNumbers(counter)=j;
%                 counter=counter+1;
%              end
%  end
%        
% size(mat)
fid = fopen(temp,'wt');
counter=1;

for ii = 1:size(mat,1)
    fprintf(fid,'%f\t',mat(counter,1));
    fprintf(fid,'%f\t',mat(counter,2));
% 	fprintf(fid,'%d',actualFrameNumbers( counter));
    fprintf(fid,'\n');
	counter=counter+1;
end
fclose(fid)