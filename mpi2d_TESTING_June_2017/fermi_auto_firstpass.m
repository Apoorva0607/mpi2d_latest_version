function firstpass_cutoff=fermi_auto_firstpass(curve,series,slice)

h = fspecial('gaussian', [1 size(curve,2)], 4);
filt_aif = imfilter(curve,h,'replicate');

diff_aif=diff(filt_aif);
step_var=diff_aif;
step_var(find(step_var<0))=-1;
step_var(find(step_var>=0))=1;

min_diff=find(diff_aif==min(diff_aif));
tmp1=find(step_var>0);

try
tmp2=tmp1(find(tmp1>=min_diff));
firstpass_cutoff=tmp2(1)-1;
catch
firstpass_cutoff=length(curve);    
end
% subplot(2,1,1)
% plot(curve);
% xlim([0 size(curve,2)])
% subplot(2,1,2)
% plot(curve(1:firstpass_cutoff));
% xlim([0 size(curve,2)])

 save(['firstpass_cutoff_',int2str(series),'_',int2str(slice),'.mat'],'firstpass_cutoff')