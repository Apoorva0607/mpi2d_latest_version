function out_movie=inter_intime(in_movie)



 int_cine=zeros(size(in_movie,1),size(in_movie,2),245);
frames1=[10 12 13 14 16 17 20 22 23 29 31 32 41 43 45 46 48 52 54 57 58 60 66 68 69 72 75 78 79 81 84 87 88 91 94 97 100 101 103 104 107 110 111 113 116 119 122 123 126 129 132 133 135 136 138 139 141 144 145 147 148 150 151 153 154 156 157 160 163 166 169 172 175 176 178 179 181 182 185 188 191 194 195 197 198 200 201 204 207 210 213 216 219 222 223 225 226 229 232 235 238 241 244 247 248 250 251 254 ]
for i=1:size(in_movie,1)
    i
    for j=1:size(in_movie,2)
        
        curv=1.*squeeze(in_movie(i,j,:));
        curv_inter=interp1(frames1,curv,10:254,'cubic');
        int_cine(i,j,:)=curv_inter;
    end
end

out_movie=int_cine(:,:,1:end);
% save('1S_rest.mat','out_movie');

