function output_img=poly_fit(input_img)

% fitting to a third order polynomial adapted from http://www.nmr.mgh.harvard.edu/~fhlin/codes/sense_recon_easy.m

mask1 = input_img > 0.1.*max(input_img(:));

[y,x]=find(mask1==1);
L = find(mask1==1); 

map01 = input_img(L);
A = [x.^3 y.^3 x.^2 x.*y y.^2 x y ones(size(x))];
b = map01;
coeffs = A\b;


ones_img= ones(size(input_img));
LL = find(ones_img == 1);
[yy,xx] = find(ones_img == 1);
AA = [xx.^3 yy.^3 xx.^2 xx.*yy yy.^2 xx yy ones(size(xx))];
cc = AA * coeffs;

output_img=input_img*0;
output_img(LL)=cc;

return;


