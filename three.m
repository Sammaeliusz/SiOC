zak = linspace(-1, 1, 1000);
function out = probe2D(k, fy)
    zak = linspace(-1, 1, 1000);
    tmp = interp1(zak, fy, linspace(-1, 1, k+2));
    tmp2 = (tmp'*tmp);
    out = tmp2(2:size(tmp2,1)-1,2:size(tmp2,1)-1)/((size(tmp2,1)-2)*2);
end
function out = filter(v, filtr)
if max(v)>1
    v=rescale(v);
end
    out = zeros(size(v));
    for i = 1:size(v,1)-size(filtr,1)
        for j = 1:size(v,2)-size(filtr,2)
          %if sum(sum(out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)))>=75
          %    out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)=(out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)+(v(i,j,1)*filtr))/2;
          %    out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)=(out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)+(v(i,j,2)*filtr))/2;
          %    out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)=(out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)+(v(i,j,3)*filtr))/2;
          %else
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)= out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)+(v(i,j,1)*filtr);
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)=out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)+(v(i,j,2)*filtr);
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)=out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)+(v(i,j,3)*filtr);
          %end
        end
    end
    if max(out)>1
        out=rescale(out);
    end
end
o = double(imread("Karnegia olbrzymia (saguaro) (2).jpg"));
a = poissrnd(4, size(o));
o = rescale((o/10 .* a)*10);
imshow(o);
filtr = probe2D(17, triangularPulse(zak));
oo = filter(o, filtr);
figure
imshow(oo);