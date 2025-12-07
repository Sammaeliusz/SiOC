zak = linspace(-1, 1, 1000);
function out = lambda(x)
    if(abs(x)<=1)
    out = 1-abs(x);
    else
        out = 0;
    end
end
function out = pi(x)
    out = zeros(size(x));
    for e = 1:size(x,2)
    if(abs(x(e))<0.5)
        out(e) = 1;
    else 
        out(e) = 0;
    end
    end
end
function out = omega(x)
    ax = abs(x);
    if(ax<=1)
        out = 1.4.*ax.^3-2.4.*ax.^2+1;
    elseif(ax<2)
        out = -0.6.*ax.^3+3*ax.^2-4.8.*ax+2.4;
    else
        out = 0;
    end
end

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
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)=out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,1)+(v(i,j,1)*filtr);
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)=out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,2)+(v(i,j,2)*filtr);
              out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)=out(i:i+size(filtr,1)-1,j:j+size(filtr,2)-1,3)+(v(i,j,3)*filtr);
        end
    end
    if max(max(out))>=1
        out=out./max(max(out));
    end
end
o = rescale(double(imread("Karnegia olbrzymia (saguaro) (2).jpg")));
function out = task(o, funk, it,k)
eta = 64;
o = poissrnd(o/eta)*eta;
filtr = probe2D(3*2^k, funk);
oo = filter(o, filtr);
ooo = imresize(oo, 1/8,it);
out = imresize(ooo, 8, it);
end
%{
for k=1:5
figure
tiledlayout(2,2)
imshow(o)
nexttile
oo = task(o, pi(zak), "nearest",k);
mse(1,k)=sum(sum(sum((o-oo).^2)))/1048576;
mae(1,k)=sum(sum(sum(abs(o-oo))))/1048576;
imshow(oo);
nexttile
oo=task(o, lambda(zak), "bilinear",k);
imshow(oo);
mse(2,k)=sum(sum(sum((o-oo).^2)))/1048576;
mae(2,k)=sum(sum(sum(abs(o-oo))))/1048576;
nexttile
oo=task(o, omega(zak), "bicubic",k);
imshow(oo);
mse(3,k)=sum(sum(sum((o-oo).^2)))/1048576;
mae(3,k)=sum(sum(sum(abs(o-oo))))/1048576;
end
mse
mae
%}
eta = 64;
o = poissrnd(o/eta)*eta;
imshow(o)