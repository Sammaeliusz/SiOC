function out = compare(x,y)
    out = zeros(size(x));
    if max(x)>1
    x=rescale(x);
    end
    if max(y)>1
    y=rescale(y);
    end
    for i = 1:min(size(x,1),size(y,1))
        for j = 1:min(size(x,2),size(y,2))
            out(i,j)=abs(x(i,j)-y(i,j));
        end
    end
    tsp = 0;
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(out(i,j)==0)
                tsp=tsp+1;
            end
        end
    end
    tsp
    tapa = tsp/(size(out,1)*size(out,2))*100
end
function out = compareC(x,y)
    out = zeros(size(x));
    for i=1:3
        out(:,:,i)=compare(x(:,:,i),y(:,:,i));
    end
end
function out = interpolate1D(vx, vy, probes, interpolationType)
out = zeros(size(probes)); % Initialize output array
if max(vy)>1
    vy=rescale(vy);
end
    if interpolationType=="Nearest"
        i=1;j=1;
        while i<=size(probes,2)
            while probes(i)>j+1
                j=j+1;
            end
            if i-probes(i)<i-probes(i)
                out(i) = vy(j);
            else
                out(i) = vy(j+1);
            end
            i=i+1;
        end
    end
    if interpolationType=="Linear"
        i = 1;j=1;
        while i <= size(probes, 2)
            while probes(i)>j+1
                j=j+1;
            end
            out(i) = vy(j)+(vy(j+1)-vy(j))/(j+1-j)*(probes(i)-j);
            i = i + 1;
        end
    end
    if interpolationType=="Cubical"
        i = 1;j=2;
        while i <= size(probes, 2)
            while probes(i)>j+1
                j=j+1;
            end
            t = (probes(i)-j);
            out(i) = (2*t^3-3*t^2+1)*vy(j)+(t^3-2*t^2+t)*(vy(j+1)-vy(j-1))/2+(-2*t^3+3*t^2)*vy(j+1)+(t^3-t^2)*(vy(j+2)-vy(j))/2;
            i = i + 1;
        end
    end
end
function out = interpolate2D (v,probes, interpolationType)
    i=1;
    while i<=size(v,1)
        tmp(i,:)=interpolate1D(size(v,2), v(i,:), probes, interpolationType);
        i=i+1;
    end
    i=1;
    while i<=size(tmp,2)
        out(i,:)=interpolate1D(size(tmp,1), tmp(:,i), probes, interpolationType);
        i=i+1;
    end
    out=out';
end
function out = resizeImageBW(v, newSize, interpolationType)
    i=1;
    tmp = zeros(size(v,1),floor(size(v,2)*newSize));
    while i<=size(v,1)
        tmp(i,:)=interpolate1D(size(v,2), v(i,:), linspace(1,size(v,2)-1,floor(size(v,2)*newSize)), interpolationType);
        i=i+1;
    end
    %out=tmp;
    i=1;
    while i<=size(tmp,2)
        out(i,:)=interpolate1D(size(tmp,1), tmp(:,i), linspace(1,size(tmp,1)-1,size(tmp,1)*newSize), interpolationType);
        i=i+1;
    end
    out=out';
end
function out = resizeImageC(v, newSize, interpolationType)
    out=cat(3,resizeImageBW(v(:,:,1), newSize, interpolationType),resizeImageBW(v(:,:,2), newSize, interpolationType),resizeImageBW(v(:,:,3), newSize, interpolationType));
end
function out = rotateImageBW(v, x, y, phi)
% 1,1->1,3 1,5->3,5 5,1->3,1 5,5 -> 3,5
    if max(v)>1
        v=rescale(v);
    end
    rot_mat=[cos(phi), sin(phi);-sin(phi), cos(phi)];
    out = zeros(size(v));
    dist = max([sqrt(x^2+y^2), sqrt((size(v,1)-x)^2+y^2), sqrt(x^2+(size(v,2)-y)^2), sqrt((size(v,1)-x)^2+(size(v,2)-y)^2)]);
    for i = 1:size(v,1)
        for j = 1:size(v,2)
            if(v(i,j)==0&&m~=(size(v,1)-1))
                new_pos = [i-x,j-y]*rot_mat;
                out(round(new_pos(1)+(x+dist)/sqrt(2)), round(new_pos(2)+(y+dist)/sqrt(2))) = v(i,j);
            end
        end
    end
    i=1;
    while (sum(out(:,i))==0)
        i=i+1;
    end
    out=out(:,i:size(out,2));
    i=1;
    while (sum(out(i,:))==0)
        i=i+1;
    end
    out=out(i:size(out,1),:);
end
function out = rotateImageC(v, x, y, phi)
    out=cat(3, rotateImageBW(v(:,:,1), x,y,phi), rotateImageBW(v(:,:,2), x,y,phi), rotateImageBW(v(:,:,3), x,y,phi));
end
function out = mosaic(v, type)
    if type == "Bayer"
      mo = cat(3, [1,0;0,0],[0,1;1,0],[0,0;0,1]);
    end
    if type == "Fujifilm"
      mo = cat(3, [0,0,0,0,1,0;1,0,1,0,0,0;0,0,0,0,1,0;0,1,0,0,0,0;0,0,0,1,0,1;0,1,0,0,0,0], [1,0,1,1,0,1;0,1,0,0,1,0;1,0,1,1,0,1;1,0,1,1,0,1;0,1,0,0,1,0;1,0,1,1,0,1],[0,1,0,0,0,0;0,0,0,1,0,1;0,1,0,0,0,0;0,0,0,0,1,0;1,0,1,0,0,0;0,0,0,0,1,0]);
    end
    mo=uint8(mo);
    for i = 1:size(mo,1):size(v,1)-size(mo,1)
        for j = 1:size(mo,2):size(v,2)-size(mo,2)
          v(i:i+size(mo,1)-1,j:j+size(mo,2)-1,:)=v(i:i+size(mo,1)-1,j:j+size(mo,2)-1,:).*mo;
        end
    end
    out = v;
end
function out = demosaic(v)
      max_len = 20;
      for i = 2:size(v,1)-1
          for j=1:size(v,2)
              if v(i,j,1)==0&&v(i-1,j,1)~=0
                  m=1;
                  while v(i+m,j,1)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,1)=v(i-1,j,1)/(m+1)*m+v(i+m,j,1)/(m+1);
              end
              if v(i,j,2)==0&&v(i-1,j,2)~=0
                  m=1;
                  while v(i+m,j,2)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,2)=v(i-1,j,2)/(m+1)*m+v(i+m,j,2)/(m+1);
              end
              if v(i,j,3)==0&&v(i-1,j,3)~=0
                  m=1;
                  while v(i+m,j,3)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,3)=v(i-1,j,3)/(m+1)*m+v(i+m,j,3)/(m+1);
              end
          end
      end
      for i = 1:size(v,1)
          for j=2:size(v,2)-1
              if v(i,j,1)==0&&v(i,j-1,1)~=0
                  m=1;
                  while v(i,j+m,1)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,1)=v(i,j-1,1)/(m+1)*m+v(i,j+m,1)/(m+1);
              end
              if v(i,j,2)==0&&v(i,j-1,2)~=0
                  m=1;
                  while v(i,j+m,2)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,2)=v(i,j-1,2)/(m+1)*m+v(i,j+m,2)/(m+1);
              end
              if v(i,j,3)==0&&v(i,j-1,3)~=0
                  m=1;
                  while v(i,j+m,3)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,3)=v(i,j-1,3)/(m+1)*m+v(i,j+m,3)/(m+1);
              end
          end
      end
      for i = 2:size(v,1)-1
          for j=1:size(v,2)
              if v(i,j,1)==0&&v(i-1,j,1)~=0
                  m=1;
                  while v(i+m,j,1)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,1)=v(i-1,j,1)/(m+1)*m+v(i+m,j,1)/(m+1);
              end
              if v(i,j,2)==0&&v(i-1,j,2)~=0
                  m=1;
                  while v(i+m,j,2)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,2)=v(i-1,j,2)/(m+1)*m+v(i+m,j,2)/(m+1);
              end
              if v(i,j,3)==0&&v(i-1,j,3)~=0
                  m=1;
                  while v(i+m,j,3)==0&&m<max_len
                    m=m+1;
                  end
                  v(i,j,3)=v(i-1,j,3)/(m+1)*m+v(i+m,j,3)/(m+1);
              end
          end
      end
  out = v;
end
tmp = [0;0.12;0.25;0.5;0.75;1];
tmp2 = [0.25;0.5;0.75;1;0;0.12];
tmp3 = [0.5;0.75;1;0;0.12;0.25];

o = cat(3, [tmp,tmp,tmp,tmp,tmp,tmp],[tmp2,tmp2,tmp2,tmp2,tmp2,tmp2], [tmp3,tmp3,tmp3,tmp3,tmp3,tmp3]);

%oo = interpolate2D(o(:,:,2), )
%o = rand(10,20);
tic
o = imread("demo4demo.png");
oob = mosaic(o, "Bayer");
ooob = demosaic(oob);
toc
tic
o = imread("demo4demo.png");
oo = mosaic(o, "Fujifilm");
ooo = demosaic(oo);
toc
imshow(o)
figure;
imshow(oo)
figure
imshow(ooo);
figure
imshow(compareC(o,ooo))
figure;
imshow(oob)
figure
imshow(ooob);
figure
imshow(compareC(o,ooob))
