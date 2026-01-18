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
function out = anscombe(x)
    out = 2*sqrt(x+3/8);
end
function out = inv_anscombe(x)
    out = 0.25*x*x-0.125;
end
function out = MDAnscombe(v, inv)
    out = zeros(size(v));
    for i = 1:size(v,1)
        for j = 1:size(v,2)
            if(inv)
                out(i,j) = inv_anscombe(v(i,j));
            else
                out(i,j) = anscombe(v(i,j));
            end
        end
    end
end
function out = MDMFHT(v)
    trans = v;
    for c = 1:2
    count_von =0;
    for i = trans
        count_von = count_von+1;
        tmp_line = zeros(size(i));
        generation = 1;
        curr_gen = i;
        while size(curr_gen,1)>=4
            steps = 0;
            next_gen = zeros(size(curr_gen,1)/4, 1);
            generation = generation + 1;
            for j = 1:4:size(curr_gen,1)-3
                steps=steps+1;
                next_gen(steps)=(curr_gen(j)+curr_gen(j+1)+curr_gen(j+2)+curr_gen(j+3))/4;
                tmp_line(steps+round(size(curr_gen,1)/4))=(curr_gen(j)+curr_gen(j+1)-curr_gen(j+2)-curr_gen(j+3))/4;
                tmp_line(2*steps-1+round(size(curr_gen,1)/2))=(curr_gen(j)-curr_gen(j+1))/2;
                tmp_line(2*steps+round(size(curr_gen,1)/2))=(curr_gen(j+2)-curr_gen(j+3))/2;
            end
            curr_gen = next_gen;
        end
        tmp_line(1:length(next_gen))=next_gen;
        trans(:,count_von) = tmp_line;
    end
    trans = trans';
    end
    out = trans;
end
function out = MDMIFHT(v)
    %a =x+y+z b=x+y-z c=x-y+q d=x-y-q
    % x = steps y=steps+size(curr_gen,1)/4 z=2*steps-1+size(curr_gen,1)/2
    % q=2*steps+size(curr_gen,1)/2
    trans = v;
    for c = 1:2
    count_von=0;
    for i = trans
        count_von = count_von+1;
        generation = 1;
        curr_gen = i(1:4);
        next_gen = [];
        while generation <= log2(size(i,1))/2
            steps = 0;
            for j = 1:min(4^(generation-1), size(i,1)-3)
                steps=steps+1;
                next_gen(4*(j-1)+1)=curr_gen(steps)+curr_gen(steps+round(size(curr_gen,1)/4))+curr_gen(2*steps-1+round(size(curr_gen,1)/2));
                next_gen(4*(j-1)+2)=curr_gen(steps)+curr_gen(steps+round(size(curr_gen,1)/4))-curr_gen(2*steps-1+round(size(curr_gen,1)/2));
                next_gen(4*(j-1)+3)=curr_gen(steps)-curr_gen(steps+round(size(curr_gen,1)/4))+curr_gen(2*steps+round(size(curr_gen,1)/2));
                next_gen(4*(j-1)+4)=curr_gen(steps)-curr_gen(steps+round(size(curr_gen,1)/4))-curr_gen(2*steps+round(size(curr_gen,1)/2));
            end
            curr_gen = zeros(4^generation,1);
            curr_gen(1:4^generation) = next_gen;
            if(size(next_gen,2)<size(i,1))
                curr_gen((4^generation)+1:min(4^(generation+1),size(i,1))) = i((4^generation)+1:min(4^(generation+1),size(i,1)));
            end
            generation = generation + 1;
        end
        trans(1:size(curr_gen,1),count_von) = curr_gen;
    end
    trans = trans';
    end
    out = trans;
end
function out = MDFHT(v)
     trans = v;
    for c = 1:2
    for k = 1:size(trans,2)
        i = trans(:,k);
        [trans(:,length(i)),aaa] = wavedec(i, wmaxlev(length(i),"haar"), "haar");
    end
    trans = trans.';
    end
    out = trans;
end
function out = MFWHT(v)
    if(isscalar(v))
        out = v;
    else
        half1 = v(1:round(length(v)/2));
        half2 = v(round(length(v)/2)+1:length(v));
        for i=1:length(half1)
            tmp1(i) = half1(i)+half2(i);
            tmp2(i) = half1(i)-half2(i);
        end
        out(1:round(length(v)/2)) = MFWHT(tmp1);
        out(round(length(v)/2)+1:length(v)) = MFWHT(tmp2);

    end
end
function out = MIFWHT(v)
    i = 0;
    in = v;
    tmp=zeros(size(i,2));
    while i < log2(size(in,2))
        for j = 1:2^(i+1):size(in,2)
            for k = 1:i+1
                tmp(j+k-1)=(in(j+k-1)+in(j+k-1+2^i))/2;
                tmp(j+k-1+2^i)=(in(j+k-1)-in(j+k-1+2^i))/2;
            end
        end
        in = tmp;
        i = i+1;
    end
    out = tmp;
end
function out = MDFWHT(v, inv)
    trans = v;
    if mod(size(trans,1),2)==1
        trans = trans(1:size(trans,2),:);
    end
    if mod(size(trans,2),2)==1
        trans = trans(:,1:size(trans,1));
    end
    for c = 1:2
    count_von =0;
    for i = trans
        count_von = count_von+1;
        if(inv)
            trans(:,count_von) = MIFWHT(i);
        else
            trans(:,count_von) = MFWHT(i);
        end
    end
    trans = trans'
    end
    out = trans;
end
function out = DCT(v)
    out = zeros(size(v));
    for i = 1:length(v)
        for j = 1:length(v)
            out(i)=out(i)+(v(j)*cos(pi*(j-1+0.5)*(i-1)/length(v)));
        end
    end
    out = out * sqrt(2/length(v));
    out(1)=out(1)/sqrt(2);
end
function out = IDCT(v)
    out = zeros(size(v));
    for i = 1:length(v)
        out(i)=sqrt(0.5)*v(1);
        for j = 2:length(v)
            out(i)=out(i)+(v(j)*cos(pi*(i-1+0.5)*(j-1)/length(v)));
        end
    end
    out = out * sqrt(2/length(v));
end
function out = MDDCT(v, inv)
    trans = v;
    for c = 1:2
    for k = 1:size(trans,2)
        i = trans(:,k);
        if(inv)
            trans(:,k) = IDCT(i);
        else
            trans(:,k) = DCT(i);
        end
    end
    trans = trans.';
    end
    out = trans;
end
function out = Threshold(v, T)
   out=zeros(size(v));
   for i = 1:size(v,1)
       for j = 1:size(v,2)
            if abs(v(i,j))>T
                out(i,j)=v(i,j);
            else
                out(i,j)=0;
            end
       end
   end
end
function out = Quantization(v, Q)
   out=zeros(size(v));
   for i = 1:size(v,1)
       for j = 1:size(v,2)
            out(i,j)=floor(Q*v(i,j)+0.5)/Q;
       end
   end
end
function out = wygladzanie(v, transfortmat, quant, TQ)
for i = 1:3
    tmp = MDAnscombe(v(:,:,i), false);
    if transfortmat == "DCT"
        tmp = MDDCT(tmp,false);
    elseif transfortmat == "WHT"
        %tmp = MDFWHT(tmp, false);
        tmp = fwht(tmp);
    elseif transfortmat == "MFHT"
        if mod(size(tmp,1),2)==1
            tmp = tmp(1:size(tmp,1)-1,:);
        end
        if mod(size(tmp,2),2)==1
            tmp = tmp(:,1:size(tmp,2)-1);
        end
        [tmp,h,vv,d] = haart2(tmp);
        
        %tmp = MDMFHT(tmp);
    elseif transfortmat == "53"
        [tmp,s] = wavedec2(tmp,wmaxlev(size(tmp),"bior2.2"),"bior2.2");
    else
        [tmp,s] = wavedec2(tmp,wmaxlev(size(tmp),"bior4.4"),"bior4.4");
    end
    if quant
        tmp = Quantization(tmp, TQ);
    else
        tmp = Threshold(tmp, TQ);
    end
    if transfortmat == "DCT"
        tmp = MDDCT(tmp,true);
    elseif transfortmat == "WHT"
        %tmp = MDFWHT(tmp, true);
        tmp = ifwht(tmp);
    elseif transfortmat == "MFHT"
        tmp = ihaart2(tmp,h,vv,d);
        %tmp = MDMIFHT(tmp);
    elseif transfortmat == "53"
        tmp = waverec2(tmp,s,"bior2.2");;
    else
        tmp = waverec2(tmp,s,"bior4.4");;
    end
    out(:,:,i) = MDAnscombe(tmp, true);
end
if(size(out,1)~=size(v,1)||size(out,2)~=size(v,2))
    out=out(1:size(v,1)-1,1:size(v,2)-1,:);
end
%out = rescale(out);
end
tst = IDCT(DCT([1,2,3,4,5,6,7,8]))
org = imread("Motur.png");
o = poissrnd(double(org));
o = rescale(o);
obraz2 = wygladzanie(obraz, "MFHT", false, 0.1);
figure
tiledlayout(2,5)
nexttile
imshow(o)
for k=1:10
    k
oo = wygladzanie(o, "WHT", false,k/500);
%mse(i,k)=sum(sum(sum((compare(org,oo)).^2)))/(size(oo,1)*size(oo,2));
%mae(i,k)=sum(sum(sum(abs(compare(org,oo)))))/(size(oo,1)*size(oo,2));
nexttile
imshow(oo);
end
mse
mae