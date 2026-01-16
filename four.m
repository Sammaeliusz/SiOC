function out = anscombe(x)
    out = 2*sqrt(x+3/8);
end
function out = rev_anscombe(x)
    out = 0.25*x*x-0.125;
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
            next_gen = zeros(floor(size(curr_gen,1)/4^generation));
            generation = generation + 1;
            for j = 1:4:size(curr_gen,1)
                steps=steps+1;
                next_gen(steps)=(curr_gen(j)+curr_gen(j+1)+curr_gen(j+2)+curr_gen(j+3))/4;
                tmp_line(steps+size(curr_gen,1)/4)=(curr_gen(j)+curr_gen(j+1)-curr_gen(j+2)-curr_gen(j+3))/4;
                tmp_line(2*steps-1+size(curr_gen,1)/2)=(curr_gen(j)-curr_gen(j+1))/2;
                tmp_line(2*steps+size(curr_gen,1)/2)=(curr_gen(j+2)-curr_gen(j+3))/2;
            end
            curr_gen = next_gen;
        end
        tmp_line(1:size(next_gen,1))=next_gen;
        trans(:,count_von) = tmp_line;
    end
    trans = trans';
    end
    out = trans;
end
function out = MDMRFHT(v)
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
                curr_gen((4^generation)+1:4^(generation+1)) = i((4^generation)+1:min(4^(generation+1),size(i,1)));
            end
            generation = generation + 1;
        end
        trans(1:size(curr_gen,1),count_von) = curr_gen;
    end
    trans = trans';
    end
    out = trans;
end
function out =FWHT(v)
    if(size(v,2)==1)
        out = v;
    else
        half1 = v(1:round(size(v,2)/2));
        half2 = v(round(size(v,2)/2)+1:size(v,2));
        for i=1:size(half1,2)
            tmp1(i) = half1(i)+half2(i);
            tmp2(i) = half1(i)-half2(i);
        end
        out(1:round(size(v,2)/2)) = FWHT(tmp1);
        out(round(size(v,2)/2)+1:size(v,2)) = FWHT(tmp2);

    end
end
function out =RFWHT(v)
    i = 0;
    in = v
    while i < log2(size(in,2))
        for j = 1:2^(i+1):size(in,2)
            j+2^i
            for k = 1:i+1
                tmp(j+k-1)=(in(j+k-1)+in(j+k-1+2^i))/2;
                tmp(j+k-1+2^i)=(in(j+k-1)-in(j+k-1+2^i))/2;
            end
        end
        tmp
        in = tmp;
        i = i+1;
    end
end
function out = MDFWHT(v, rev)
    trans = v;
    for c = 1:2
    count_von =0;
    for i = trans
        count_von = count_von+1;
        if(rev)
            trans(:,count_von) = RFWHT(i);
        else
            trans(:,count_von) = FWHT(i);
        end
    end
    trans = trans'
    end
    out = trans;
end
MDMRFHT(MDMFHT([2,4,6,3,1,2,5,9,8,13,15,7,8,9,0,4;2,4,6,3,1,2,5,9,8,13,15,7,8,9,0,4;2,4,6,3,1,2,5,9,8,13,15,7,8,9,0,4;2,4,6,3,1,2,5,9,8,13,15,7,8,9,0,4]'));
RFWHT(FWHT([1,2,3,4,5,6,7,8]))