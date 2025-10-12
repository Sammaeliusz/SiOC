blades = 6
x = 0:pi/100:4*pi;
frame = getframe(gcf).cdata;
for s = 120:10:200
    for l = 1:1:16
        o=0;
        for m = 2:1:size(frame, 1)
            f = sin((blades/2)*x+m*pi/s).^2;
            polarplot(x,f);
            hold on
            f = sin(-1*(blades/2-1)*x+m*pi/s*1.75+pi/3).^2;
            polarplot(x,f,'g');
            hold off
            for n = 1:1:l
                if o+n<size(frame, 1)
                    frame(o+n,:,:)=getframe(gcf).cdata(o+n,:,:);
                end
            end
            o=o+l;
        end
        %imwrite(frame,"35."+l+"."+s+".png");
    end
end