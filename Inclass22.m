%Inclass 22

%1. Consider the case of the auto-activating gene that we discussed in class
%today. Make a bifurcation diagram for this system by varying the
%activated transcription rate for three cases - in which 1, 4, or 8 copies of the
%transcripton factor are necessary to activate transciption. Comment on your
%results. 

figure; hold on;
ku=0;
n=1;
gxfunc2 = @(x,kb) kb*x.^n./(1+x.^n)-x;
for kb=0:0.05:5
    gxfunc = @(x) gxfunc2(x,kb);
    for x0 = 0:0.1:3
        [rt,~,exitflag] = fzero(gxfunc,x0);
        if exitflag == 1
            plot(kb,rt,'k.');
        end
    end
end
hold off;
xlabel('k_b');ylabel('Fixed points');
set(gca,'Fontsize',24);

figure; hold on;
ku=0;
n=4;
gxfunc2 = @(x,kb) kb*x.^n./(1+x.^n)-x;
for kb=0:0.05:5
    gxfunc = @(x) gxfunc2(x,kb);
    for x0 = 0:0.1:3
        [rt,~,exitflag] = fzero(gxfunc,x0);
        if exitflag == 1
            plot(kb,rt,'k.');
        end
    end
end
hold off;
xlabel('k_b');ylabel('Fixed points');
set(gca,'Fontsize',24);

figure; hold on;
ku=0;
n=8;
gxfunc2 = @(x,kb) kb*x.^n./(1+x.^n)-x;
for kb=0:0.05:5
    gxfunc = @(x) gxfunc2(x,kb);
    for x0 = 0:0.1:3
        [rt,~,exitflag] = fzero(gxfunc,x0);
        if exitflag == 1
            plot(kb,rt,'k.');
        end
    end
end
hold off;
xlabel('k_b');ylabel('Fixed points');
set(gca,'Fontsize',24);

% 2. Make a similar diagram for the case of an autorepressing gene in the
% case that 2 copies are need to turn off the gene.


figure;hold on;
ku = 0;
for kb = 0:0.5:5
    polycoeff = [1 kb 1 -ku];
    rts = roots(polycoeff);
    rts = rts(imag(rts) == 0);
    plot (kb*ones(length(rts),1),rts,'r.');
   
end