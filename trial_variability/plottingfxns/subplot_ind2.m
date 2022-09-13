function [sp1,sp2] = subplot_ind2(Ns)

sp1 = ceil(sqrt(Ns));
while mod(Ns,sp1) ~= 0
    sp1 = sp1 - 1;    
end
sp2 = Ns / sp1;

if (sp1 == 1 || sp2 == 1) & Ns > 5
    sp1 = ceil(sqrt(Ns));
    while mod(Ns,sp1) ~= 0
        Ns = Ns + 1;
    end
    sp2 = Ns / sp1;
end