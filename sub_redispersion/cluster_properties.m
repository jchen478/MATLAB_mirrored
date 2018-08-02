clc;
clear;
close all;

% read file
filename = 'Cluster_resultsTotal.txt';
if (exist(filename, 'file') ~= 0)
    File = fopen(filename,'r');
    data = fscanf(File,'%f',[20 Inf])';
    fclose(File);
end

% conditions
muArr = caseArr('$\mu$',[0 5 10 15], 1);
attArr = caseArr('$A_N$',[0 9 20 30 50], 2);

% data arrays
nC = zeros(muArr.ndata,attArr.ndata);
S = zeros(muArr.ndata,attArr.ndata);
Rg = zeros(muArr.ndata,attArr.ndata);
lamroot1 = zeros(muArr.ndata,attArr.ndata);
lamroot2 = zeros(muArr.ndata,attArr.ndata);
lamroot3 = zeros(muArr.ndata,attArr.ndata);
I1 = zeros(muArr.ndata,attArr.ndata);
I2 = zeros(muArr.ndata,attArr.ndata);
I3 = zeros(muArr.ndata,attArr.ndata);
p11 = zeros(muArr.ndata,attArr.ndata);
p12 = zeros(muArr.ndata,attArr.ndata);
p13 = zeros(muArr.ndata,attArr.ndata);
p21 = zeros(muArr.ndata,attArr.ndata);
p22 = zeros(muArr.ndata,attArr.ndata);
p23 = zeros(muArr.ndata,attArr.ndata);
p31 = zeros(muArr.ndata,attArr.ndata);
p32 = zeros(muArr.ndata,attArr.ndata);
p33 = zeros(muArr.ndata,attArr.ndata);
cond1 = zeros(muArr.ndata,attArr.ndata);
cond2 = zeros(muArr.ndata,attArr.ndata);
xproj = zeros(muArr.ndata,attArr.ndata);
yproj = zeros(muArr.ndata,attArr.ndata);
zproj = zeros(muArr.ndata,attArr.ndata);

% place data into storage
for i=1:muArr.ndata
    for j=1:attArr.ndata
        ind = (i-1)*attArr.ndata + j;
        cond1(i,j) = data(ind,1);
        cond2(i,j) = data(ind,2);
        nC(i,j) = data(ind,3);
        S(i,j) = data(ind,4);
        Rg(i,j) = data(ind,5);
        lamroot1(i,j) = data(ind,6);
        lamroot2(i,j) = data(ind,7);
        lamroot3(i,j) = data(ind,8);
        I1(i,j) = data(ind,9);
        I2(i,j) = data(ind,10);
        I3(i,j) = data(ind,11);
        p11(i,j) = data(ind,12);
        p12(i,j) = data(ind,13);
        p13(i,j) = data(ind,14);
        p21(i,j) = data(ind,15);
        p22(i,j) = data(ind,16);
        p23(i,j) = data(ind,17);
        p31(i,j) = data(ind,18);
        p32(i,j) = data(ind,19);
        p33(i,j) = data(ind,20);
        
        vec = [p31(i,j),p32(i,j),p33(i,j)];
        vec = abs(vec);
        [pmax,I] = max(vec); 
        if (I == 1)
            zproj(i,j) = lamroot1(i,j)/pmax;
        elseif (I == 2)
            zproj(i,j) = lamroot2(i,j)/pmax;
        else
            zproj(i,j) = lamroot3(i,j)/pmax;
        end
        
        vec = [p21(i,j),p22(i,j),p23(i,j)];
        vec = abs(vec);
        [pmax,I] = max(vec); 
        if (I == 1)
            yproj(i,j) = lamroot1(i,j)/pmax;
        elseif (I == 2)
            yproj(i,j) = lamroot2(i,j)/pmax;
        else
            yproj(i,j) = lamroot3(i,j)/pmax;
        end
        
        vec = [p11(i,j),p12(i,j),p13(i,j)];
        vec = abs(vec);
        [pmax,I] = max(vec);
        if (I == 1)
            xproj(i,j) = lamroot1(i,j)/pmax;
        elseif (I == 2)
            xproj(i,j) = lamroot2(i,j)/pmax;
        else
            xproj(i,j) = lamroot3(i,j)/pmax;
        end
    end
end

%%
crit = (S>1280*0.9);

Rg = Rg.*crit;
Rg (Rg == 0) = nan;
S = S.*crit;
S (S == 0) = nan;
I1 = I1.*crit;
I1 (I1 == 0) = nan;
I2 = I2.*crit;
I2 (I2 == 0) = nan;
I3 = I3.*crit;
I3 (I3 == 0) = nan;
lamroot1 = lamroot1.*crit;
lamroot1 (lamroot1 == 0) = nan;
lamroot2 = lamroot2.*crit;
lamroot2 (lamroot2 == 0) = nan;
lamroot3 = lamroot3.*crit;
lamroot3 (lamroot3 == 0) = nan;
xproj = xproj.*crit;
xproj (xproj == 0) = nan;
yproj = yproj.*crit;
yproj (yproj == 0) = nan;
zproj = zproj.*crit;
zproj (zproj == 0) = nan;

RgObj = data2D('$R_g/L_x$',Rg/600);
I1Obj = data2D('$I_1$',I1);
I2Obj = data2D('$I_2$',I2);
I3Obj = data2D('$I_3$',I3);
lamroot1Obj = data2D('$\lambda_1$',lamroot1);
lamroot2Obj = data2D('$\lambda_2$',lamroot2);
lamroot3Obj = data2D('$\lambda_3$',lamroot3);
lamroot123Obj = data2D('$8\lambda_1\lambda_2\lambda_3 / V_{box}$',...
                8*lamroot1.*lamroot2.*lamroot3/600^3);
xprojObj = data2D('x',xproj);
zprojObj = data2D('z',zproj);
yprojObj = data2D('y',yproj);
yzprojObj = data2D('yz',4*yproj.*zproj/600^2);
xyzprojObj = data2D('xyz',8*xproj.*yproj.*zproj/600^3);
volfracObj = data2D('$\phi_{cluster} (v/v\%)$',100*S./(8*xproj.*yproj.*zproj)*150*pi);

plot2dim1(RgObj,attArr,muArr);
plot2dim1(I1Obj,attArr,muArr);
plot2dim1(I2Obj,attArr,muArr);
plot2dim1(I3Obj,attArr,muArr);
plot2dim1(lamroot1Obj,attArr,muArr);
plot2dim1(lamroot2Obj,attArr,muArr);
plot2dim1(lamroot3Obj,attArr,muArr);
plot2dim1(xprojObj,attArr,muArr);
plot2dim1(yprojObj,attArr,muArr);
plot2dim1(zprojObj,attArr,muArr);
plot2dim1(yzprojObj,attArr,muArr);
plot2dim1(xyzprojObj,attArr,muArr);
plot2dim1(lamroot123Obj,attArr,muArr);
plot2dim1(volfracObj,attArr,muArr);