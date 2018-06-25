%mirror
imax = 2304;
tmax = 6280;

%firehose
imax = 4096;
tmax = 3333;

% for mirror secular phase
js=1000;
je=5500;

% for mirror linear phase
js=1;
je=499;

% for mirror saturated phase
js=5500;
je=6279;

% for firehose saturated phase: 0.05
js=1333;
je=3333;

% for firehose secular phase: 0.013
js=266;
je=400;

% for firehose transition phase: 0.038
js=500;
je=1333;

js=8800;
je=9990;

clear b;

n=1;
for i=1:imax
j = js;

muss = smooth(mus(:,i),30);

while (j < je) 

%while (vprls(j+1,i)*vprls(j,i) > 0 && j < je-1)
%    j=j+1;
%end
jstart = j;
muo = muss(j);
while (muss(j)/muo < exp(1) && muss(j)/muo > exp(-1) && j < tmax)
    j=j+1;
end

b(n) = j-jstart;
n=n+1;
j=j+1;

end

end

[nem,cen] = hist(b,100);
plot(cen,log(nem),'.')


js=1;je=3333;
for i=1000:2304
    i
    
subplot(3,1,1); plot(t(js:je),vprls(js:je,i)); subplot(3,1,2); plot(t(js:je),mus(js:je,i)); subplot(3,1,3); plot(t(js:je),B(js:je,i));
pause(0.5);
end


clear b;

n=1;
for i=1:imax
j = js;

muss = smooth(mus(:,i),3);

while (j < je-10) 

jstart = j;
muo = muss(j);
while (muss(j)/muo < exp(1) && muss(j)/muo > exp(-1) && j < tmax ...
        && vprls(j+1,i)*vprls(j,i) > 0)
    j=j+1;
end

b(n) = j-jstart;
n=n+1;
j=j+1;

end

end

[nem,cen] = hist(b,100);
