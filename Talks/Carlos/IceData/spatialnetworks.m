%RGG %Quick code Carlos J Melian
%November 2013
J = 1000;r = 0.1;%r = unifrnd(0.01,1);
D = zeros(J,J);

%Asymptotic behavior
mu = J*(e^(-pi * r^2 * J))
MA = log(J) - log(mu);
MB = pi*J;
rc = sqrt(MA/MB);


n = unifrnd(0,1,J,2);
for i = 0:J-1;
     for j = i+1:J;
         A = (n(i,1) - n(j,1))^2;%Euclidean distance
         B = (n(i,2) - n(j,2))^2;
         d(i,j) = sqrt(A + B);
         if d(i,j) < r;
            D(i,j) = 1;
         else
            D(i,j) = 0;
         end
     end
end
D1=D+D';

%plot network
gplot(D1,n, "r.-")
set (get (gca, ("children")), "markersize", 12);

%giant component
[blocks,dag] = components(D1);AT = sort(blocks);
connectivity = [ find(AT(1:end-1) ~= AT(2:end)) length(AT) ];
numberclusters = AT(connectivity);
sizeclusters = diff([0 connectivity]);
max(sizeclusters)
