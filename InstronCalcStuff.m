% Instron calc
% Surface area:
d = 33.45
h = 9.71
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Instron(:,2)/SA;
STRAIN = Instron(:,1)/h;

figure; plot(STRESS/STRAIN)

firstIndex = 70
secondIndex = 40
(STRESS(firstIndex)-STRESS(secondIndex)) / (STRAIN(firstIndex)-STRAIN(secondIndex))

% Units: N/m^2

STRESSm = Instron(:,2)/SAm
STRAINm = Instron(:,1)/hm;

(STRESSm(firstIndex)-STRESSm(secondIndex)) / (STRAINm(firstIndex)-STRAINm(secondIndex))

% Strain pct
figure; plot(STRAIN,STRESS)

% Sliding window slopes
windowLength = 10;
for k = windowLength:length(STRAIN)
    slopeInstron(k-(windowLength-1)) = (STRESSm(k)-STRESSm(k-(windowLength-1))) / (STRAINm(k)-STRAINm(k-(windowLength-1)))
end

figure; plot(slopeInstron)

