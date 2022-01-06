clearvars;
CSVFiles = natsortfiles(rdir('E:\GelatinHalfHalf20211221\gelation_phantom_RIP_20211221.is_comp_RawData\*.csv'))
for k = 13:16
opts = delimitedTextImportOptions("NumVariables", 3);
    
    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Var1", "Extension", "Load"];
    opts.SelectedVariableNames = ["Extension", "Load"];
    opts.VariableTypes = ["string", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    
    % Import the data
    SpecimenRawData1 = readtable(CSVFiles(k).name, opts);
    
    Extension(k-12,:) = table2array(SpecimenRawData1(:,1));
    Load(k-12,:) = table2array(SpecimenRawData1(:,2));
end

figure;
for k = 1:size(Extension,1)
    plot(Extension(k,:),Load(k,:));
    hold on;
end
hold off; 

Extension = mean(Extension);
Load = mean(Load);
figure; plot(Extension,Load)
d = 32.51
h = 7.99
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Load/SA;
STRAIN = Extension/h;

figure; plot(STRAIN,STRESS)

STRAIN = STRAIN(62:end);
STRESS = STRESS(62:end);
figure; plot(STRAIN,STRESS)

firstIndex = 70
secondIndex = 40
(STRESS(firstIndex)-STRESS(secondIndex)) / (STRAIN(firstIndex)-STRAIN(secondIndex))

% Units: N/m^2

STRESSm = Load/SAm
STRAINm = Extension/hm;

(STRESSm(firstIndex)-STRESSm(secondIndex)) / (STRAINm(firstIndex)-STRAINm(secondIndex))

% Strain pct
figure; plot(STRAIN,STRESS)

% Sliding window slopes
windowLength = 10;
for k = windowLength:length(STRAIN)
    slopeInstron(k-(windowLength-1)) = (STRESSm(k)-STRESSm(k-(windowLength-1))) / (STRAINm(k)-STRAINm(k-(windowLength-1)))
end

figure; plot(slopeInstron)

STRESSm(/STRAINm

