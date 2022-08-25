function [TargetMineralVolumes,mineralWF,mineralVF, ElementWF, TargetMineralWeights, Rhob] = CreateSyntheticData() 
    %use randfixedsum from MATLAB central file exchange

    %Create a table of the elemental composition of pure minerals
    [mineralTable, uncertainties, Elements] = CreateMineralTable();

    %Create 20,000 samples of random compositions
    [mineralWF, ElementWF] = createCompositionRand(mineralTable,  1:23, 20000, uncertainties);

    %Create a table of known densities of subject minerals
    [densityTable]=CreateMineralDensitiesTable();

    %Create volume fractions
    [mineralVF,Rhob]=CreatemineralVF(densityTable,mineralWF);

    %Sum minerals weight fractions to groups of similar compositions
    [TargetMineralWeights, MineralNames]= createTargetWeights(mineralWF);

    %Sum minerals volume fractions to groups of similar compositions
    [TargetMineralVolumes]=createTargetVolumes(mineralVF);
 
end

function [mineralTable, uncertainties, Elements] = CreateMineralTable()

        % Elements =    [Al,	    Ca,	    Fe,	    Gd,	     Mg,	     S,	    Si,	    Ti,         K]
        datat=[	        0,	        40.04,	0,	    0,	     0,	         0,	    0,	    0,          0;	...      % Calcite 1
                        0,          21.73,  0,	    0,	     13.18,	     0,	    0,	    0,          0;	...      % Dolomite 2
                        0,          29.44,	0,      0,       0,          23.56, 0,      0,          0;	...      % Anhydrite 3
                        0,          0,      0,      0,       0,          0,     46.75,	0,          0;	...      % Quartz 4
                        9.69,       0,      0,      0,       0,          0,     30.27,  0,          14.05;	...  % Kspar/Orthoclase 5
                        10.29,      0,      0,      0,       0,          0,     32.13,  0,          0; ...       % Plagioclase 6
                        10.77,      0.76,   0,      0,       0,          0,     31.50,	0,          0;	...      % Plagioclase/Albite 7
                        12.2,       3.02,	0,      0,       0,          0,     29.63,	0,          0;	...      % Plagioclase/Oligoclase 8
                        19.4,       14.41,	0,      0,       0,          0,     20.19,	0,          0;	...      % Plagioclase/Anorthite 9
                        0,          0,      48.2,	0,       0,          0,     0,      0,          0;	...      % Siderite 10               
                        20.9,       0,  	0,      0,       0,          0,     21.76,	0,          0;	...      % Kaolinite 11
                        9.07,   	0,      11.73,	0,       15.31,      0,     14.16,	0,          0;	...      % Clinochlore/Chlorite MgChlor 12
                        8.12,       0,      29.43,	0,       5.49,       0,     12.69,	0,          0; ...	     % Chamosite/Chlorite FeChlor 13
                        9.98,       0,      18.85,  0,       10.94,      0,     12.64,	0,          0;	...      % Chlorite 14
                        9.83,   	0.73,	0,      0,       0,          0,     20.46,	0,          0;	...      % Smectite 15
                        11.53,   	0,      2.98,   0,       1.95,       0,     29.26,	0,          0;...        % Smectite2/montmorillonite 16
                        0,          0,      0,      0,       0,          0,     0,      59.93,      0;	...      % Anatase 17
                        0,          0,      46.55,	0,       0,          53.45,	0,      0,          0;	...      % Pyrite 18
                        0,          0,      69.94,	0,       0,          0,     0,      0,          0;	...      % Hematite 19
                        8.11,       0,      7.27,   0,       3.41,       0,     26.73,  0,          5.88;	...  % IlliteLK 20
                        11.16,      0,      5.46,   0,       2.56,       0,     25.34,  0,          6.86;	...  % IlliteMK 21
                        17.26,      0,      1.82,   0,       0.85,       0,     22.55,  0,          8.83; ...    % IlliteHK 22
                        14.21,      0,      3.65,   0,       1.71,       0,     23.95,  0,          7.84];	...  % Illite 23                    

                    

% Elements =    {  'Al',    'Ca',   'Fe',   'Gd',   'Mg',   'S',    'Si',   'Ti',   'K'}; 
uncertainties=  [         5,          5,      5,      1,       1,			 5,     3,      5,      3]	;


Elements =     {        'Al',       'Ca',   'Fe',   'Gd',   'Mg',        'S',   'Si',   'Ti', 'K'}; 
Minerals = {'Calcite';'Dolomite';  'Anhydrite';'Quartz';'Kspar';'Plagioclase';'Albite';'Oligoclase';'Anorthite';'Siderite';...%
    'Kaolinite';'Clinochlore';'Chamosite';'Chlorite';'Smectite';...
    'Smectite/montmorillonite';'Anatase';'Pyrite';'Hematite';...
    'IlliteLK';'IlliteMK';'IlliteHK';'Illite'};

mineralTable = array2table(datat,'RowNames',Minerals,'VariableNames',Elements);
 
end
   
 
function [mineralWF, ElementWF] = createCompositionRand(mineralTable,  mineralPos, num, uncertainties)
mineralWF = zeros(size(mineralTable,1),num); % preallocation

for b=1:num

    posi=randperm(18,1);
    rselection=[0.1:0.1:1,1,1,1,1,1,1,1,1]; 
    rselect=rselection(1,posi);

    if rselect==1 %pure mineral samples
       mineralWF(randperm(size(mineralTable,1),1),b)=1;

    else %random composition samples 
        %create a composition that sums to 1
        h=randfixedsum(6,1,(1-rselect),0,(1-rselect));
        fix=[rselect h'];

        %randomly assign composition to minerals
        mineralWF(randperm(size(mineralTable,1),7),b)=fix;
      
    end
end

if nargout>1
        %Calculate corresponding elemental composition
        ElementWF = ElementWFFromMinerals(mineralTable, mineralWF, uncertainties); % calculate elemental weights from mineral weights
    else
        ElementWF = nan;
    end
    
    %Transpose weight fraction matrix
    mineralWF=mineralWF';
end


                
function [ElementWF] = ElementWFFromMinerals(mineralTable, mineralWF, uncertainties)
    mineralTable = mineralTable{:,:}; % mineral definition table converted
    sy = size(mineralWF,2); %get the number of elements
    ElementWF = zeros(sy,9); %pre-allocation
    for i=1:sy
        %Calculate random noise based on elemental uncertainties
        err = (-1 + (2)*rand(1,length(uncertainties))).*uncertainties;

        %Calculate element weigh fraction plus noise
        ElementWF(i,:) = (mineralTable'*mineralWF(:,i) + err')./100;

    end  
    
    %to set the negative values to zero
    ElementWF(ElementWF<0) = 0;

end



function [TargetMineralWeights, MineralNames]= createTargetWeights(mineralWF)


%Sum the mineral weight fractions to create groups of correlated compositions
TargetMineralWeights = [mineralWF(:,1:4) sum(mineralWF(:,[5 20:23]),2) sum(mineralWF(:,[6:9 11 15:16]),2)  ...
    mineralWF(:,10) sum(mineralWF(:,12:14),2) mineralWF(:,17:19)];

%Final Minerals
MineralNames = ['Calcite' 'Dolomite'  'Anhydrite' 'Quartz' 'Kspar_Illite' 'Plagioclase_Kaolinite_Smect' 'Siderite'   ...
    'Chlorite'  'Anatase' 'Pyrite' 'Hematite' ];

end

function [densityTable]=CreateMineralDensitiesTable()
%Create a table of known mineral densities
MinDen=[	2.71;... %Calcite
            2.84;... %Dolomite
            2.97;... %Anhydrite
            2.66;... %Quartz
            2.56;... %Kspar/Orthoclase
            2.63; ...%Plagioclase
            2.62;... %Albite
            2.65;... %Oligoclase
            2.73;... %Anorthite
            3.95;... %Siderite
            2.65;... %Kaolinite
            2.71;... %MgChlor
            3.21;... %FeChlor
            2.94;... %Chlorite
            2.35;... %Smectite/montmorillonite
            2.74;... %Smectite2
            3.95;... %Anatase
            5.05;... %Pyrite
            5.37;... %Hematite
            2.82;... %IlliteLK
            2.82;... %IlliteMK
            2.81;... %IlliteHK
            2.82];... %Illite

Density={'Density'};
    %Mineral names
    Minerals = {'Calcite';'Dolomite'; 'Anhydrite';'Quartz';'Kspar';...
   'Plagioclase';'Albite';'Oligoclase';'Anorthite';'Siderite'; ...
   'Kaolinite';'Clinochlore';'Chamosite';'Chlorite';...
    'Smectite/montmorillonite';'Smectite2';'Anatase';'Pyrite'; ...
    'Hematite';'IlliteLK';'IlliteMK';'IlliteHK';'Illite'};

densityTable= array2table(MinDen,'RowNames',Minerals,'VariableNames',Density);

end


function [mineralVF,Rhob]=CreatemineralVF(densityTable,mineralWF)
    %Convert the mineral density table to array and transpose
    densityTable=table2array(densityTable);
    densityTable=densityTable';

    %number of minerals and number of samples
    nominer=size(mineralWF,2);
    nosamples=size(mineralWF,1);

    %pre-allocation
    mineralVol=zeros(nosamples,nominer);
    Vmatrix=zeros(nosamples,1);
    Rhoma=zeros(nosamples,1);  
    Rhob=zeros(nosamples,1);  
    Poro=zeros(nosamples,1);
    Vrock=zeros(nosamples,1);
    mineralVF=zeros(nosamples,nominer);

    for k=1:nosamples
        %assume 1 g mass of matrix
        for u=1:nominer %minerals
            %Calculate the mineral volume in 1 g of matrix
            mineralVol(k,u) = (mineralWF(k,u)/densityTable(1,u));            
        end
        
        %Calculate volume of matrix
        Vmatrix(k,1)=sum(mineralVol(k,:));

        %Calculate matrix density
        Rhoma(k,1)=1./Vmatrix(k,1);

        %Calculate fluid density error
        Rhoferror(k,1)=(1-(rand(1)*2)).*0.25;

        %Calculate fluid density
        Rhof(k,1)=1+Rhoferror(k,1);

        %Calculate the bulk density 
        Rhob(k,1)=rand(1).*Rhoma(k,1);

        %Calculate corresponding porosity
        Poro(k,1)=(Rhoma(k,1)-Rhob(k,1))/(Rhoma(k,1)-Rhof(k,1));

        %ensure sensible porosity-otherwise calculate density again
        while Poro(k,1) > 0.4
            Rhob(k,1)=rand(1).*Rhoma(k,1);
            Poro(k,1)=(Rhoma(k,1)-Rhob(k,1))/(Rhoma(k,1)-Rhof(k,1));
            if Poro(k,1) < 0.4
                break
            end
        end
        
        %Calculate rock volume
        Vrock(k,1)=Vmatrix(k,1)/(1-Poro(k,1));

        %Calculate mineral volume fractions
        for r=1:nominer
            mineralVF(k,r) = mineralVol(k,r)/Vrock(k,1);
        end
    end
end


function [TargetMineralVolumes]=createTargetVolumes(mineralVF)

    %Sum the mineral volume fractions to create groups of correlated compositions
    TargetMineralVolumes = [mineralVF(:,1:4) sum(mineralVF(:,[5 20:23]),2) sum(mineralVF(:,[6:9 11 15:16]),2)  ...
    mineralVF(:,10) sum(mineralVF(:,12:14),2) mineralVF(:,17:19)];% mineralVF(:,24)

%Final minerals
MineralNames = ['Calcite' 'Dolomite'  'Anhydrite' 'Quartz' 'Kspar_Illite' 'Plagioclase_Kaolinite_Smect' 'Siderite'   ...
    'Chlorite'  'Anatase' 'Pyrite' 'Hematite'];

end


