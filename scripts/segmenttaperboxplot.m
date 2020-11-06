% Script by Ashkan Pakzad to compare intra/intertapering across lobes

%% Load data

AirQuantDir = AirQuantAddPath();
results_dir = fullfile(AirQuantDir,'results');
casenames = ["N1"; "S1"];

cases = cell(size(casenames, 1),1);
for ii = 1:size(casenames, 1)
    dataname = [char(casenames(ii)), '_segmenttaper_0prune.mat'];
    data = load(dataname, 'SegmentTaperResults');
    cases{ii,1} = data.SegmentTaperResults;
end
    

%% Boxplot

caseidx = 2;

figure
taperdata = [cases{caseidx,1}.inner_intra];
typetxt = 'Inner lumen intratapering';
datalabels = [cases{caseidx,1}.lobe];

boxplot(taperdata, datalabels)
xlabel('Lobe')
ylabel('Intratapering (%)')
title(typetxt)

figure
taperdata = [cases{caseidx,1}.inner_inter];
typetxt = 'Inner lumen intertapering';
datalabels = [cases{caseidx,1}.lobe];

boxplot(taperdata, datalabels)
xlabel('Lobe')
ylabel('Intertapering (%)')
title(typetxt)

%% intrataper lumen per lobe - up to gen 6

tiledlayout (1, 2)
ax1 = nexttile;
N1_intramean = intrataper_gen6(cases{1, 1}, 'N1 Inner lumen intratapering');
ax2 = nexttile;
S1_intramean = intrataper_gen6(cases{2, 1}, 'S1 Inner lumen intratapering');
linkaxes([ax1 ax2],'y')



%% intertaper lumen per lobe - up to gen 6

tiledlayout (1, 2)
ax1 = nexttile;
N1_intermean = intertaper_gen6(cases{1, 1}, 'N1 Inner lumen intertapering');
ax2 = nexttile;
S1_intermean = intertaper_gen6(cases{2, 1}, 'S1 Inner lumen intertapering');
linkaxes([ax1 ax2],'y')
ax1.YLim = [-20 60];

%% functions
function tapermeanperlobe = intrataper_gen6(caseinfo, typetxt)
    taperdata = [caseinfo.inner_intra];
    datalabels = [caseinfo.lobe];
    genlabels = [caseinfo.generation];

    % process generations
    genlabels = genlabels <= 6;
    taperdata = taperdata(genlabels);
    datalabels = datalabels(genlabels);

    boxplot(taperdata, datalabels)
    xlabel('Lobe')
    ylabel('Intratapering (%)')
    title(typetxt)
    
    % compute table of avgs
    lobes = {'LL'; 'LU'; 'LUlin'; 'RL'; 'RM'; 'RU'};
    taperperlobe = cell(size(lobes));
    for ii = 1:size(lobes, 1)
    taperperlobe{ii,1} = taperdata(ismember(datalabels, lobes{ii,1}));
    end
    tapermeanperlobe = cellfun(@nanmean, taperperlobe);
    tapermeanperlobe(:,2) = cellfun(@nanstd, taperperlobe);

end

function tapermeanperlobe = intertaper_gen6(caseinfo, typetxt)
    taperdata = [caseinfo.inner_inter];
    datalabels = [caseinfo.lobe];
    genlabels = [caseinfo.generation];

    % process generations
    genlabels = genlabels <= 6;
    taperdata = taperdata(genlabels);
    datalabels = datalabels(genlabels);

    boxplot(taperdata, datalabels)
    xlabel('Lobe')
    ylabel('Intertapering (%)')
    title(typetxt)
    
    % compute table of avgs
    lobes = {'LL'; 'LU'; 'LUlin'; 'RL'; 'RM'; 'RU'};
    taperperlobe = cell(size(lobes));
    for ii = 1:size(lobes, 1)
    taperperlobe{ii,1} = taperdata(ismember(datalabels, lobes{ii,1}));
    end
    tapermeanperlobe = cellfun(@nanmean, taperperlobe);
    tapermeanperlobe(:,2) = cellfun(@nanstd, taperperlobe);
    
end




