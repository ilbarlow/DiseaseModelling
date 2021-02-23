% script to identify candidate genes for developing disease models in worms

% load data
ortholist = readtable('ortholist_master.txt');
wormbase = readtable('wormbase_simplemine_results.txt');

% check probability of having an OMIM phenotype based on how many orthology
% programs agree on ortholog identification.
phenotypeFlag = ~cellfun(@isempty, ortholist.OMIMPhenotypes);
hitProb = NaN(max(ortholist.No_OfPrograms), 1);
for ii = 1:max(ortholist.No_OfPrograms)
    % get number of phenotypes for current program number
    hitProb(ii) = sum(phenotypeFlag(ortholist.No_OfPrograms == ii)) / ...
        sum(ortholist.No_OfPrograms == ii); 
end

% figure
% plot(hitProb)

% find conserved genes that have associated disease mutations in OMIM and
% match particular keywords
% keywords = {'channel', 'glutamate', 'serotonin', 'dopamine', ...
%     'neuropeptide', 'coupled receptor', 'behavior', 'behaviour', ...
%     'acetylcholine'};
keywords = {'neural', 'neuron', 'muscle', 'muscular'};
geneInds = false(size(ortholist, 1), 1);

for ii = 1:numel(keywords)
    % get the matches to the current keywords
    keywordMatchInds = ...
        ~cellfun(@isempty, strfind(wormbase.DescriptionText, keywords{ii}));
    
    % get the corresponding gene IDs
    keywordGeneIDs = wormbase.WormBaseGeneID(keywordMatchInds);
    
    % find the gene IDs in ortholist (note, may be longer because multiple
    % orthologs are found for the same gene)
    ortholistInds = ismember(ortholist.WormBaseID, keywordGeneIDs);
    
    % filter by presence of OMIM phenotype
    ortholistInds = ortholistInds & phenotypeFlag;
    
    % filter by number of orthology programs in agreement
    ortholistInds = ortholistInds & ortholist.No_OfPrograms >= 2;
    
    % add to list of hits
    geneInds = geneInds | ortholistInds;
end

% get subset of ortholist table for selected genes
selectedGenes = ortholist(geneInds, :);

% get the wormbase gene descriptions of the selected genes
[~, wbInds] = ismember(selectedGenes.WormBaseID, wormbase.WormBaseGeneID);
wbDesciption = wormbase.DescriptionText(wbInds);

selectedGenes = [selectedGenes, wbDesciption];

writetable(selectedGenes,'selectedGenesOrtholist2.csv','Delimiter',',') 

%% select genes for specific diseases
disease1 = {'Epilepsy', 'Epileptic', 'Autism', 'Rett',...
    'mental\w*', 'McArdle', 'Cerebral palsy', 'Schizophrenia',...
    'Major depressive disorder', 'personality', 'Bipolar', 'ADHD', ...
    'Attention deficit-hyperactivity disorder', 'seizure', 'parkinson\w*'}

geneticDiseaseInds = false(size(selectedGenes, 1), 1);

for ii = 1:numel(disease1)
    % get the matches to the current keywords
    keywordMatchInds = ...
        ~cellfun(@isempty, regexp(selectedGenes.OMIMPhenotypes, disease1{ii}));
    
    
    
    % add to list of hits
    geneticDiseaseInds = geneticDiseaseInds | keywordMatchInds;
end

%make a refined table of these hits
geneticDiseaseList = selectedGenes(geneticDiseaseInds,:);

writetable(geneticDiseaseList,'geneticDiseases.csv', 'Delimiter', ',');

%% muscular dystrophy filter
disease2 = {'muscular dystrophy'};

muscularDiseaseInds = false(size(selectedGenes, 1), 1);

for ii = 1:numel(disease2)
    % get the matches to the current keywords
    keywordMatchInds = ...
        ~cellfun(@isempty, strfind(selectedGenes.OMIMPhenotypes, disease2{ii}));
    
    % add to list of hits
    muscularDiseaseInds = muscularDiseaseInds | keywordMatchInds;
end

%make a refined table of these hits
muscularDiseaseList = selectedGenes(muscularDiseaseInds,:);

writetable(muscularDiseaseList, 'muscularDystrophy.csv', 'Delimiter', ',');

%% dopamine, serotonin, and drug filter - on OMIM and wormbase description
Monoaminelist = {'dopamin\w*', 'serotonin\w*', '5-HT',  ...
    'GPCR', 'antidepressant', 'antipsychotic'}

MonoamineInds = false(size(selectedGenes, 1), 1);

for ii = 1:numel(Monoaminelist)
    % get the matches to the current keywords
    keywordMatchInds = ...
        ~cellfun(@isempty, regexp(selectedGenes.Var12, Monoaminelist{ii}));
    
    keywordMatchIndsOMIM = ...
        ~cellfun(@isempty, regexp(selectedGenes.OMIMPhenotypes, Monoaminelist{ii}));
    
    % add to list of hits
    MonoamineInds = MonoamineInds | keywordMatchInds | keywordMatchIndsOMIM;
end

%make a refined table of these hits
MonoamineGeneList = selectedGenes(MonoamineInds,:);

writetable(MonoamineGeneList, 'Monoamines.csv', 'Delimiter', ',');

%% ions
ionList = {'sodium', 'potassium', 'chloride', 'proton', 'Kv', 'NMDA', 'AMPA',...
    'iGluR'} %check how case sensitive these searches are 
       
ionInds = false(size(selectedGenes, 1), 1);

for ii = 1:numel(Monoaminelist)
    % get the matches to the current keywords
    keywordMatchInds = ...
        ~cellfun(@isempty, strfind(selectedGenes.Var12, ionList{ii}));
    
    % add to list of hits
    ionInds = ionInds | keywordMatchInds;
end

%make a refined table of these hits
ionGeneList = selectedGenes(ionInds,:);

writetable(ionGeneList, 'Ions.csv', 'Delimiter', ',');

%% neuro and feeding filter

neuroFeedFilter = {'neur\w*', 'synap\w*', 'axon\w*', 'chemotaxis\w*', 'feed\w*'}

neuroFeedInds = false(size(selectedGenes,1),1);
for ii =1:length(neuroFeedFilter)
    neuroFilterWB = ~cellfun(@isempty, regexp(selectedGenes.Var12, neuroFeedFilter{ii}));
    neuroFilterOMIM = ~cellfun(@isempty, regexp(selectedGenes.OMIMPhenotypes, neuroFeedFilter{ii}));
    
    neuroFeedInds = neuroFeedInds | neuroFilterWB | neuroFilterOMIM
end

neuroFeedList = selectedGenes(neuroFeedInds, :);

writetable(neuroFeedList, 'neuroFeed.csv', 'Delimiter', ',');

