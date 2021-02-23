%script to select the myosin-worm homologs for cardiac disease modelling
%with James Ware.

% load data
ortholist = readtable('ortholist_master.txt');
wormbase = readtable('wormbase_simplemine_results.txt');

%find all the MYH6 and MYH7 homologues
myo  = {'MYH6' 'MYH7', 'MYO6', 'MYO7'};
myoMatchInd= false(size(ortholist,1),1);
for i =1:size(myo,2)
    keyInd = ~cellfun(@isempty, strfind(ortholist.HGNCSymbol, myo{i}));
    
    myoMatchInd = myoMatchInd | keyInd;
    clear keyInd
end

myoMatches = ortholist(myoMatchInd,:);

%add in the wormbase information by finding index of wormbaseIDs that match
[~, wbInds] = ismember(myoMatches.WormBaseID, wormbase.WormBaseGeneID);

wormbaseMyo = wormbase.DescriptionText(wbInds);

wormMyo = [myoMatches, wormbaseMyo];

writetable(wormMyo, 'wormMyosinHomologue.xls');
