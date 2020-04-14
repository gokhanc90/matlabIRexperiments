tw={'BM25' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'LGD' 'PL2'};
for t=1:8
    dataName = strcat('CW12B_NDCG100_',tw{t});
    data = eval(dataName);

    fileID = fopen('oracleSignificant.txt','a');

    stemmers= {'KStemQBS'    'SnowballEngQBS'    'HPS'    'Gupta'    'KStem'    'SnowballEng'};
    for i=1:6
        [p,isSig,oracle] = getOracle(data.NoStem, data.(stemmers{i}) );
        fprintf(fileID,'%.2e\t%d\t%s\n',p,isSig,strcat(dataName,'.NoStem-', dataName, '.', stemmers{i}) );    
    end
    fprintf(fileID,'--------------\n');
end
