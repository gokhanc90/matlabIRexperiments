%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 24);

% Specify sheet and range
opts.Sheet = "MQ07Transposed";
opts.DataRange = "B3:Y1526";

% Specify column names and types
opts.VariableNames = ["MAP_BM25","MAP_DFIC","MAP_DLH13","MAP_LGD","MAP_DLM","MAP_PL2","MAP_DPH","MAP_DFRee",...
    "NDCG20_BM25","NDCG20_DFIC","NDCG20_DLH13","NDCG20_LGD","NDCG20_DLM","NDCG20_PL2","NDCG20_DPH","NDCG20_DFRee",...
    "NDCG100_BM25","NDCG100_DFIC","NDCG100_DLH13","NDCG100_LGD","NDCG100_DLM","NDCG100_PL2","NDCG100_DPH","NDCG100_DFRee",...
    ];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double", "double",...
    ];

% Import the data
MQ07 = readtable("/home/ubuntu/Desktop/Data/Sonuçlar/LovinsRunTopic.xlsx", opts, "UseExcel", false);

MQ07_MAP_BM25.Lovins=MQ07.MAP_BM25;
MQ07_MAP_DFIC.Lovins=MQ07.MAP_DFIC;
MQ07_MAP_DLH13.Lovins=MQ07.MAP_DLH13;
MQ07_MAP_LGD.Lovins=MQ07.MAP_LGD;
MQ07_MAP_DLM.Lovins=MQ07.MAP_DLM;
MQ07_MAP_PL2.Lovins=MQ07.MAP_PL2;
MQ07_MAP_DPH.Lovins=MQ07.MAP_DPH;
MQ07_MAP_DFRee.Lovins=MQ07.MAP_DFRee;

MQ07_NDCG20_BM25.Lovins=MQ07.NDCG20_BM25;
MQ07_NDCG20_DFIC.Lovins=MQ07.NDCG20_DFIC;
MQ07_NDCG20_DLH13.Lovins=MQ07.NDCG20_DLH13;
MQ07_NDCG20_LGD.Lovins=MQ07.NDCG20_LGD;
MQ07_NDCG20_DLM.Lovins=MQ07.NDCG20_DLM;
MQ07_NDCG20_PL2.Lovins=MQ07.NDCG20_PL2;
MQ07_NDCG20_DPH.Lovins=MQ07.NDCG20_DPH;
MQ07_NDCG20_DFRee.Lovins=MQ07.NDCG20_DFRee;


MQ07_NDCG100_BM25.Lovins=MQ07.NDCG100_BM25;
MQ07_NDCG100_DFIC.Lovins=MQ07.NDCG100_DFIC;
MQ07_NDCG100_DLH13.Lovins=MQ07.NDCG100_DLH13;
MQ07_NDCG100_LGD.Lovins=MQ07.NDCG100_LGD;
MQ07_NDCG100_DLM.Lovins=MQ07.NDCG100_DLM;
MQ07_NDCG100_PL2.Lovins=MQ07.NDCG100_PL2;
MQ07_NDCG100_DPH.Lovins=MQ07.NDCG100_DPH;
MQ07_NDCG100_DFRee.Lovins=MQ07.NDCG100_DFRee;

%% Clear temporary variables
clear opts