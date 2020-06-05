function [Y,Xf,Af] = myNeuralNetworkFunction(X,~,~)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Auto-generated by MATLAB, 29-May-2020 01:36:24.
%
% [Y] = myNeuralNetworkFunction(X,~,~) takes these arguments:
%
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = Qx26 matrix, input #1 at timestep ts.
%
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = Qx1 matrix, output #1 at timestep ts.
%
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1.59970545852425;-1.94533321813139;-1.87640511136435;-0.693129748315317;-1.63355134852371;-0.681385273102391;-1.79432926240141;-1.32576915895429;-0.523088725529093;-0.744573223105749;-2.56953766253936;-1.4170675903003;-0.567169156811807;-0.329296100347013;-0.439484567310285;-0.375597235952408;-0.438257819588213;-1.70793521071043;-0.508807128733938;-0.416503887412968;-0.495453628821827;-1.43343209312919;-3.37675325300613;-1.53803155241897;-0.641748450824542;-2.96311337285318];
x1_step1.gain = [0.715092326593172;0.442096386928316;0.442630947000331;0.312287335678364;0.462550392196808;0.332161804475055;0.46586655953911;0.346771860392247;0.223344255641314;0.477130217466173;0.486580570359142;0.318280384146983;0.591789837449498;0.174584067920577;0.224174603261069;0.218284828639936;0.228408596813643;0.440190736843305;0.21196383117905;0.275471492088977;0.217953976815707;0.560845745648242;0.474225242100022;0.470118304706508;0.280344543547928;0.550936499629357];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.6712633138573507185;-1.0917252687960048441;-0.81198622349552929123;-0.53358560434041391041;0.07269632143675290592;-0.058531769070092970475;0.52680911181116185649;0.93878400299449726862;1.0788328067197749949;-1.5252120349196107707];
IW1_1 = [-0.32784940647282201365 -0.0081326628731234768682 -0.074697948073842335148 -0.44930398619188738341 -0.070025327951656279879 0.11826090409978214091 -0.48837453954650017485 0.044888478828170499946 0.068612253235677295171 0.33588076233536151571 0.26102313514094160896 0.28340231161654977354 -0.16064640024091142445 -0.26703906894382489767 -0.46776766298680183542 -0.47178246445096927442 0.34091215236255106857 0.10843969843376602413 0.31438158976054875193 -0.10948186907251589162 0.13062290987917613005 -0.33345434443978561045 -0.15727250131045081294 0.17191903859039925351 0.036237441147340072878 0.62605547598382438323;0.029794368115344430081 -0.54990605103578371615 -0.11454432351383847521 -0.021280329162950875621 -0.54574690086850485837 -0.10504069258672101006 0.050778829505721684301 0.014158832280385533903 -0.60371197197262560508 0.16945675560262049575 -0.090250773500319167453 -0.46477169239827326175 0.70098985825617021739 -0.40372703973507523267 -0.087246850951725934964 0.45358652140558053389 -0.22936753104756552046 -0.061576598469139723346 0.17124747570600939994 0.34346229574458658762 0.46587625857499381787 0.69360258744324243541 -0.15312719525968956913 0.30942140489662844915 -0.41498871957697675006 0.11853379855919665697;0.67661526723420040152 -0.19973764215108197329 -0.18949054102420359591 0.086707406690115867987 -0.21309216109522371996 0.31158039961408556806 -0.40215383138632326832 -0.41560794614760504517 0.36712796930643259996 0.1012407439130461051 -0.27829031053985736621 0.32388993536659582695 -0.54512683471466283347 0.31388246233827332698 -0.090767901885376200966 -0.16247618694329346578 0.08087548550050993601 -0.40692275937971866817 -0.32309755681651292081 0.54664831600673868817 0.26145118445077092373 -0.49592158134123215385 0.27002314669826232851 0.11157139350027964197 -0.35560540802293361118 -0.25785610917886531146;-0.0042837108958952620913 -0.36076565430990253924 0.21914662152238351123 -0.5256712442128765872 0.42246031335463624101 -0.39685208372177682712 -0.53250396951788525612 -0.06757009064570251855 0.12924142007380950981 0.13684187215030277351 0.47724398890759556213 0.19714275519156973515 -0.085433973725809422839 0.34128343953194656413 -0.32321810315242427425 -0.27303223249929337024 0.32166400907103664952 0.33425987995421158061 0.26229630832546158015 0.24454073727977981845 -0.43887180294661615187 0.072881651234055266908 -0.11357273033321825528 -0.28818208663344369036 -0.077033180914394192018 -0.056521661318949652231;-0.57993595621294058429 -0.012771416713148329297 -0.69934120612275452178 -0.25051668307960189841 -0.73431975528865756075 -0.38892098501772814911 -0.30104793582979016975 -0.51397089231677217658 -0.20311982052750646166 0.69777069948167269597 -0.81679780780631972004 -0.19556585103506918255 -0.21706595917343005087 -0.49863491500843182269 0.21166676953622001744 -0.10108984276401805635 -0.15002885617169064991 -0.37714923142130052369 -0.15344515282086174102 0.18715924497251332226 0.25533068337995323827 0.083577038423569519066 -0.60008333934516577823 -0.092905448636006832119 0.46350394195244926099 -0.067624035111399985465;-0.19317103328780760307 0.24117628252868811289 0.014925047841254070191 0.0087698035586911970984 -0.26362903444727836533 -0.21676202832919655439 0.25606506427162928352 -0.35951389511783632402 0.3703277272559585831 0.030317837155669324911 0.77164968053346172372 0.049737295823668696559 0.063672680898827827001 0.28365976066403680855 0.36606984200394221407 -0.27884417263699673439 -0.47138057316999509849 0.26807905262066539231 0.27270932432427291703 0.53916643008723874608 0.45370237881962077742 0.41314173362843303883 -0.055791057850779930682 0.32478434427884456115 -0.35732099627147934084 0.21745528467548388973;0.034639114692953706931 -0.091269447883429208335 0.4586569658146153472 -0.41336095164794317558 0.068571264976867274399 0.49514174631787483127 0.40889334618766304263 0.10096670616089267747 -0.10210887075789720768 -0.24979891947206386837 -0.28870962394527965866 -0.32176457936169622354 -0.23606092834922562385 -0.0049900759851464835531 0.46142925250843697116 -0.2733154730040741609 0.15603342052364302384 -0.048375872737038111415 0.13863730961969406019 -0.37230752645525688749 -0.39945774617038903864 -0.44969020942089721338 0.097977390936309805003 0.18233586583689381788 0.35027862024134170937 -0.39588216657856101222;0.44780804476194774333 0.38517748592358691528 0.19421542696589053079 0.35971129992252404151 0.0079356405807095437011 -0.45624071157828810508 -0.30944159781871788395 -0.43346722972158469034 -0.46318650583752424277 0.00010246861132306951322 -0.31114815558167124854 0.35435021654089182697 -0.078590965960162670334 -0.49823364820289250954 -0.081706765193707842077 -0.12007876115787269322 -0.22924419611617768977 0.42097618016162746679 0.33491689286310960538 -0.079543234384635466694 0.29234595162363413712 0.43479992292396685283 0.050309128024401891799 -0.19990428476379490386 -0.13908373452117400726 -0.23543523776159872019;0.31482900553730269744 -0.35030257630893668619 -0.12602375578885829932 -0.27102455066768577074 0.083880344136745527139 0.6928497276582665787 0.46266172206142958112 0.11228342306302419029 0.0081363896377381964009 0.50814269248663768153 -0.15505499177023696067 -0.29474442773649695004 0.1571343208294083138 0.30570572824843106341 0.0523713786106151033 -0.35494787588382953336 0.46583893171310175596 0.11728851244958411748 0.59648764710219281948 0.12056067829306163497 0.080912991749497314986 -0.26966775404333637178 -0.566578965776425858 0.46499383249613773783 0.054536882633688089916 -0.058900700614699992574;-0.26927343414790838594 0.35885420146232827054 0.39669819163345293456 0.47407986813218377931 0.21264449824108461184 0.28620798156741655838 -0.12167611237270928304 0.35903936078079201621 0.43030524450820550886 0.28927529326370515772 0.38988104405754397952 -0.16790359017136624598 -0.20113684701444256575 -0.25777854189805360363 0.34392484162961828309 0.02822488140642203866 0.43231411003656949887 -0.18358992879829036382 -0.19813596457018056718 0.14769857293310761914 0.28549714756623356671 0.33088704435838350859 -0.40670874010920154973 -0.0026829975979123813976 -0.23315416229461261022 -0.34842441212381308757];

% Layer 2
b2 = -0.12690698796491389766;
LW2_1 = [-0.761906832128421585 -0.93896408520251417595 -1.0042022546169413477 0.10292720882815623029 1.229282436776809373 -1.4543771204351314541 0.26232361289291872275 -0.51142641646139497258 1.0945849247840278018 0.94883511805656817195];

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
    X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
    Q = size(X{1},1); % samples/series
else
    Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS
    
    % Input 1
    X{1,ts} = X{1,ts}';
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = logsig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Output 1
    Y{1,ts} = a2;
    Y{1,ts} = Y{1,ts}';
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

% Format Output Arguments
if ~isCellX
    Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Positive Transfer Function
function a = logsig_apply(n,~)
a = 1 ./ (1 + exp(-n));
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end