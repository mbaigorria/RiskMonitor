%% Implied Risk-Neutral Distribution
% getImpliedDistribution.m (11/06/2015)
%
% This script gets data of CME options on futures and calculates the 
% risk-neutral implied distribution of an asset by using a slight variation 
% of figlewski's method, which is described in 'Estimation the Implied Risk 
% Neutral Density for the U.S. Market Portfolio'. It also generates some
% .txt files which are used later on when compiling a pdf report in LaTeX.
%
% This code can easily be adapted to other type of options.
%
% Dependencies
% 1. JSON toolbox     - loads object from JSON.
% 2. matrix2latex.m   - converts matrix to latex format.
% 3. cell2csv.m       - converts a cell array to csv.
% 4. getDataFromCme.m - gets data from CME.
% 5. mtit.m           - sets general title to multiple subgraphs
%
% Programmed by Martin Baigorria

%% Prepare workspace
clear all;
close all;

addpath('jsonlab'); % load JSONlab toolbox.

% edit these parameters!
selectPrice = 5;  % select which price to use to compute the implied
                  % distribution.
                  % 2: best bid, 3: best offer, 4: average price, 5: prior
                  % settle.

R = 0.01; % anual risk-free rate
blendThreshold = 30; % threshold to define which IVs are going to be blended
dS = 0.5; % for numeric integration in testing. (use below 0.5) (x*dS = 1)

% degree of polinomial fit via OLS (2 or 4 works)
fitDegree = 4;

% gev tail fit parameters
leftDif = 0.03;
rightDif = 0.03;

% prob tables
priceRange_metricTon = [10, 20, 30, 40, 50, 60, 80, 100];
priceRange_bushels   = [25, 50, 75, 100, 125, 150, 200, 250];

% verification
dataVerification = true; % outputs verification files
priceCheck = true;       % automatic price verification
% end edit parameters

mkdir('cache')
mkdir('validation');
mkdir('report');
mkdir('LaTeX/corn/images/');
mkdir('LaTeX/corn/tables/');
mkdir('LaTeX/corn/data/');
mkdir('LaTeX/soybean/images/');
mkdir('LaTeX/soybean/tables/');
mkdir('LaTeX/soybean/data/');
mkdir('LaTeX/wheat/images/');
mkdir('LaTeX/wheat/tables/');
mkdir('LaTeX/wheat/data/');

delete('validation/*'); % remove old validation files
delete('report/*');
delete('LaTeX/corn/images/*');
delete('LaTeX/corn/tables/*');
delete('LaTeX/corn/data/*');
delete('LaTeX/soybean/images/*');
delete('LaTeX/soybean/tables/*');
delete('LaTeX/soybean/data/*');
delete('LaTeX/wheat/images/*');
delete('LaTeX/wheat/tables/*');
delete('LaTeX/wheat/data/*');

% legacy code (everything is computed in bushels anyway)
selectUnit  = 2;  % 1: metric ton, 2: bushels

% table price range & graph fixed domain limit.
% if x_limit is not chosen correctly, the pdf will be truncated, meaning
% the numeric integration will fail.
if selectUnit == 1 % usd/ton
    x_limit = 1000;
    unit = 'metric ton';
else               % cents/bu
    x_limit = 3000;
    unit = 'bushel';
end

% month codification
key = {'F', 'G', 'H', 'J', 'K', 'M', 'N', 'Q', 'U', 'V', 'X', 'Z'};
v = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
monthCodes = containers.Map(key, v);

%% loop through assets
for z = {'soybean','wheat','corn'},

tic;
    
assetType = z{1};

display(['Processing ', assetType]);

% bushels to metric ton conversion (see CME metric guide)
if strcmp(assetType, 'soybean')
    convertToMetricTon = 36.7437;
elseif strcmp(assetType, 'corn')
    convertToMetricTon = 39.3682;
elseif strcmp(assetType, 'wheat')
    convertToMetricTon = 36.7437;
else
    error('Invalid asset type.');
end

for k = 2:13 % loop through options
selectOption = k;
    
% to test a specific contract
%assetType = 'corn';
%selectOption = 6;

% load CME data.
% format: SP | Best Bid | Best Offer | Average Price | Prior Settle | Volume
try 
    optionData = getDataFromCme(assetType, selectOption);
catch exeption
    disp(exeption.message);
    continue;
end

callData          = optionData.callData; 
putData           = optionData.putData;
monthCode         = optionData.monthCode;
expirationDate    = optionData.expirationDate;
expirationMonth   = optionData.expirationMonth;
daysForExpiration = optionData.daysForExpiration;

month = lower(expirationMonth(1:3));

if daysForExpiration <= 20
    display('Contract expires in less than 20 days. Skipping!');
    continue;
end

% display if the contract is important
importantContract = false;
if strcmp(assetType, 'soybean')
    if strcmp(month, 'may') || strcmp(month, 'nov')
        display('Important contact.');
        importantContract = true;
    end
elseif strcmp(assetType, 'corn')
    if strcmp(month, 'may') || strcmp(month, 'dec')
        display('Important contact.');
        importantContract = true;
    end
else % wheat
    if strcmp(month, 'may') || strcmp(month, 'dec')
        display('Important contact.');
        importantContract = true;
    end
end

% skip serial options. codes not to be skipped are
% encoded in a string inside the first parameter of strfind.
if strcmp(assetType, 'corn')    && isempty(strfind('HKNUZ', monthCode(1))) || ...
  (strcmp(assetType, 'soybean') && isempty(strfind('FHKNQUX', monthCode(1)))) || ...
  (strcmp(assetType, 'wheat')   && isempty(strfind('HKNUZ', monthCode(1))))
    display('Serial option detected. Skipping!');
    continue;
end

if daysForExpiration > 300 % 250 (thought it was working days?)
   display('Contract expires in over a year. Skipping!');
   continue;
end

if isempty(callData) || isempty(putData)
   display('Not enough data!');
   continue;
end

% set prior settle underlying price
if selectPrice == 5
    underlyingPrice = optionData.underlyingPriorSettle;
else
    underlyingPrice = optionData.underlyingPrice;
end

% consecutive values shouldn't have the same price (cme data bug?)
for i = 0:length(callData)-2
	if (callData(end - (i+1), selectPrice) - callData(end - i, selectPrice)) <= 0
        callData(end - i, selectPrice) = NaN;
	end
end
    
for i = 1:length(putData)-1
	if putData(i, selectPrice) - putData(i+1, selectPrice) >= 0
        putData(i, selectPrice) = NaN;
	end
end

% save callData and putData (historic purposes)
optionJson = struct('callData', callData, 'putData', putData);
optionJson = savejson('optionData', optionJson);
optionJson = regexprep(optionJson, '/\s+/S', '');
optionJson = strrep(optionJson, sprintf('\t'), '');
optionJson = strrep(optionJson, sprintf('\n'), '');
optionJson = strrep(optionJson, ' ', '');
optionJson = strrep(optionJson, '_NaN_', '');

% convert and identify units
if selectUnit == 1
    blendThreshold = blendThreshold   * convertToMetricTon / 100;
    underlyingPrice = underlyingPrice * convertToMetricTon / 100;
    callData(:,1:5) = callData(:,1:5) * convertToMetricTon / 100;
    putData(:,1:5) = putData(:,1:5)   * convertToMetricTon / 100;
elseif selectUnit ~= 2
    error('Invalid unit selection.');
end

T = daysForExpiration/250;

%% calculate implied volatility
for i = 1:length(callData)
	% call implied volatility
	if isnan(callData(i,selectPrice))
        callData(i,7) = NaN;
    else
        %callData(i,7) = blsimpv(underlyingPrice, callData(i,1), R, ...
        %    daysForExpiration/250, callData(i,selectPrice), 10,  0,  1e-6, true);
        callData(i,7) = blsimpv(underlyingPrice, callData(i,1), 0, T, ...
            callData(i,selectPrice)*exp(R*T), 10,  0,  1e-6, true);
	end
    % put implied volatility
	if isnan(putData(i,selectPrice)) %|| putData(i,1) > underlyingPrice * 1.3
        putData(i,7) = NaN;
    else
        %putData(i,7) = blsimpv(underlyingPrice, putData(i,1), R, ...
        %    daysForExpiration/250, putData(i,selectPrice), 10,  0,  1e-6, false);
        putData(i,7) = blsimpv(underlyingPrice, putData(i,1), 0, T, ...
            putData(i,selectPrice)*exp(R*T), 10,  0,  1e-6, false);
	end
end

%% implied volatility blending & interpolation
% first we define x_high and x_low, which are used to identify 
% the domain in which we will blend the put and call IVs.
x_low = callData(callData(:,1) >= underlyingPrice - blendThreshold,  1);
x_low = x_low(1);

x_high = callData(callData(:,1) <= underlyingPrice + blendThreshold, 1);
x_high = x_high(end);

blended_iv = zeros(length(callData),2);   % pre-allocate memory.

% start blending
for i = 1:length(callData)
    blended_iv(i,1) = callData(i,1);      % strike price
    if callData(i,1) < x_low              % out of the money puts.
        blended_iv(i,2) = putData(i,7);
    elseif callData(i,1) > x_high         % out of the money calls.
        blended_iv(i,2) = callData(i,7);
    else                                  % blend iv.
        if isnan(callData(i,6))
            blended_iv(i,2) = putData(i,7); % you cannot blend with nans.
        elseif isnan(putData(i,6))
            blended_iv(i,2) = callData(i,7);
        else
            w = (x_high - callData(i,1)) / (x_high - x_low);
            blended_iv(i,2) = w*putData(i,7) + (1-w)*callData(i,7);
        end
    end
end

% volatility smile interpolation
blended_iv   = blended_iv(blended_iv(:,2) > 0,:); % remove nan rows for polyfit.
domain = (0:dS:x_limit)';

[p,S,mu] = polyfit(blended_iv(:,1), blended_iv(:,2), fitDegree);
[smileFit,~] = polyval(p,domain,S,mu);

% only fit smile in available data range (can be made more efficient)
smileFit(domain < blended_iv(1,1) | domain > blended_iv(end, 1), 1) = NaN;

% the containers are larger than the length of the actual data to make
% space to later fit the distribution tails
callPrices = zeros(length(domain), 1);
putPrices  = zeros(length(domain), 1);

%% from implied volatility to price space
% iterate through the implied volatilities to go from the volatility
% to the price space.
for j = 1:length(domain)  % iterate through the generated strike prices.
    
    if isnan(smileFit(j,1))
        callPrices(j,1) = NaN;
        putPrices(j,1) = NaN;
        continue;
    end
    
    [callPrices(j,1), putPrices(j,1)] = blsprice(underlyingPrice, ...
        domain(j,1), 0, daysForExpiration/250, max(smileFit(j,1), 0.0001));
    
    callPrices(j,1) = callPrices(j,1) * exp(-R*T);
    putPrices(j,1)  = putPrices(j,1)  * exp(-R*T);
    
    % calls should have a lower value as the strike price increases
    % j > 1 is to avoid an illegal access.
	if j > 1 && (callPrices(j-1,1) - callPrices(j,1)) < 0 % one price effect beats the other.
        callPrices(j,1) = NaN;
        break;
	end
    
end

% build underlying distribution from the new call and put prices.
callDist = zeros(length(callPrices) - 2, 2);
%putDist  = zeros(length(domain) - 2, 2);

% iterate prices and generate pdf and cdf in our data range by discrete 
% estimation of first and second order partial derivatives. this could be
% improved by only iterating the range where we have data.
for j = 1:length(callPrices) - 2; 

    if isnan(callPrices(j, 1))
        callDist(j,1) = NaN;
        callDist(j,2) = NaN;
    end
    
    callDist(j,1) = exp(daysForExpiration/250*R)*(callPrices(j+2,1) - ...
        2*callPrices(j+1,1) + callPrices(j,1))/(dS^2); % PDF from calls
	callDist(j,2) = exp(daysForExpiration/250*R)*(callPrices(j+2,1) - ...
        callPrices(j,1))/(2*dS) + 1;  % CDF from calls.
    
	%putDist(j,1) = exp(daysForExpiration/250*R)*(putPrices(j+2,1) - ...
    %    2*putPrices(j+1,1) + putPrices(j,1))/dS^2;  % PDF from Puts
 	%putDist(j,2) = exp(daysForExpiration/250*R)*(putPrices(j+2,1) - ...
    %    putPrices(j,1))/2;  % CDF from puts.
       
end

pdf = horzcat(domain(2:length(callPrices) - 1,1), callDist(:,1));
cdf = horzcat(domain(2:length(callPrices) - 1,1), callDist(:,2));

% cut pdf if data behaves badly
cleanRight = pdf(1:end-1,2)-pdf(2:end,2) < 0 & pdf(1:end-1,1) > underlyingPrice * 1.05;
cleanRight = find(cleanRight > 0) - 1;

if ~isnan(cleanRight)
    display('Removed badly behaved right tail from empiric data.');
    cleanRight = cleanRight(1) - 1;
    pdf(cleanRight:end, 2) = NaN;
    cdf(cleanRight:end, 2) = NaN;
end

% careful with finite arithmetic!
cleanLeft = (pdf(2:end,2)-pdf(1:end-1,2) < -10e-5 | pdf(1:end-1,2) < 0) & pdf(1:end-1,1) < underlyingPrice * 0.95;
cleanLeft = find(cleanLeft > 0) + 1;

if ~isnan(cleanLeft)
    display('Removed badly behaved left tail from empiric data.');
    cleanLeft = cleanLeft(end) + 1;
    pdf(1:cleanLeft, 2) = NaN;
    cdf(1:cleanLeft, 2) = NaN;
end

%testing
%figure('Visible','on');
%plot(pdf(:,1), pdf(:,2))
%plot(blended_iv(:,1), blended_iv(:,2), domain, smileFit)

%% GEV (generalized extreme value) distribution fitting.

% this version doesnt work, guess vpasolve doesnt identify the built-in gev
% functions as symbolic, so i will re-write them.
%syms e sigma mu
%S = vpasolve([...
%gevcdf(alpha0R(1,1),e,sigma,mu) == alpha0R(1,2),...
%gevpdf(alpha0R(1,1),e,sigma,mu) == alpha0R(1,3),...
%gevpdf(alpha1R(1,1),e,sigma,mu) == alpha1R(1,3)], [e, sigma, mu]);

% right tail gev

% gev fit parameters
% round is used to avoid issues with unit conversion.
comp = cdf(cdf(:,2) > 0, 2);
comp = comp(end) - 0.01;

upperFitR = min(0.95, comp);
lowerFitR = upperFitR - rightDif;

% format: SP - CDF - PDF
[r0,~] = find(cdf(:,2) < upperFitR);
r0 = max(r0);
alpha0R = horzcat(cdf(r0, :), pdf(r0,2));

[r1,~] = find(cdf(:,2) < lowerFitR);
r1 = max(r1);
alpha1R = horzcat(cdf(r1, :), pdf(r1,2));

syms e sigma mu
S = vpasolve([...
    exp(-(1+e*((alpha0R(1,1)-mu)/sigma))^(-1/e)) == alpha0R(1,2),...
    (1/sigma) * (1+e*((alpha0R(1,1)-mu)/sigma))^(-(1/e)-1) * ...
    exp(-(1+e*((alpha0R(1,1)-mu)/sigma))^(-(1/e))) == alpha0R(1,3),...
    (1/sigma) * (1+e*((alpha1R(1,1)-mu)/sigma))^(-(1/e)-1) * ...
    exp(-(1+e*((alpha1R(1,1)-mu)/sigma))^(-(1/e))) == alpha1R(1,3)], [e, sigma, mu]);

try
    e     = double(S.e);
    sigma = double(S.sigma);
    mu    = double(S.mu);

    pdf(length(domain),:) = 0;
    cdf(length(domain),:) = 0;

    pdf(:,1) = domain;
    cdf(:,1) = domain;

    for j = r0:length(domain) % (pdf(:,1)); %mem issues, testing! len(prices) < len(dom)
        pdf(j,2) = gevpdf(pdf(j,1), e, sigma, mu);
        cdf(j,2) = gevcdf(cdf(j,1), e, sigma, mu);
    end
    
    rightFit = true;
catch
    rightFit = false;
    display('Right tail fit failed');
end

% left tail gev
comp = cdf(cdf(:,2) > 0, 2);
comp = comp(1) + 0.01; % add an extra 1% to param search
lowerFitL = max(0.05, comp);
upperFitL = lowerFitL + leftDif;

% format: SP - CDF - PDF
[r0,~] = find(cdf(:,2) <= lowerFitL);
r0 = r0(end);
alpha0L = horzcat(cdf(r0, :), pdf(r0,2));

[r1,~] = find(cdf(:,2) <= upperFitL);
r1 = r1(end);
alpha1L = horzcat(cdf(r1, :), pdf(r1,2));

syms e sigma mu
S = vpasolve([...
    exp(-(1+e*((-alpha0L(1,1)-mu)/sigma))^(-1/e)) == (1 - alpha0L(1,2)),...
    (1/sigma) * (1+e*((-alpha0L(1,1)-mu)/sigma))^(-(1/e)-1) * ...
    exp(-(1+e*((-alpha0L(1,1)-mu)/sigma))^(-(1/e))) == alpha0L(1,3),...
    (1/sigma) * (1+e*((-alpha1L(1,1)-mu)/sigma))^(-(1/e)-1) * ...
    exp(-(1+e*((-alpha1L(1,1)-mu)/sigma))^(-(1/e))) == alpha1L(1,3)], [e, sigma, mu]);

try
    e     = double(S.e);
    sigma = double(S.sigma);
    mu    = double(S.mu);

    for j = 1:r0;
        pdf(j,2) = gevpdf(-pdf(j,1), e, sigma, mu);
        cdf(j,2) = 1 - gevcdf(-cdf(j,1), e, sigma, mu);
    end
    
    leftFit = true;
catch
    leftFit = false;
    display('Left tail fit failed');
end

pdf = real(pdf); % vpasolve also optimizes in complex axis.
cdf = real(cdf);

%% data verification is done on original CME units (cents per bushel)
if dataVerification == true % after tail fit

    callFile = fopen('validation/callPriceValidation.txt', 'a');
    putFile =  fopen('validation/putPriceValidation.txt' , 'a');

    lower_cutoff = find(~isnan(pdf(2:end,2)));% en algun lado le agregue una primera fila
    lower_cutoff = lower_cutoff(1);

    upper_cutoff = find(isnan(pdf(2:end,2)) & pdf(1:end-1,1) > underlyingPrice); % en algun lado le agregue una primera fila
    if isempty(upper_cutoff)
        upper_cutoff = length(pdf);
    else
        upper_cutoff = upper_cutoff(1);
    end

    % to test cutoffs & tails (useful to find bugs)
    %{
    figure('Visible','on');
    hold on
    title(strcat(assetType, {' '}, num2str(k), {' - '}, expirationMonth));
    xlabel(strcat({'Future Price (USD per '},unit,')'));
    plot(pdf(:,1), pdf(:,2));
    line([callData(1,1) callData(1,1)], get(gca, 'ylim'), 'Color', 'green');
    line([callData(end,1) callData(end,1)], get(gca, 'ylim'), 'Color', 'green');
    %refline(0,0);
    %hline.Color = 'black';
    line([pdf(lower_cutoff,1) pdf(lower_cutoff,1)], get(gca, 'ylim'), 'Color', 'red');
    line([pdf(upper_cutoff,1) pdf(upper_cutoff,1)], get(gca, 'ylim'), 'Color', 'red');
    hold off
    %}
    %k = waitforbuttonpress 
    
    % data verification
    % calls
    compCall = callData(~isnan(callData(:, 5)), [1,5]);
    
    % for each contract, inner loop calculates empiric price.
    for i = 1:length(compCall)

        callSum = 0;

        for j = lower_cutoff:upper_cutoff
            callSum = callSum + (exp(-R*daysForExpiration/250) * ...
                max(pdf(j,1) - compCall(i,1), 0) * max(pdf(j,2),0)) * dS;
        end

        compCall(i, 3) = callSum;

    end
    
    % puts
    compPut = putData(~isnan(putData(:, 5)), [1,5]);
    for i = 1:length(compPut)

        putSum = 0;

        for j = lower_cutoff:upper_cutoff
            putSum = putSum + (exp(-R*daysForExpiration/250) * ...
                max(compPut(i,1) - pdf(j,1), 0) * max(pdf(j,2),0)) * dS;
        end

        compPut(i, 3) = putSum;

    end
    
    if selectUnit == 1 % restore original units (bushels) (legacy)
        underlyingPrice = underlyingPrice / convertToMetricTon * 100; % restore underlying price
        
        compCall = compCall ./ convertToMetricTon .* 100;
        compCall(:,1) = floor(compCall(:,1));
        
        compPut = compPut ./ convertToMetricTon .* 100;
        compPut(:,1) = floor(compPut(:,1));
    end
    
    % can be improved pre-allocating memory.
    compCall(:,4) = compCall(:,3) - compCall(:,2);
    compCall(:,5) = abs(compCall(:,2)./compCall(:,3)-1);
    
    compPut(:,4) = compPut(:,3) - compPut(:,2);
    compPut(:,5) = abs(compPut(:,2)./compPut(:,3)-1);
    
    % check martingale
    checkMartingale = 0;
    for j = lower_cutoff:upper_cutoff
        checkMartingale = checkMartingale + pdf(j,1)*max(pdf(j,2),0) * dS;
    end;
    
    if selectUnit == 1
        checkMartingale = checkMartingale / convertToMetricTon * 100;
    end;
    
    %thresh = 100 / convertToMetricTon * 100; % 100 dolar range to cents per bushel
    thresh = 100;
    
    % filter prices
    compCall = compCall(compCall(:,1) > underlyingPrice - thresh & ...
        compCall(:,1) < underlyingPrice + thresh, :);
    compPut = compPut(compPut(:,1) > underlyingPrice - thresh & ...
        compPut(:,1) < underlyingPrice + thresh, :);
    
    fwrite(callFile, [upper(assetType(1)) assetType(2:end) ': ' expirationMonth  ' (' monthCode ')' char(10)]);
    fprintf(callFile, 'underlyingPrice: %5.0f\n', underlyingPrice);
    
    fwrite(putFile, [upper(assetType(1)) assetType(2:end) ': ' expirationMonth  ' (' monthCode ')' char(10)]);
    fprintf(putFile, 'underlyingPrice: %5.0f\n', underlyingPrice); 
    
    fprintf(callFile, 'checkMartingale: %5.0f\n', checkMartingale);
    fprintf(putFile, 'checkMartingale: %5.0f\n', checkMartingale);
    
    if rightFit == false
        fprintf(callFile, 'Right tail fit failed.\n');
    end
    
	if leftFit == false
        fprintf(putFile, 'Left tail fit failed.\n');
    end
    
    container = compCall(compCall(:,1) > underlyingPrice - thresh & ...
        compCall(:,1) < underlyingPrice + thresh, :);
    
    % you might wonder why you do not get prices for every strike price
    % in the original data. this is because the prices are calculated
    % straight from the volatility smile fit. because we only consider
    % out of the money options (both calls and pulls) for this fit, this
    % does not imply we will have a volatility value for every available
    % strike price from a call or put. this means we wont be able to
    % calculate the price for every strike price in the original data.
    
    % get corresponding calculated prices from fitted vol smile
    prices = zeros(length(container),1);
    
    for i = 1:length(prices)
        
        index = find(domain == container(i, 1));
        index = index(1);
        
        prices(i) = callPrices(index);
    end
    
    fwrite(callFile, ['Strike      Original       P_smile  Recalculated   Delta      %' char(10)]);
    fprintf(callFile, '%6.0f  %12.2f  %12.2f  %12.2f  %6.2f  %3.2f%% \n', [container(:,1:2), prices,container(:,3:5)]');
   
    fwrite(callFile, char(10));
    
    container = compPut(compPut(:,1) > underlyingPrice - thresh & ...
        compPut(:,1) < underlyingPrice + thresh, :);
    
    prices = zeros(length(container),1);
    
    for i = 1:length(prices)
        
        index = find(domain == container(i, 1));
        index = index(1);
        
        prices(i) = putPrices(index);
    end
    
    fwrite(putFile, ['Strike      Original       P_smile  Recalculated   Delta      %' char(10)]);
    fprintf(putFile, '%6.0f  %12.2f  %12.2f  %12.2f  %6.2f  %3.2f%% \n', [container(:,1:2), prices,container(:,3:5)]');
    fwrite(putFile, char(10));

    fclose(callFile);
    fclose(putFile);
    
    % restore underlying price
    if selectUnit == 1
        underlyingPrice = underlyingPrice * convertToMetricTon / 100;
    end

end

%% price checks are done on original CME units (bushels)
if priceCheck == true
    
    % validate martingale
    if priceCheck == true
        if abs(checkMartingale/underlyingPrice -1) > 0.1
            display('Martingale check failed.');
            continue
        end
    end
    
	stop = false;

    n = 0;
    avg = 0;
    
	for i = 1:length(compCall)

        if (compCall(i,1) >= underlyingPrice*0.8 && compCall(i,1) <= underlyingPrice*1.2) && compCall(i,2) > 10
            n = n + 1;
            avg = avg + abs(compCall(i,2) / compCall(i,3) - 1);
        end

	end

    avg = avg / n;
    
    display(['Call average distance %: ', num2str(avg)]);
    
    % right fit failed, we have more tolerance on average error since not
    % all the integral could be computed.
    if rightFit == false
        bound = 0.05;
    else
        bound = 0.04;
    end
    
    if avg > bound && ~importantContract
        display('Control diagnosis failed for Calls. Skipping!');
        stop = true;
    end

	if stop == true
        continue; % go to next asset
	end

    n = 0;
    avg = 0;
    
	for i = 1:length(compPut)

        if (compPut(i,1) >= underlyingPrice*0.8 && compPut(i,1) <= underlyingPrice*1.2) && compPut(i,2) > 10
            n = n + 1;
            avg = avg + abs(compPut(i,2) / compPut(i,3) - 1);
        end

	end

    avg = avg / n;
    
    display(['Put average distance %: ', num2str(avg)]);
    
    if leftFit == false
        bound = 0.05;
    else
        bound = 0.03;
    end
    
    if avg > bound && ~importantContract
        display('Control diagnosis failed for Puts. Skipping!');
        stop = true;
    end
    
	if stop == true
        continue;
 	end
end

%% Compute probability tables
%expectedValue = pdf(pdf(:,2) == max(pdf(:,2)),1);
expectedValue = underlyingPrice;

% fixed price variation
tableIncrease_bushels = zeros(length(priceRange_bushels), 2);
tableIncrease_bushels(:, 1) = priceRange_bushels;
tableDrop_bushels = zeros(length(priceRange_bushels), 2);
tableDrop_bushels(:, 1) = priceRange_bushels;
for j=1:length(tableIncrease_bushels); % prob increase
	aux = cdf(cdf(:,1) <= expectedValue + tableIncrease_bushels(j,1),2);
	tableIncrease_bushels(j,2) = 100 - aux(end)*100;
	aux = cdf(cdf(:,1) >= expectedValue - tableDrop_bushels(j,1),2);
	if isnan(aux(1)) % quick fix, change this later to do < x %
        aux(1) = -1/100;
	end
    tableDrop_bushels(j,2) = aux(1)*100;
end

% fixed price variation
tableIncrease_metricTon = zeros(length(priceRange_metricTon), 2);
tableIncrease_metricTon(:, 1) = priceRange_metricTon;
tableDrop_metricTon = zeros(length(priceRange_metricTon), 2);
tableDrop_metricTon(:, 1) = priceRange_metricTon;
for j=1:length(tableIncrease_metricTon); % prob increase
	aux = cdf(cdf(:,1) <= expectedValue + tableIncrease_metricTon(j,1)/convertToMetricTon*100,2);
	tableIncrease_metricTon(j,2) = 100 - aux(end)*100;
	aux = cdf(cdf(:,1) >= expectedValue - tableDrop_metricTon(j,1)/convertToMetricTon*100,2);
	if isnan(aux(1))
        aux(1) = -1/100;
	end
	tableDrop_metricTon(j,2) = aux(1)*100;
end

%% save graphs, data & results.
columnLabels = {'Increase greater than (in USD)', 'Probability'};
matrix2latex(tableIncrease_bushels, strcat('LaTeX/', assetType, '/tables/suba', ...
    int2str(selectOption),'.txt'), 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f');

columnLabels = {'Drop greater than (in USD)', 'Probability'};
matrix2latex(tableDrop_bushels, strcat('LaTeX/', assetType, '/tables/baja', ...
    int2str(selectOption),'.txt'), 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f');

% save graph for latex/word report generation.
h = figure('Visible','off');
subplot(2,1,1);
plot(pdf(:,1),pdf(:,2),'-r');
ylim([0 max(ylim)]);
xlabel(strcat({'Future Price (USD per '},unit,')'));
ylabel('Density');
%legend('Underlying Density','Location','northwest')
title('Risk Neutral Density');

subplot(2,1,2);
plot(cdf(:,1),cdf(:,2),'-r');
xlabel(strcat({'Future Price (USD per '},unit,')'));
ylabel('Probability');
%legend('Underlying Distribution','Location','northwest')
title('Risk Neutral Distribution');
saveas(h, strcat('LaTeX/', assetType,'/images/dist',int2str(selectOption), '.eps'), 'eps');
saveas(h, strcat('report/', assetType, '-', lower(expirationMonth), '.eps'), 'eps');

%% density graph with horizontal indicators (metric tons)
% comment: using the eps or the png engines to generate the graphs give
% you completely different results.

% change of variable to the pdf/cdf
% pdf is used to find the height of the label.
% cdf is used to find strike price cutoff.
pdf(:,1) = pdf(:,1) * convertToMetricTon / 100;
pdf(:,2) = pdf(:,2) / convertToMetricTon * 100;
cdf(:,1) = cdf(:,1) * convertToMetricTon / 100;

h = figure('Visible','off');
hold on; % build the whole figure first

distFilter = zeros(length(pdf), 1);
distFilter(lower_cutoff:upper_cutoff) = 1;
filterRows = pdf(:,2) >= max(pdf(:,2)) * 0.01 & distFilter;

% get ticks
aux = cdf(distFilter & cdf(:,2) >= 0.1, 1);
strike_10  = aux(1);
density_10 = pdf(pdf(:,1) == strike_10, 2);

aux = cdf(distFilter & cdf(:,2) >= 0.25, 1);
strike_25  = aux(1);
density_25 = pdf(pdf(:,1) == strike_25, 2);

aux = cdf(distFilter & cdf(:,2) >= 0.75, 1);
strike_75  = aux(1);
density_75 = pdf(pdf(:,1) == strike_75, 2);

aux = cdf(distFilter & cdf(:,2) >= 0.9, 1);
strike_90  = aux(1);
density_90 = pdf(pdf(:,1) == strike_90, 2);

% area fill
%{
filter1 = cdf(:,1) <= strike_10 & filterRows;
area(pdf(filter1, 1), pdf(filter1, 2), 'EdgeColor', 'none', 'FaceColor', [170/255 240/255 152/255]);

filter2 = cdf(:,1) > strike_10 & cdf(:,1) <= strike_25 & filterRows;
area(pdf(filter2, 1), pdf(filter2, 2), 'EdgeColor', 'none', 'FaceColor', [150/255 230/255 152/255]);

filter3 = cdf(:,1) > strike_25 & cdf(:,1) <= strike_75 & filterRows;
area(pdf(filter3, 1), pdf(filter3, 2), 'EdgeColor', 'none', 'FaceColor', [120/255 220/255 152/255]);

filter4 = cdf(:,1) > strike_75 & cdf(:,1) <= strike_90 & filterRows;
area(pdf(filter4, 1), pdf(filter4, 2), 'EdgeColor', 'none', 'FaceColor', [150/255 230/255 152/255]);

filter5 = cdf(:,1) > strike_90 & filterRows;
area(pdf(filter5, 1), pdf(filter5, 2), 'EdgeColor', 'none', 'FaceColor', [170/255 240/255 152/255]);
%}

% save cutpoints in original units (bushels) -> (x,y) = (strike, density)
cut_points = [strike_10 / convertToMetricTon * 100  density_10 * convertToMetricTon / 100; ...
              strike_25 / convertToMetricTon * 100  density_25 * convertToMetricTon / 100; ...
              strike_75 / convertToMetricTon * 100  density_75 * convertToMetricTon / 100; ...
              strike_90 / convertToMetricTon * 100  density_90 * convertToMetricTon / 100];

plot(pdf(filterRows,1),pdf(filterRows,2),'-blue', ...
    [strike_10, strike_10], [0, density_10], '-black', ...
    [strike_25, strike_25], [0, density_25], '-black', ...    
    [strike_75, strike_75], [0, density_75], '-black', ...
    [strike_90, strike_90], [0, density_90], '-black');

% plot config (failed to hide y axis, I did try!)
ylim([0 max(ylim)]);
%ylim([0 0.011]);
set(gca,'yTick',[]);

try
    set(gca,'xTick', [str2num(sprintf('%.0f',strike_10)) str2num(sprintf('%.0f',strike_25)) ...
        floor(underlyingPrice*convertToMetricTon/100) str2num(sprintf('%.0f',strike_75)) str2num(sprintf('%.0f',strike_90))]);
catch
    display('Low variance distribution caused tick values to overlap. Ticks ignored.');
end

set(gca,'LooseInset',get(gca,'TightInset')); 
%set(gca, 'FontName', 'Arial')
set(gca,'box','off');
set(gca,'FontSize',10);
set(gca,'YColor','w');
set(gca,'TickLength',[ 0 0 ])
%set(findall(h, 'Type','text'), 'FontSize', 25);

% to fix font size issue in Linux (Matlab bug?)
% sudo apt-get install xfonts-base xfonts-100dpi xfonts-75dpi
% restart after!

% end plot config
labelSize = 13;

% add labels
text(strike_10*0.99, (density_10 / 5) , ['10%' char(10) 'chance'],'HorizontalAlignment','right','FontSize',labelSize)
%text((strike_10 + strike_20) / 2, (density_10 + density_20) / 2, '10% \rightarrow','HorizontalAlignment','right','FontSize',labelSize)
text((strike_10 + strike_25) / 2 , max(pdf(filterRows,2))/5 , ['15%' char(10) 'chance'],'HorizontalAlignment','center','FontSize',labelSize)
text((strike_75 + strike_90) / 2 , max(pdf(filterRows,2))/5 , ['15%' char(10) 'chance'],'HorizontalAlignment','center','FontSize',labelSize)
text((strike_25 + strike_75) / 2,  max(pdf(filterRows,2))/2 , ['50%' char(10) 'chance'],'HorizontalAlignment','center','FontSize',labelSize)
%text((strike_80 + strike_90) / 2, (density_80 + density_90) / 2, '\leftarrow 10%','HorizontalAlignment','left','FontSize',labelSize)
text(strike_90 + 5, (density_90 / 4) , ['10%' char(10) 'chance'],'HorizontalAlignment','left','FontSize',labelSize)

% language config. to generate these graphs in english some changes are
% required. we do some language translations as well.
%feature('DefaultCharacterSet', 'UTF-8'); % try UTF-8 Windows-1250

if strcmp(assetType, 'corn')
    assetSp = 'maiz';
elseif strcmp(assetType, 'soybean')
    assetSp = 'soja';
else
    assetSp = 'trigo';
end

%xlabel(strcat([upper(assetType(1)),assetType(2:end),' price at expiration date (USD cents per ',unit,')']));
xlabel(strcat(['D',char(243),'lares por tonelada de ',assetSp,' al vencimiento']),'interpreter','tex');

key = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
v = {'Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio', ...
  'Julio', 'Agosto', 'Septiembre', 'Octubre', 'Noviembre', 'Diciembre'};
months = containers.Map(key, v);

dayData = strsplit(expirationMonth);

% matlab 2014 doesn't have good date functions, this is to get yesterday's
% date. pretty ugly but it works...
currentDateData = strsplit(date,'-');

%title({[upper(assetType(1)) assetType(2:end) ': ' expirationMonth], 'Risk Neutral Density'})
title({['Escenarios de precios de ' assetSp ' en ' months(dayData{1}) ...
    ' ' dayData{2} ' de acuerdo al mercado de Chicago el ' datestr(datenum(str2num(currentDateData{3}),find(strcmp(key, currentDateData{2})),str2num(currentDateData{1}))-1) ]});
print(h, strcat('report/', assetType, '-', int2str(k), '-', ...
    lower(expirationMonth), '_density'), '-dpng', '-r300');

% UNIX command to generate gif.
% convert -delay 100 -loop 0 *_density.eps animation.gif

% restore pdf/cdf
pdf(:,1) = pdf(:,1) / convertToMetricTon * 100;
pdf(:,2) = pdf(:,2) * convertToMetricTon / 100;
cdf(:,1) = cdf(:,1) / convertToMetricTon * 100;

%% Save reports/data
% save option variables for word report
fileID = fopen('report/data.txt', 'a');
fwrite(fileID, [
    upper(assetType(1)) assetType(2:end) ': ' expirationMonth  ' (' monthCode ')' char(10)]);
fprintf(fileID, 'Expiration date: %s\n', expirationDate);
fprintf(fileID, 'Expected price per %s: %3.0f USD\n', unit, expectedValue);
fprintf(fileID, 'Working days till expiration: %3.0f\n', daysact(date,expirationDate));
fprintf(fileID, 'Risk-free rate: %6.3f\n', R);
fprintf(fileID, 'underlyingPrice: %5.0f\n\n', underlyingPrice);
fprintf(fileID, 'Strike     P(inc)    P(drop)\n');
fprintf(fileID, '%6.0f  %8.2f%% %8.2f%% \n', horzcat(tableIncrease_bushels, tableDrop_bushels(:,2))');
fprintf(fileID, '\n');
fclose(fileID);

% save option variables for latex report
fileID = fopen(strcat('LaTeX/', assetType,'/data/report.txt'), 'a');
fprintf(fileID, strcat('\\input{data/data', int2str(selectOption),'.txt \n}', ...
'\\begin{center} \n', ...
'\\includegraphics[scale=0.8]{images/dist', int2str(selectOption),'.eps} \n', ...
'\\end{center} \n', ...
'\\begin{table}[h] \n', ...
'\\begin{tabular}{ll} \n', ...
'Price Increase Scenario & Price Drop Scenario \\\\ \n', ... 
'\\input{tables/suba', int2str(selectOption),'.txt} & \\input{tables/baja', int2str(selectOption),'.txt} \n', ...
'\\end{tabular} \n', ...
'\\end{table} \n', ...
'All changes are computed taking the expected future price as a start point. \n', ...
'\\newpage \n\n'));
fclose(fileID);

% save files for latex report
fileID = fopen(strcat('LaTeX/', assetType,'/data/data',int2str(selectOption),'.txt'), 'w');
fprintf(fileID, '\\section{ %s } Expiration date: %s \\\\ Expected price per %s: %3.0f USD \\\\ Working days till expiration: %3.0f \\\\ Risk-free rate: %6.3f', ...
    expirationMonth, expirationDate, unit, expectedValue, daysact(date,expirationDate), R);
fclose(fileID);

% save diagnosis graphs
h = figure('Visible','off');

subplot(4,2,1);
plot(callData(:,1),callData(:,7),'--xg',putData(:,1),putData(:,7),'-ob',domain,smileFit,'--b');
xlabel('Strike Price');
ylabel('Implied Volatility');
%legend('Call IVs','Put IVs','2nd degree polynomial','Location','southoutside','Orientation','horizontal');
title('Empirical Volatility Smile Fit');

subplot(4,2,2);
scatter(callData(:,1), callData(:, selectPrice));
hold on
scatter(putData(:,1), putData(:, selectPrice));
hold on
plot(domain(1:length(callPrices)),callPrices(:,1),'-r', ...
    domain(1:length(putPrices)),putPrices(:,1),'-b');
hold off
xlabel('Strike Price');
ylabel('Option Price');
%legend('Call Price','Put Price','Location','southoutside','Orientation','horizontal');
title('Black-Scholes Prices');

subplot(4,2,3);
plot(pdf(:,1),pdf(:,2),'-r');
xlabel('Future Price');
ylim([0 max(ylim)]);
ylabel('Density');
%legend('Underlying Density');
title('Risk Neutral Density');

subplot(4,2,4);
plot(cdf(:,1),cdf(:,2),'-r');
xlabel('Future Price');
ylabel('Probability');
%legend('Underlying Distribution');
title('Risk Neutral Distribution');

try
        
    subplot(4,2,5);
    if selectUnit == 1
        filter = compCall(:,1) >= underlyingPrice/ convertToMetricTon * 100 * 0.9;
    else
        filter = compCall(:,1) >= underlyingPrice * 0.9;
    end

    plot(compCall(filter,1), ...
         compCall(filter,2),'-r', ...
         compCall(filter,1), ...
         compCall(filter,3),'-b');
    xlim([floor(min(compCall(filter,1))) ceil(max(compCall(filter,1)))]);
    ylim([floor(min(min(compCall(filter,[2,3])))) ceil(max(max(compCall(filter,[2,3]))))]);
    xlabel('Strike Price');
    ylabel('Price');
    %legend('Empiric Prices', 'Numeric Prices');
    title('Call Prices');

catch
    
end
 
subplot(4,2,6);
if selectUnit == 1
    filter = compPut(:,1) <= underlyingPrice/ convertToMetricTon * 100 * 1.1;
else
    filter = compPut(:,1) <= underlyingPrice * 1.1;
end
plot(compPut(filter,1), ...
     compPut(filter,2),'-r', ...
     compPut(filter,1), ...
     compPut(filter,3),'-b');
xlim([floor(min(compPut(filter,1))) ceil(max(compPut(filter,1)))]);
ylim([floor(min(min(compPut(filter,[2,3])))) ceil(max(max(compPut(filter,[2,3]))))]);
xlabel('Strike Price');
ylabel('Price');
%legend('Empiric Prices', 'Numeric Prices');
title('Put Prices');

mtit(h, strcat(assetType,'-', int2str(k), '-', monthCode));
set(h, 'Visible', 'off'); % to solve mtit visibility 'bug'

saveas(h, strcat('validation/', assetType,'-', int2str(k), '-', monthCode, '.eps'), 'epsc');

toc

% Calculate Volatility, Skewness and Kurtosis.
callData = callData(callData(:,selectPrice)>0, [1, selectPrice]);
putData  = putData(putData(:,selectPrice)>0, [1, selectPrice]);
%[volatility, skewness, kurtosis] = ...
%    getImpliedMoments(underlyingPrice, R, daysForExpiration/250, callData, putData);

index = find(callData(:, 2) < underlyingPrice); % get call data ATM.
index = index(end);
volatility = blsimpv(underlyingPrice, callData(index, 1), R, ... 
    daysForExpiration/250, callData(index, 2), 10,  0,  1e-6, true);

% less data is required to show the whole pdf in SVG on the website
pdf = pdf(~isnan(pdf(:,2)), :);
pdf = pdf(pdf(:,2) > 10e-5, :);
pdf = pdf(1:3:length(pdf),:);

% Save JSON file
indicatorData = struct('currentDate', datestr(date,'dd-mmm-yyyy'), ...
              'expirationDate', expirationDate, ...
              'underlyingPrice', underlyingPrice, ...
              'expectedPrice', expectedValue, ... % wrong, not symetrical
              'daysForExpiration', daysForExpiration, ...
              'riskFreeRate', R, ...
              'volatility', volatility);
tables = struct('tableIncrease_bushels', tableIncrease_bushels, ...
               'tableDrop_bushels', tableDrop_bushels, ...
               'tableIncrease_metricTon', tableIncrease_metricTon, ...
               'tableDrop_metricTon', tableDrop_metricTon);
jsonContainer = struct('data', indicatorData, ...
                       'pdf', pdf, ...
                       'cut_points', cut_points, ...
                       'tables', tables);
finalJson = savejson('indicator', jsonContainer);

% remove all line-breaks
finalJson = regexprep(finalJson, '/\s+/S', '');
finalJson = strrep(finalJson, sprintf('\t'), '');
finalJson = strrep(finalJson, sprintf('\n'), '');
finalJson = strrep(finalJson, ' ', '');

if exist('storeRows','var')
    newRow = cell(1,4);
    newRow{1,1} = datestr(now,'mm-dd-yyyy'); % format for google spreadsheets.
    newRow{1,2} = datestr(expirationDate,'dd-mmm-yyyy');
    newRow{1,3} = underlyingPrice;
    newRow{1,4} = monthCode;
    newRow{1,5} = assetType;
    newRow{1,6} = finalJson;
    newRow{1,7} = optionJson;
    storeRows = vertcat(storeRows, newRow);
else
    storeRows = cell(1,4);
    storeRows{1,1} = datestr(now,'mm-dd-yyyy');
    storeRows{1,2} = datestr(expirationDate,'dd-mmm-yyyy');
    storeRows{1,3} = underlyingPrice;
    storeRows{1,4} = monthCode;
    storeRows{1,5} = assetType;
    storeRows{1,6} = finalJson;
    storeRows{1,7} = optionJson;
end

clearvars -except assetType selectUnit selectOption selectPrice R ...
    blendThreshold leftDif rightDif k storeRows ...
    dataVerification priceCheck priceRange_bushels priceRange_metricTon ...
    dS x_limit unit convertToMetricTon fitDegree monthCodes

end % end option looping

end % end asset looping

%end % end unit looping

cell2csv('loadData.csv', storeRows, ';'); % csv must be split then by ;

% save option history
append = fileread('loadData.csv');
history = fopen('history.csv', 'a');
fprintf(history, append);
fclose(history);
