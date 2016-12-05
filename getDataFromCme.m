% Get Data From CME
% getDataFromCme.m (4/03/2015)
%
% This program extracts option data from the CME.
% It is really easy to modify this program to work 
% with other type of options.
%
% Units: Cents per bushel.
%
% Programmed by Martin Baigorria

function optionData = getDataFromCme(asset, select_option)

    TIMEOUT = 20;
    addpath('jsonlab');

    file_path = strcat('cache/', datestr(now, 1), asset, num2str(select_option), '.mat');
    
    if exist(file_path, 'file')
        load(file_path);
        display(sprintf('Getting data for %s', optionData.expirationMonth));
        display('Found data in cache!');
        return
    end

    % input validation
    %if nargin <= 1
    %    error('Format: getDataFromCme(asset, select_option)');
    %end
    
    % select_option is an integer used to select the option month to 
    % load the JSON string from.
    if ~isnumeric(select_option)
        error('Error 1: Invalid parameter: Not numeric.');
    end
    
    asset = lower(asset);
    
    [product_group, product_subgroup] = getAssetData(asset);
    
    % Get page html from cme
    try
        url = strcat('http://www.cmegroup.com/trading/', ...
            product_group,'/',product_subgroup,'/',asset,'_quotes_globex_options.html');
        getHtml = urlread(url, 'Timeout', TIMEOUT);
    catch
        disp(url);
        error('Error 2: CME seems to be down, failed to retrieve data!');
    end
    
    if strfind(getHtml, 'experiencing issues')
         error('Error 3: We are currently experiencing issues with our weekly options quotes showing stale data. We are working to have this fixed as soon as possible.');
    end

    % In order to request the JSON price data, we must first build a 
    % sintactically corrrect request URL. This is built from the JSON
    % data in the original HTML.
    from = 'component.categories = ';
    to   = ';';
    pos1 = strfind(getHtml, from);
    pos2 = pos1 + strfind(getHtml(pos1:end), to);
    jsObject = getHtml((pos1 + length(from)):(pos2(1) - length(to) - 1));

    % convert js object into JSON (remove format first)
    jsObject = sscanf(jsObject,'%s');
    jsObject = strrep(jsObject, ',}}', '}}');
    json = strrep(jsObject, '},}', '}}');

    % the site sometimes behaves badly. changes in the html or js object ruin
    % the parsing.
    try
        rawStruct = loadjson(json);
    catch exception
        error('Error 4: Unable to load option data from JSON.');
    end

    % To unify the criteria to get the data from CME, we have an
    % issue with field names. This is from jsonlab toolbox:
    % From MATLAB doc: field names must begin with a letter, which may be
    % followed by any combination of letters, digits, and underscores.
    % Invalid characters will be converted to underscores, and the prefix
    % "x0x[Hex code]_" will be added if the first character is not a letter.
    %
    % This is used to get the correct data from the option dropdown.
    
    sections = fieldnames(rawStruct);
	path = rawStruct.(sections{1}).expirations;
    
    % get data from the new struct.
    optionNames = fieldnames(path);
    availableOptions = cell(length(optionNames),2);
    optionProductId = path.(optionNames{1}).productId; % from cmegroup

    for i = 1:length(fieldnames(path));
        availableOptions{i,1} = cellstr(path.(optionNames{i}).label);
        availableOptions{i,2} = cellstr(path.(optionNames{i}).expiration);
    end
    
    if select_option > length(fieldnames(path))
        error('Error 5: Non-existent contract selected.');
    end
    
    % Now that we have the data to build the URL, we must pick a specific
    % month to do the request.
    monthCode = char(availableOptions{select_option,2});
    
    % Get price data from the available options via id.
    % URL Format
    % 1. Base url: http://www.cmegroup.com/CmeWS/mvc/Quotes/Option/
    % 2. optionProductId
    % 3. G: Grains
    % 4. Future month code.
    % 5. Option product id.
   % jsonString = urlread(strcat('http://www.cmegroup.com/CmeWS/mvc/Quotes/Option/321/G/',char(availableOptions{select_option,2}),'/ALL?optionProductId=',num2str(optionProductId)));
   
    %display(strcat('http://www.cmegroup.com/CmeWS/mvc/Quotes/Option/', num2str(optionProductId),'/G/',char(availableOptions{select_option,2}),'/ALL'));
   
    try
        url = strcat('http://www.cmegroup.com/CmeWS/mvc/Quotes/Option/', ...
            num2str(optionProductId),'/G/',char(availableOptions{select_option,2}), '/ALL');
        jsonString = urlread(url, 'Timeout', TIMEOUT);
    catch
        disp(url);
        error('Error 6: CME seems to be down, failed to retrieve data!');
    end
    
    rawStruct = loadjson(jsonString);

    underlyingFuture = rawStruct.underlyingFutureContractQuotes{1, 1}; % rename rawStruct for the underlying.

    numberOfStrikes = length(rawStruct.optionContractQuotes);

    % pre-allocate memory for call and put data.
    callData = zeros(numberOfStrikes, 6); %format: Strike Price | Best Bid | Best Offer | Prior Settle | Volume
    putData  = zeros(numberOfStrikes, 6);

    % this is the contract identifier, from the expiration of the
    % underlying.
    expiration = underlyingFuture.expirationMonth;
    expiration = lower(expiration);
    expiration(1,1) = upper(expiration(1,1));
    
    % get contact expiration date (which differs from underlying expiration
    % date). since i'm not extracting any JSON (not available), this might
    % not be 100% failsafe. i dont think the CME will change this table, so
    % to make it simple (not to have to convert from JS object format into
    % JSON) i will just process it straight from the HTML table.
    try
        url = strcat('http://www.cmegroup.com/trading/', ...
            product_group,'/', product_subgroup,'/',asset,'_product_calendar_options.html?optionProductId=', num2str(optionProductId));
        expirationDate = urlread(url, 'Timeout', TIMEOUT);
    catch
        disp(url);
        error('Error 7: CME seems to be down, failed to retrieve data!');
    end
    
    from = strcat('<th scope="row">', expiration,'</td>');
    to   = '</tr>';
    pos1 = strfind(expirationDate, from);
    pos2 = pos1 + strfind(expirationDate(pos1:end), to);
    
	expirationDate = expirationDate((pos1 + length(from)):(pos2(1) - length(to)));
	expirationDate = strrep(expirationDate, '</td>', ''); % can be avoided
    expirationDate = expirationDate(length(expirationDate)-10:length(expirationDate));
    expirationDate = strrep(expirationDate, ' ', '-');
    
    daysForExpiration = wrkdydif(date, expirationDate, length(holidays(date, expirationDate)));

    if daysForExpiration < 0
        error('Error: You have loaded an expired contract.');
    end;
    
    display(sprintf('Getting data for %s', expiration));
    
    % parse underlying price
    underlyingPrice = parsePrice(underlyingFuture.last);
    underlyingPriorSettle = parsePrice(underlyingFuture.priorSettle);

    % process raw data
    for i = 1:numberOfStrikes

        % strike price
        callData(i,1) = str2num(rawStruct.optionContractQuotes{1, i}.strikePrice);
        putData(i,1)  = str2num(rawStruct.optionContractQuotes{1, i}.strikePrice);

        % best bid
        callData(i,2) = parsePrice(rawStruct.optionContractQuotes{1, i}.call.low);
        putData(i,2)  = parsePrice(rawStruct.optionContractQuotes{1, i}.put.low);

        % best offer
        callData(i,3) = parsePrice(rawStruct.optionContractQuotes{1, i}.call.high);
        putData(i,3)  = parsePrice(rawStruct.optionContractQuotes{1, i}.put.high);

        % prior settle
        callData(i,5) = parsePrice(rawStruct.optionContractQuotes{1, i}.call.priorSettle);
        putData(i,5)  = parsePrice(rawStruct.optionContractQuotes{1, i}.put.priorSettle);
        
        % volume
        callData(i,6) =  str2num(strrep(rawStruct.optionContractQuotes{1, i}.call.volume, ',', ''));
        putData(i,6)  =  str2num(strrep(rawStruct.optionContractQuotes{1, i}.put.volume, ',', ''));

    end
    
    % average price
    callData(:,4) = (callData(:,2) + callData(:,3))/2;
    putData(:,4)  = (putData(:,2) + putData(:,3))/2;
    
    %if sum(~isnan(putData(:,3))) < 5
    %    error('Insuficient amount of data. Probably an issue with CME.');
    %end
    
    % filter zero values (CME bug?)
    callData = callData(callData(:,5) > 0, :);
    putData = putData(putData(:,5) > 0, :);
    
    optionData = struct('type', product_group, ...
                        'monthCode', monthCode, ...
                        'expirationMonth', expiration, ...
                        'expirationDate', expirationDate, ...
                        'daysForExpiration', daysForExpiration, ...
                        'underlyingPrice', underlyingPrice, ...
                        'underlyingPriorSettle', underlyingPriorSettle, ...
                        'callData', callData, ...
                        'putData', putData);
                    
    save(file_path, 'optionData');
                    
end

function [product_group, product_subgroup] = getAssetData(asset)

    if strcmp(asset, 'soybean') || strcmp(asset, 'corn') || ...
       strcmp(asset, 'wheat')   || strcmp(asset, 'kc-wheat')
        product_group    = 'agricultural';
        product_subgroup = 'grain-and-oilseed';
    else
        error('Error 8: Invalid asset.');
    end

end

function price = parsePrice(rawPrice)

	price = strrep(rawPrice, 'a', '');
    price = strrep(price, 'b', '');
    price = strrep(price, char(39), '.');
    price = str2double(price);
    
end