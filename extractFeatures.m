classdef extractFeatures < readimzML
    % Extracting and aligning relevant features from SIMS spectra 

    properties
        featureList     
        uniqueFeatures
    end

    properties (Access = private)        
        peakList
    end

    methods
        function obj = getFeatureList(obj)
           clc 
           obj.uniqueFeatures = []; obj.featureList = [];
           validateIntensityInput(obj); 
           fprintf('Performing peak picking \n');
           if ~iscell(obj.files)
               tempMS = []; tempPeaks = []; tempList = [];
               tempList = obj.spectra{1};
               for j = 1:length(tempList)
                   tempMS = cell2mat(tempList(j,1));
                   tempPeaks{j,1} = mspeaks(tempMS(:,1),tempMS(:,2),...
                   'HeightFilter',obj.options.thresholdIntensity,...
                   'Denoising',false);
               end  
               obj.peakList = tempPeaks;
               fprintf('Peak picking completed \n');
               fprintf('Finding unique features \n');
               obj = getUniqueFeatures(obj);
               fprintf('Aligning features \n');
               obj = alignFeatures(obj);
           else
               tempPeakList = [];
               for j = 1:length(obj.files)
                   fprintf('File: %d \n',j);
                   tempMS = []; tempPeaks = []; tempList = []; 
                   tempList = obj.spectra{j};
                   for n = 1:length(tempList)
                       tempMS = cell2mat(tempList(n,1));
                       tempPeaks{n,1} = mspeaks(tempMS(:,1),tempMS(:,2),...
                       'HeightFilter',obj.options.thresholdIntensity,...
                       'Denoising',false);
                   end
                   obj.peakList = tempPeaks;
                   obj = getUniqueFeatures(obj,j);
                   obj = alignFeatures(obj,j);
                   fprintf('Peak picking completed \n');
               end            
           end
           fprintf('Done \n');
        end

        function obj = validateIntensityInput(obj)
            if isequal(obj.options.thresholdType,'basepeak')
                if obj.options.thresholdIntensity <= 0 || obj.options.thresholdIntensity >= 100
                   error('Threshold intensity for type "basepeak" should be between 0 and 100'); 
                end
            elseif isequal(obj.options.thresholdType,'absolute')
                if obj.options.thresholdIntensity < 0
                   error('Threshold intensity for type "absolute" should be 0 or larger'); 
                end
            elseif isequal(obj.options.thresholdType,'snr')
                if obj.options.thresholdIntensity <= 0
                   error('Threshold intensity for type "snr" should be larger than 0'); 
                end
            else
               error('Invalid input for threshold'); 
            end
            fprintf('Threshold type:    %s \n',obj.options.thresholdType);
            fprintf('Threshold value:   %s \n',num2str(obj.options.thresholdIntensity));
        end

       function obj = getUniqueFeatures(obj,iteration)
            peakVector = cell2mat(obj.peakList);
            tolerance = 0.001;
            if nargin < 2
               iteration = 1; 
            end
            r = [];
            tempPeaks = [];
            for j = 1:length(peakVector) 
              if ~isempty(r) 
                 if ~isempty(find(r(:,1)==j,1)) 
                    continue
                 end
              end

              matchIons = find(peakVector(:,1) >= peakVector(j,1)-tolerance & ... 
                                peakVector(:,1) <= peakVector(j,1)+tolerance);
              if numel(matchIons) > 1
                  r = [r;matchIons];
                  tempPeaks = [tempPeaks;median(peakVector(matchIons,1))];
              else 
                  r = [r;matchIons];
                  tempPeaks = [tempPeaks;peakVector(matchIons,1)];
              end
            end
            obj.uniqueFeatures{iteration} = tempPeaks;
       end

       function obj = alignFeatures(obj,iteration)
           if nargin < 2
               iteration = 1;
           end
           tolerance = 0.001;
           featureMatrix = zeros(length(obj.uniqueFeatures),obj.pixelRows*obj.pixelColumns);
           tempFeatures = cell2mat(obj.uniqueFeatures(iteration));
            for k = 1:length(tempFeatures)
                index = cellfun(@(x) find(x(:,1) > tempFeatures(k,1)-tolerance & x(:,1) < tempFeatures(k,1)+tolerance),...
                    obj.peakList, 'UniformOutput', false);
                emptyIndex = cellfun(@(x) isempty(x), index, 'UniformOutput', true);
                for p = 1:length(emptyIndex)
                    if emptyIndex(p,1) == true
                        continue
                    else
                        tempIDX = index{p,1};
                        if numel(tempIDX) > 1
                            tempPeaks = cell2mat(obj.peakList(p,1));
                            tempMZ = tempPeaks(tempIDX,1);
                            difference = tempFeatures(k,1)-tempMZ;
                            minDiff = find(difference==min(difference));
                            featureMatrix(k,p) = tempPeaks(tempIDX(minDiff),2);
                        else
                            tempPeaks = cell2mat(obj.peakList(p,1));
                            featureMatrix(k,p) = tempPeaks(tempIDX,2);
                        end
                    end
                end
            end
            obj.featureList{iteration} = featureMatrix;
       end
    end

end