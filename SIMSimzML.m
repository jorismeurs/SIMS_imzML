classdef SIMSimzML < readimzML & customisePlot & extractFeatures
    % Parsing .imzML files exported from Surface Lab for generating ion
    % images on mass spectrometry imaging data
    %
    % Written in MATLAB R2017a 
    %
    % Required toolboxes:
    % - Bioinformatics Toolbox
    % - Statistics & Machine Learning Toolbox
    
    properties 
        version = '0.2.2'
        developer = 'Joris Meurs, MSc'
        matlabVersion = 'R2017a'
        dependencies = {'Bioinformatics Toolbox'}
    end
    
    properties
       heterogeneityData
       options
       file
       mz
       intensityData
    end
    
    methods
        
        function obj = SIMSimzML
            obj.options.plotimages = 'all';
            obj.options.saveimage = 'true';
            obj.options.exportformat = 'tif';
            obj.options.tolerance = 0.001;
            obj.options.featureSelection = 'peaks';
            obj.options.thresholdIntensity = 1000;
            obj.options.title = 'false';
            obj.options.thresholdType = 'absolute';
        end
        
        function obj = selectFile(obj)
           clc
           for j = 1:length(obj.files)
               fprintf('(%d) %s \n',j,obj.files{j});
           end
           obj.file = input('Select file: ');
        end
        
        function obj = selectMZ(obj)
            clc
            obj.mz = input('Enter m/z value: ');
        end
        
        function obj = ionImage(obj)

           if isequal(obj.options.plotimages,'all')
               for j = 1:length(obj.files)
                  mzInt = constructImage(obj,j); 
                  obj.intensityData{j} = mzInt;
                  f = figure;
                  reconstructedIntensities = reshape(mzInt,...
                  obj.pixelRows,obj.pixelColumns);
                  colormap(obj.CMAP);
                  imagesc(reconstructedIntensities);
                  drawnow;
                  setPlot(obj,j)  
                  if isequal(obj.options.saveimage,'true')
                      currentFolder = cd;
                      exportFolder = [cd '\images\'];
                      if ~exist(exportFolder,'dir')
                         mkdir images
                      end
                      try
                        saveas(f,[currentFolder '\images\' obj.files{j} '.tif']);
                      catch
                        saveas(f,[currentFolder '\images\' obj.files '.tif']);  
                      end
                  end
                  close(f);
               end
           end
        end
        
        function [mzInt,obj] = constructImage(obj,iteration)
           if isempty(obj.mz)
               warning('No m/z value input');
               return
           end
           if isequal(obj.options.featureSelection,'peaks')
              mzInt = []; 
              clc
              pixelMS = cell2mat(obj.uniqueFeatures(iteration)); 
              potentialIDX = find(pixelMS(:,1) > obj.mz-obj.options.tolerance & pixelMS(:,1) < obj.mz+obj.options.tolerance);
              for k = 1:length(potentialIDX)
                 fprintf('(%d) m/z %.4f \n',potentialIDX(k),pixelMS(potentialIDX(k),1)); 
              end
              ionIDX = input('Select m/z value index: ');
              tempMat = cell2mat(obj.featureList(iteration)); 
              tempTIC = cell2mat(obj.totIonCount{iteration});
              mzInt = tempMat(ionIDX,:);
              mzInt = mzInt'./tempTIC;
              mzInt = mzInt';
           end
           if isequal(obj.options.featureSelection,'profile')
               fileSpectra = obj.spectra{iteration};
               tempTIC = cell2mat(obj.totIonCount{iteration});
               mzInt = zeros(length(fileSpectra),1);
               for j = 1:length(fileSpectra)
                   pixelMS = cell2mat(fileSpectra(j,1));
                   ionIDX = find(pixelMS(:,1) > obj.mz-obj.options.tolerance & ...
                       pixelMS(:,1) < obj.mz+obj.options.tolerance);
                   if ~isempty(ionIDX)
                       mzInt(j,1) = sum(pixelMS(ionIDX,2));
                   end
               end
               mzInt = mzInt./tempTIC;
               %disp(mean(mzInt));
               mzInt = mzInt';
           end
        end
        
        function setPlot(obj,j)
              if isequal(obj.options.title,'true')
                  try
                      title(obj.files{j},'interpreter','none'); 
                  catch
                      title(obj.files,'interpreter','none'); 
                  end
              end
              xlabel(obj.XLabel);
              ylabel(obj.YLabel);          
              colorbar();
              if ~isempty(obj.XTick)
                  set(gca,'XTick',obj.XTick);
              end
              if ~isempty(obj.YTick)
                  set(gca,'YTick',obj.YTick);
              end    
              if ~isempty(obj.XTickLabels)
                 set(gca,'XTickLabels',obj.XTickLabels) 
              end
              if ~isempty(obj.YTickLabels)
                  set(gca,'YTickLabels',obj.YTickLabels)
              end
              if ~isempty(obj.colorbarMin) && ~isempty(obj.colorbarMax)
                 caxis(gca,[obj.colorbarMin obj.colorbarMax]) 
              end
              set(gcf,'Color','white');
              set(gca,'FontName',obj.fontName,...
                  'FontSize',obj.fontSize);
        end
        
        function obj = ticImage(obj)
              for j = 1:length(obj.totIonCount)
                  tempTIC = cell2mat(obj.totIonCount{j});
                  tempTIC = reshape(tempTIC,obj.pixelRows,obj.pixelColumns);
                  f = figure;
                  colormap(obj.CMAP);
                  imagesc(tempTIC);
                  setPlot(obj,j);
                  if isequal(obj.options.saveimage,'true')
                      currentFolder = cd;
                      exportFolder = [cd '\images\'];
                      if ~exist(exportFolder,'dir')
                         mkdir images
                      end
                      try
                        saveas(f,[currentFolder '\images\TIC_' obj.files{j} '.tif']);
                      catch
                        saveas(f,[currentFolder '\images\TIC_' obj.files '.tif']);  
                      end
                  end
                  close(f);
              end
        end
        
        function obj = imageHeterogeneity(obj)
           obj.heterogeneityData = []; 
           for n = 1:length(obj.spectra) 
               fileSpectra = obj.spectra{n};
               tempTIC = cell2mat(obj.totIonCount{n});
               mzInt = zeros(length(fileSpectra),1);
               for j = 1:length(fileSpectra)
                   pixelMS = cell2mat(fileSpectra(j,1));
                   ionIDX = find(pixelMS(:,1) > obj.mz-obj.options.tolerance & ...
                       pixelMS(:,1) < obj.mz+obj.options.tolerance);
                   if ~isempty(ionIDX)
                       mzInt(j,1) = sum(pixelMS(ionIDX,2));
                   end
               end
               mzInt = mzInt./tempTIC;
               %mzInt(mzInt==0) = NaN;
               obj.heterogeneityData = [obj.heterogeneityData;ones(length(mzInt),1).*(21+n),mzInt];
           end
        end
    end
    
end

