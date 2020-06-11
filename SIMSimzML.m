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
        version = '0.0.7'
        developer = 'Joris Meurs, MSc'
        matlabVersion = 'R2017a'
        dependencies = {'Bioinformatics Toolbox'}
    end
    
    properties
       options
       file
       mz
    end
    
    methods
        
        function obj = SIMSimzML
            obj.options.plotimages = 'all';
            obj.options.saveimage = 'true';
            obj.options.exportformat = 'tif';
            obj.options.tolerance = 0.001;
            obj.options.featureSelection = 'peaks';
            obj.options.thresholdIntensity = 1000;
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
           clc
           if isempty(obj.mz)
               warning('No m/z value selected');
               return
           end

           if isequal(obj.options.plotimages,'all')
               for j = 1:length(obj.spectra)
                  spectralData = obj.spectra{j};
                  fileTIC = cell2mat(obj.totIonCount{j});
                  mzInt = constructImage(obj,spectralData,fileTIC);                
                  figure(j)
                  reconstructedIntensities = reshape(mzInt,...
                  obj.pixelRows,obj.pixelColumns);
                  colormap(obj.CMAP);
                  imagesc(reconstructedIntensities);
                  drawnow;
                  setPlot(obj,j)  
                  if isequal(obj.options.saveimage,'true')
                      
                  end
               end
           end
        end
        
        function mzInt = constructImage(obj)

           if isequal(obj.options.featureSelection,'peaks')
              mzInt = []; 
              for n = 1:length(obj.uniqueFeatures)
                 clc
                 pixelMS = cell2mat(obj.uniqueFeatures(n)); 
                 for k = 1:length(pixelMS)
                    fprintf('(%d) m/z %.4f \n',k,pixelMS(k,1)); 
                 end
                 ionIDX = input('Select m/z value index: ');
                 tempMat = cell2mat(obj.featureList(n)); 
                 tempTIC = cell2mat(obj.totIonCount{1});
                 mzInt = tempMat(ionIDX,:);
                 mzInt = mzInt'./tempTIC;
                 mzInt = mzInt';
              end
           end
        end
        
        function setPlot(obj,j)
              title(obj.files{j},'interpreter','none'); 
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
              set(gcf,'Color','white');
              set(gca,'FontName',obj.fontName);
        end
    end
    
end

