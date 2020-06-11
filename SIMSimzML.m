classdef SIMSimzML < readimzML & customisePlot
    % Parsing .imzML files exported from Surface Lab for generating ion
    % images on mass spectrometry imaging data
    %
    % Written in MATLAB R2017a 
    %
    % Required toolboxes:
    % - Bioinformatics Toolbox
    % - Statistics & Machine Learning Toolbox
    
    properties 
        version = '0.0.6'
        developer = 'Joris Meurs, MSc'
        matlabVersion = 'R2017a'
        dependencies = {'Bioinformatics Toolbox','Statistics & Machine Learning Toolbox'}
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
               warning('No m/z value selected. Executing selectedMZ');
               obj = selectMZ(obj);
           end
           try
               if isequal(obj.options.plotimages,'single')
                   spectralData = obj.spectra{obj.file};
                   fileTIC = cell2mat(obj.totIonCount{obj.file});
                   mzInt = [];
                   for j = 1:length(spectralData)
                      pixelMS = cell2mat(spectralData(j,1)); 
                      ionIDX = find(pixelMS(:,1) > obj.mz-0.001 & pixelMS(:,1) < obj.mz+0.001);
                      if ~isempty(ionIDX)
                         mzInt = [mzInt;sum(pixelMS(ionIDX,2))/fileTIC(j,1)]; 
                      else
                         mzInt = [mzInt;0]; 
                      end
                   end
                   reconstructedIntensities = reshape(mzInt,...
                       obj.pixelRows,obj.pixelColumns);
                   colormap(obj.CMAP);
                   imagesc(reconstructedIntensities);
                   xlabel(obj.XLabel);
                   ylabel(obj.YLabel);
                   colorbar();
                   if ~isempty(obj.XTick)
                        set(gca,'XTick',XTick);
                   end
                   if ~isempty(obj.YTick)
                        set(gca,'YTick',YTick);
                   end 
                   if ~isempty(obj.XTickLabels)
                        set(gca,'XTickLabels',num2str(XTickLabels)) 
                   end
                   if ~isempty(obj.YTickLabels)
                        set(gca,'XTickLabels',num2str(XTickLabels))
                   end
               end
               if isequal(obj.options.plotimages,'all')
                   for j = 1:length(obj.spectra)
                      spectralData = obj.spectra{j};
                      fileTIC = cell2mat(obj.totIonCount{j});
                      mzInt = [];
                      for n = 1:length(spectralData)
                          pixelMS = cell2mat(spectralData(n,1));
                          ionIDX = find(pixelMS(:,1) > obj.mz-0.001 & pixelMS(:,1) < obj.mz+0.001);
                          if ~isempty(ionIDX)
                             mzInt = [mzInt;sum(pixelMS(ionIDX,2))/fileTIC(j,1)]; 
                          else
                             mzInt = [mzInt;0]; 
                          end
                      end
                      figure(j)
                      reconstructedIntensities = reshape(mzInt,...
                      obj.pixelRows,obj.pixelColumns);
                      colormap(obj.CMAP);
                      imagesc(reconstructedIntensities);
                      xlabel(obj.XLabel);
                      ylabel(obj.YLabel);
                      title(obj.files{j},'interpreter','none');
                      colorbar();
                      if ~isempty(obj.XTick)
                          set(gca,'XTick',XTick);
                      end
                      if ~isempty(obj.YTick)
                          set(gca,'YTick',YTick);
                      end    
                      if ~isempty(obj.XTickLabels)
                         set(gca,'XTickLabels',num2str(XTickLabels)) 
                      end
                      if ~isempty(obj.YTickLabels)
                          set(gca,'XTickLabels',num2str(XTickLabels))
                      end
                   end
               end
           catch e
               rethrow(e)
           end
        end
        
    end
    
end

