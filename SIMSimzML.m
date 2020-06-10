classdef SIMSimzML < readimzML & customisePlot
    % Parsing .imzML files exported from Surface Lab for generating ion
    % images and performing multivariate analysis on mass spectrometry
    % imaging data
    %
    %
    % Written in MATLAB R2017a 
    %
    % Required toolboxes:
    % - Bioinformatics Toolbox
    % - Statistics & Machine Learning Toolbox
    
    properties 
        version = '0.0.5'
        developer = 'Joris Meurs, MSc'
        matlabVersion = 'R2017a'
        dependencies = {'Bioinformatics Toolbox','Statistics & Machine Learning Toolbox'}
    end
    
    properties
       components = 2;
       file
       mz
    end
    
    methods
        
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
           catch e
               rethrow(e)
           end
        end
        
    end
    
end

