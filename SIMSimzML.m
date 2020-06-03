classdef SIMSimzML < readimzML & extractFeatures
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
        version = '0.0.2'
        developer = 'Joris Meurs, MSc'
        matlabVersion = 'R2017a'
        dependencies = {'Bioinformatics Toolbox','Statistics & Machine Learning Toolbox'}
    end
    
    properties
       mvaType = 'NMF'
       normalisedIntensity
       components = 2;
       selectedMZ
    end
    
    methods
        function obj = normalise(obj)
            normInt = [];
            if iscell(obj.totIonCount)
                obj.totIonCount = cell2mat(obj.totIonCount);
            end
            
            for j = 1:size(obj.featureList,2)
                normInt(:,j) = obj.featureList(:,j)./obj.totIonCount(j,1);
            end
            obj.normalisedIntensity = normInt'; 
        end
        
        function obj = selectMZ(obj)
            clc
            for j = 1:length(obj.uniqueFeatures)
                fprintf('(%d) m/z %.4f \n',j,obj.uniqueFeatures(j));
            end
            obj.selectedMZ = input('Select row number: ');
        end
        
        function obj = ionImages(obj)
           if isempty(obj.selectedMZ)
               warning('No m/z value selected. ');
           end
        end
        
        function obj = multivariate(obj)
            if isempty(obj.normalisedIntensity)
                warning('Normalising intensities first');
                obj = normalise(obj);
            end

            rng(1);
            [W,H] = nnmf(obj.normalisedIntensity,obj.components);
            for j = 1:obj.components
                f = figure;
                endmemberImage = reshape(W(:,j),obj.pixelRows,obj.pixelColumns);
                colormap(hot);
                imagesc(endmemberImage)
                colorbar();
                title(sprintf('Endmember %d',j));
                xlabel('X');
                ylabel('Y');
                set(gcf,'Color','white');
            end
        end
    end
    
end

