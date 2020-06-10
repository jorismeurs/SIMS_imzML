classdef SIMSimzML < readimzML & extractFeatures & customisePlot
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
       mvaType = 'NMF'
       normalisedIntensity
       components = 2;
       selectedFile
       selectedMZ
    end
    
    methods
        
        function obj = selectFile(obj)
           clc
           for j = 1:length(obj.files)
               fprintf('(%d) %s \n',j,obj.files{j});
           end
           obj.selectedFile = input('Select file: ');
        end
        
        function obj = normalise(obj)
            clc
            obj.normalisedIntensity = [];            
            tempTIC = cell2mat(obj.totIonCount{obj.selectedFile});
            tempFeatures = obj.featureList{obj.selectedFile};
            for j = 1:size(tempFeatures,2)
                obj.normalisedIntensity(:,j) = tempFeatures(:,j)./tempTIC(j,1);
            end
        end
        
        function obj = selectMZ(obj)
            clc
            tempFeatures = obj.uniqueFeatures{obj.selectedFile};
            for j = 1:length(tempFeatures)
                fprintf('(%d) m/z %.4f \n',j,tempFeatures(j));
            end
            obj.selectedMZ = input('Select row number: ');
        end
        
        function obj = ionImage(obj)
           clc
           if isempty(obj.selectedMZ)
               warning('No m/z value selected. Executing selectedMZ');
               obj = selectMZ(obj);
           end
           try
               reconstructedIntensities = reshape(obj.normalisedIntensity(obj.selectedMZ,:),...
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
        
        function obj = multivariate(obj)
            if isempty(obj.normalisedIntensity)
                warning('Normalising intensities first');
                obj = normalise(obj);
            end
            if isequal(obj.mvaType,'NMF')
                rng(1);
                [W,H] = nnmf(obj.normalisedIntensity',obj.components);
                for j = 1:obj.components
                    figure;
                    endmemberImage = reshape(W(:,j),obj.pixelRows,obj.pixelColumns);                   
                    subplot(1,2,1)
                    colormap(obj.CMAP);
                    imagesc(endmemberImage)
                    colorbar();
                    title(sprintf('Endmember %d',j));
                    xlabel(obj.XLabel);
                    ylabel(obj.XLabel);
                    axis square
                    set(gca,'FontName',obj.fontName);
                    
                    subplot(1,2,2)
                    stem(obj.uniqueFeatures{obj.selectedFile},H(j,:),'Marker','none','Color','b');
                    xlabel('m/z');
                    ylabel('H');
                    set(gca,'FontName',obj.fontName);
                    
                    set(gcf,'Color','white');
                end
            elseif isequal(obj.mvaType,'PCA')
                [loadings,scores,~,~,explained] = pca(obj.normalisedIntensity');
                
                for j = 1:obj.components
                    figure;
                    subplot(1,2,1) % Scores
                    colormap(obj.CMAP);
                    pcaImage = reshape(scores(:,j),obj.pixelRows,obj.pixelColumns);
                    imagesc(pcaImage);
                    xlabel(obj.XLabel);
                    ylabel(obj.YLabel);
                    title(sprintf('PC%d (%.1f%%)',j,explained(j)));
                    set(gca,'FontName',obj.fontName);
                    axis square
                                        
                    subplot(1,2,2) % Loadings
                    stem(obj.uniqueFeatures{obj.selectedFile},loadings(:,j),'Marker','none','Color','b');
                    xlabel('m/z');
                    ylabel(sprintf('Loadings PC%d',j));
                    set(gca,'FontName',obj.fontName);
                    set(gcf,'Color','white');
                end
                
            else
               error(sprintf('MVA type %s not supported',obj.mvaType)); 
            end
        end
    end
    
end

