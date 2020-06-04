classdef readimzML < customisePlot
    %READIMZML uses the Java classes from the imzML Converter [1] to parse 
    % imzML files into MATLAB
    %
    % [1] Alan M. Race, Iain B. Styles, Josephine Bunch, Journal of Proteomics, 75(16):5111-5112, 2012. http://dx.doi.org/10.1016/j.jprot.2012.05.035 
    
    properties
        pixelRows
        pixelColumns
        totIonCount
        spectra
        filePath
        files
        imzMLConverterFolder
    end
    
    methods 
        function obj = locateFiles(obj)
            if isempty(obj.imzMLConverterFolder)
                [fileName,pathName] = uigetfile('.jar',...
                    'Locate imzMLConverter.jar');
                if isequal(fileName,0)
                    return
                end
                obj.imzMLConverterFolder = fullfile(pathName,fileName);
                javaclasspath(obj.imzMLConverterFolder);
            end
           [fileName,pathName] = uigetfile('.imzML',...
               'Select .imzML files',...
               'MultiSelect','on');
           if isequal(fileName,0)
               disp('No .imzML files selected');
               return
           end
           obj.filePath = pathName;
           obj.files = fileName;
        end
        
        function obj = parseFiles(obj)
            clc
            validateIMZML(obj);
            findIBD(obj);
            if isempty(obj.files)
               disp('No files selected');
            end
            
            if ~iscell(obj.files)
               fileCount = 1; 
            else
               fileCount = numel(obj.files); 
            end
            
            for j = 1:fileCount
               
               if fileCount == 1
                  fprintf('Parsing file: %s \n',obj.files); 
                  imzML = imzMLConverter.ImzMLHandler.parseimzML(fullfile(obj.filePath,obj.files)); 
               else
                  fprintf('Parsing file: %s \n',obj.files{j});  
                  imzML = imzMLConverter.ImzMLHandler.parseimzML(fullfile(obj.filePath,obj.files{j}));  
               end
               
               if j == 1
                  obj.pixelRows = imzML.getHeight();
                  obj.pixelColumns = imzML.getWidth();
               end
               
               count = 0; tempSpectra = []; tempTIC = [];
               for c = 1:obj.pixelRows 
                    for r = 1:obj.pixelColumns
                        count = count+1;
                        mzs = imzML.getSpectrum(c,r).getmzArray();
                        int = imzML.getSpectrum(c,r).getIntensityArray();
                        tempSpectra{count,1} = [mzs,int];
                        tempTIC{count,1} = sum(int);
                    end
               end
               obj.totIonCount{j} = tempTIC;
               obj.spectra{j} = tempSpectra;
            end
            fprintf('Parsing completed \n');
        end
    end
    
    methods
        function validateIMZML(obj) 
            fileLoc = obj.files;
            if ~iscell(fileLoc)
                tf = ~isempty(regexpi(fileLoc,'imzml','match'));
                if tf == true
                    fprintf('File is .imzML \n');
                else
                    error('File is not .imzML');
                end
            else
                for j = 1:length(fileLoc)
                    tf = ~isempty(regexpi(fileLoc{j},'imzml','match'));
                    if tf == true
                        fprintf('File %d is .imzML \n',j);
                    else
                        error(sprintf('File % d is not .imzML',j));
                    end
                end
            end
        end
        
        function findIBD(obj)
            fileLoc = obj.files;
            if ~iscell(fileLoc)
               fileExtensionLoc = find(fileLoc=='.');
               fileNameShort = fileLoc(1:fileExtensionLoc-1);
               ibdName = [fileNameShort '.ibd'];
               if exist(fullfile(obj.filePath,ibdName),'file') == 2
                  fprintf('.ibd file found \n'); 
               else
                  error('.ibd file not found'); 
               end
            else
                for j = 1:length(fileLoc)
                   fileExtensionLoc = find(fileLoc{j}=='.');
                   fileNameShort = fileLoc{j}(1:fileExtensionLoc-1);
                   ibdName = [fileNameShort '.ibd'];
                   if exist(fullfile(obj.filePath,ibdName),'file') == 2
                      fprintf('.ibd file found \n'); 
                   else
                      error('.ibd file not found'); 
                   end 
                end
            end
        end
        
        function checkImage(obj)
            fileToCheck = input('Select file: ');
            try
               fileTIC = cell2mat(obj.totIonCount{fileToCheck}); 
            catch e
               rethrow(e)
            end
            reconstructedImage = reshape(fileTIC,obj.pixelRows,obj.pixelColumns);
            colormap(obj.CMAP)
            imagesc(reconstructedImage);
            xlabel(obj.XLabel);
            ylabel(obj.YLabel);
            colorbar();
        end
    end
    
end

