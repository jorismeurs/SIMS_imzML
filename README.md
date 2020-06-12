# SIMS_imzML
Processing SIMS imzML files exported from SurfaceLab

## Main functionality
1. ```obj = SIMSimzML``` Initiate class with default values </br>
2. ```obj = obj.locateFiles``` Locate the imzML Converter .jar file ```imzMLConverter.jar``` followed by the .imzML files of interest </br>
3. ```obj = obj.parseFiles``` Obtain spectra per pixel per file along with the image dimensions. All spectra are stored in N X 1 cell array </br>
Next, the spectra can be further processed to generate the ion images

## Processing options
To generate an ion image, this class allows two methods. The default setting ```obj.options.featureSelection``` is ```'peaks'``` which 
performs peak picking first to select all unique features throughout all spectra. By default the threshold intensity for peaks is set
to 1000. This can be altered to e.g. 100 as follows ```obj.options.thresholdIntensity = 100```. Feature selection is then executed by
```obj = obj.getFeatureList```
</br>
Alternatively, features are directly taken from the profile spectra. First the *m/z* value of interest has to be set to e.g. 146.0610 
```obj.mz = 146.0610``` (recommended when peak intensities are low).
</br>
To generate the ion image(s), execute ```obj.ionImage```. Per file, an ion image is generated. When ```obj.options.featureSelection``` is  set to
```'peaks'```, the ion index has to be used as input. A number of potential ion will be shown in the command window which are found within the
```obj.options.tolerance``` window. The default export image format is ```.tif``` and can be altered via ```obj.options.exportimage```. By
default, images will saved. Type ```obj.options.saveimage = 'false'``` if export of images is not required 

## Customising plot
- Change label on x-axis:     ```obj.XLabel``` 
- Change label on y-axis:     ```obj.YLabel```
- Change tick steps x-axis:   ```obj.XTick```
- Change tick steps y-axis:   ```obj.YTick```
- Change tick labels x-axis:  ```obj.XTickLabels```
- Change tick labels y-axis:  ```obj.YTickLabels```
- Change font:                ```obj.fontName```
- Change colormap:            ```obj.CMAPval``` 
  - Options:                  ```hot``` or ```cividis```
  - Run ```obj = obj.setCMAP```
