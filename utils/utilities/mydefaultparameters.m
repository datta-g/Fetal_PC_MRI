function [Param,FixedImage,MovingImage,Mask] = mydefaultparameters(N)

Param{1} = elxDefaultParameters('TranslationTransform', 2);
Param{1}.FixedImagePyramid = 'FixedSmoothingImagePyramid';
Param{1}.MovingImagePyramid = 'MovingSmoothingImagePyramid';
Param{1}.Metric = 'NormalizedMutualInformation';
Param{1}.AutomaticScalesEstimation = true;
Param{1}.AutomaticTransformInitialization = true;
Param{1}.AutomaticTransformInitializationMethod = 'CenterOfGravity';
Param{1}.NumberOfHistogramBins = 32;
Param{1}.UseRandomSampleRegion = true;
Param{1}.NewSamplesEveryIteration = true;
Param{1}.NumberOfResolutions = 2;
Param{1}.ImageSampler = 'RandomCoordinate';
Param{1}.MaximumNumberOfIterations = 300;
Param{1}.DefaultPixelValue = 0;
Param{1}.SampleRegionSize = [50 50];
FixedImage.x{1,1}=0:N-1;
FixedImage.x{1,2}=0:N-1;
MovingImage.x{1,1}=0:N-1;
MovingImage.x{1,2}=0:N-1;
Mask.x{1,1}=0:N-1;
Mask.x{1,2}=0:N-1;

end