%% RMS CONTRAST IMAGE ANALYSIS 

function rmsContrast = computeRMSContrast(image)

    image = double(image);

    pixels = image(:);

    mu = mean(pixels);

    rmsContrast = sqrt(mean((pixels - mu).^2));

end 