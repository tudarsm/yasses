# YASSES
Yet another Siemens Star evaluation script (YASSES) is a code to determine the object side MTF curve from Siemens Star data.

# How it works
Instead of resampling the siemens star in concentric circles, the entire image is transformed into polar coordinates (R,Theta), basically **unwrapping the spokes into vertical lines**, making it a lot easier to handle.
This step involves some bilinear interpolation which should not change the resolution result, but helps the algorithm to avoid sampling artefacts.
The resulting image is than reshaped into a stack of linepairs and averaged in this direction. By this operation, the **average** resolution on a concentric circle is computed.
The Michelson contrast is computed by (I_max-I_min)/(I_max+I_min) for every radial position.
The result is normalized by the maximum contrast in the field (close to DC, i.e. very low spatial frequencies).
The radial axis is converted to a millimeter axis based on pixel size, radial supersampling and optical magnification (user specified).

# How to use it
Run ```example_minimal.m```, should be straightforward.
For a more in-depth example that extracts PSFs and verifies the implementation, check out ```example_extended.m```

# Tipps for good siemens star recordings
- Illuminate the siemens star from behind with a homogeneous light source to provide even illumination
- For nearly monochromatic applications (e.g. narrowband LIF detection schemes), use a lightsource that has spectral lines (for OH, e.g. HgAr Penray has some lines around 300-320nm) and place your filter in front of your lens. This prevents chromatic abberation, thus preventing your resolution to look worse than it actually is.
- Acquire a dark frame (everything in place, illumination off)
- Acquire a flat field, ideally with the same light-source
- Correct image with dark frame and flat field
- Crop image to a region, where only the spokes of the star are visible. Use as much of the star as possible, because the outer parts are used as contrast==1 reference by the script!

# Use cases for this script
- Get a reliable figure of your resolution for a publication (some use 10% contrast as resolution limit, some 50%, some 17%....)
- If your camera software supports hooking in MATLAB scripts (like e.g. LabVIEW does), this allows focusing on a live MTF