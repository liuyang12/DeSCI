-------------------------------------------------------------------

                V-BM4D software for video denoising
            Public release ver. 1.0  (11 December 2014)

-------------------------------------------------------------------

Copyright (c) 2010-2014 Tampere University of Technology. 
All rights reserved.
This work should be used for nonprofit purposes only.

Authors:                     Matteo Maggioni
                             Alessandro Foi


V-BM4D web page:             http://www.cs.tut.fi/~foi/GCF-BM3D


-------------------------------------------------------------------
 Contents
-------------------------------------------------------------------

The package contains these files

*) demo_denoising.m         : denoising demo script
*) vbm4d.m                  : V-BM4D video denoising filter
*) read_video.m             : reads data from video file

-------------------------------------------------------------------
 Installation & Usage
-------------------------------------------------------------------

Unzip VBM4D.zip (contains codes) in a folder that is in the MATLAB 
path. Execute the script "demo_denoising.m" to run a video denoising
demo. You can freely modify the parameters involved in the filtering
at the beginning of the demo.

-------------------------------------------------------------------
 Requirements
-------------------------------------------------------------------

*) MS Windows 64 bit, Linux 64 bit or Mac OS X 64 bit
*) Matlab R2011b or later with installed:
   -- Image Processing Toolbox (only for visualization with "implay")
   -- Signal Processing Toolbox (only for non-default transforms in VBM4D)
   -- Wavelet Toolbox (only for non-default transforms in VBM4D)

-------------------------------------------------------------------
 Change log
-------------------------------------------------------------------
v1.0   (11 December 2014)
 + initial version

-------------------------------------------------------------------
 References
-------------------------------------------------------------------

[1] M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising 
    Using Separable 4D Nonlocal Spatiotemporal Transforms",  Proc.
	SPIE Electronic Imaging 2011, Image Processing: Algorithms and
	Systems IX, 7870-2, San Francisco (CA), USA, January 2011.
	doi:10.1117/12.872569

[2] M. Maggioni, G. Boracchi, A. Foi, K. Egiazarian, "Video Denoising, 
    Deblocking and Enhancement Through Separable 4-D Nonlocal 
    Spatiotemporal Transforms", IEEE Trans. Image Proc., 
    vol. 21, no. 9, Sep. 2012.   doi:10.1109/TIP.2012.2199324
 
-------------------------------------------------------------------
 Disclaimer
-------------------------------------------------------------------

Any unauthorized use of these routines for industrial or profit-
oriented activities is expressively prohibited. By downloading 
and/or using any of these files, you implicitly agree to all the 
terms of the TUT limited license, as specified in the document
Legal_Notice.txt (included in this package) and online at
http://www.cs.tut.fi/~foi/GCF-BM3D/legal_notice.html

-------------------------------------------------------------------
 Feedback
-------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact    Matteo Maggioni   at  matteo.maggioni<at>tut.fi


