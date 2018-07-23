# denoiser
This directory contains the denoisers used in this paper. Total variation (TV) denoiser for GAP-TV, Weighted nuclear norm minimization (WNNM) denoiser for DeSCI, and video block-matching and 4-D filtering (V-BM4D) denoiser for DeSCI-VBM4D.
## TV
The TV denoiser is adapted from the GAP-TV code by [Xin Yuan](https://www.bell-labs.com/usr/x.yuan "Xin Yuan, Bell labs"). Please cite the GAP-TV paper if you use the TV denoiser for GAP-TV.
```
[1] X. Yuan, "Generalized alternating projection based total variation minimization for compressive sensing," in IEEE International Conference on Image Processing (ICIP), 2016, pp. 2539-2543.
```

## WNNM
The WNNM denoiser for video denoising is adapted from the WNNM denoiser for image denoising proposed by [Shuhang Gu](https://sites.google.com/site/shuhanggu/home), [Lei Zhang](http://www4.comp.polyu.edu.hk/~cslzhang/), Wangmeng Zuo, and Xiangchu Feng. We download the original WNNM image denoiser from the author(s) personal website (available on http://www4.comp.polyu.edu.hk/~cslzhang/code/WNNM_code.zip).
Please cite the WNNM paper if you use the WNNM image denoiser.
```
[2] S. Gu, L. Zhang, W. Zuo, and X. Feng, "Weighted Nuclear Norm Minimization with Application to Image Denoising," in IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014, pp. 2862-2869.
[3] S. Gu, Q. Xie, D. Meng, W. Zuo, X. Feng, and L. Zhang, "Weighted Nuclear Norm Minimization and Its Applications to Low Level Vision," International Journal of Computer Vision, vol. 121, no. 2, pp. 183-208, 2017.
```
Please refer to the DeSCI paper if you use the WNNM video denoiser. 

## VBM4D
The VBM4D denoiser is downloaded from the [BM3D](http://www.cs.tut.fi/~foi/GCF-BM3D/ "block-matching and 3-D filtering") website (available on http://www.cs.tut.fi/~foi/GCF-BM3D/VBM4D_v1.zip). We use VBM4D instead of VBM3D nor BM4D because VBM4D contains noise estimation, which eliminates the need to estimate the noise level of the alternating projection process compared with VBM3D, and VBM4D outperforms BM4D in terms of videos.
Please cite the VBM4D paper is you use the VBM4D denoiser.
```
[4] M. Maggioni, G. Boracchi, A. Foi, and K. Egiazarian, "Video Denoising, Deblocking, and Enhancement Through Separable 4-D Nonlocal Spatiotemporal Transforms," IEEE Transactions on Image Processing, vol. 21, no. 9, pp. 3952-3966, May 15 2012.
```



