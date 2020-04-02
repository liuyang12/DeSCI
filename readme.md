# Decompress Snapshot Compressive Imaging (DeSCI)
This repository contains the MATLAB code for the paper **Rank Minimization for Snapshot Compressive Imaging** (***IEEE Transactions on Pattern Analysis and Machine Intelligence*** 2019) by [Yang Liu*](https://liuyang12.github.io/), [Xin Yuan*](https://www.bell-labs.com/usr/x.yuan), [Jinli Suo](https://sites.google.com/site/suojinli/), [David J. Brady](https://ece.duke.edu/faculty/david-brady), and [Qionghai Dai](http://media.au.tsinghua.edu.cn/) (*Equal contributions).
[[pdf]](https://arxiv.org/pdf/1807.07837.pdf "arXiv preprint")   [[github]](https://github.com/liuyang12/DeSCI "github repository")  [[data]](https://drive.google.com/open?id=1d2uh9nuOL5Z7WnEQJ5HZSDMWK2VAT9sH "real data and code")  [[arXiv]](https://arxiv.org/abs/1807.07837 "arXiv preprint")   [[doi]](https://doi.org/10.1109/TPAMI.2018.2873587 "10.1109/TPAMI.2018.2873587")  

**[New] Real data and associated code are available at [this link on Google Drive](https://drive.google.com/open?id=1d2uh9nuOL5Z7WnEQJ5HZSDMWK2VAT9sH).** Note that the code for real data is not tested but with raw results as in the paper. Please refer to the readme file for the original source(s) of the real data.

<p align="center">
<img src="https://github.com/liuyang12/DeSCI/blob/master/results/video/desci_gmm_gaptv_kobe32.gif?raw=true">
</p>

Figure 1. Reconstructed `Kobe` video using DeSCI compared with the state-of-the-art methods, *i.e.*, GMM-TP (TIP'14), MMLE-GMM (TIP'15), MMLE-MFA (TIP'15), and GAP-TV (ICIP'16). Here, 8 video frames are encoded in a single measurement and a total of 32 frames are presented by reconstructing 4 snapshot measurements. The `Kobe` video is used in the MMLE-GMM [paper](https://doi.org/10.1109/TIP.2014.2365720 "Compressive Sensing by Learning a Gaussian Mixture Model From Measurements, TIP'15").

<p align="center">
<img src="https://github.com/liuyang12/DeSCI/raw/master/results/spectral/desci_cassi_toy_spectra.png?raw=true" width="500">
</p>

Figure 2. Reconstructed spectra of `toy` hyperspectral images using DeSCI compared with GAP-TV (ICIP'16). Here, 32 spectral frames are encoded in a single measurement. The `toy` hyperspectral images are from the CAVE multispectral image [database](http://www.cs.columbia.edu/CAVE/databases/multispectral/ "Multispectral Image Database | CAVE | Columbia University"). 

## Snapshot compressive imaging (SCI)
Snapshot compressive imaging (SCI) refers to encoding the three- or higher- dimensional data in a snapshot with distinct mask (or coded aperture) for each slice of the data. Decompress snapshot compressive imaging (DeSCI) then exploits the nonlocal self-similarity of natural scenes and applys an alternating minimization algorithm to solve the ill-posed problem. 

Performance boost compared with the state-of-the-art methods includes more than 4 dB in terms of PSNR for simulated data and significant improvement for real data, which is addressed in the DeSCI paper. A video comparison of the proposed DeSCI method with the state-of-the-art algorithms is shown in Figure 1. 

This code contains the simulated high-speed video `Kobe` dataset encoded with shifting masks following the coded aperture compressive temporal imager (CACTI), and the simulated hyperspectral `toy` dataset following the coded aperture shapshot spectral imager (CASSI). A brief review of the SCI systems, including CASSI and CACTI, and compressive light-field imaging is shown in Section 2 of the DeSCI paper. DeSCI could also be adapted for other compressive imaging systems with minor modifications, since we only need to change the sensing matrix for various coding strategies and the nonlocal self-similarity always hold for natural scenes. Code for the CACTI and CASSI data from real systems are available upon request. 

## Usage
### Download the DeSCI repository
0. Requirements are MATLAB(R) with Parallel Computing Toolbox (`parfor` for multi-CPU acceleration).
1. Download this repository via git
```
git clone https://github.com/liuyang12/DeSCI
```
or download the [zip file](https://github.com/liuyang12/DeSCI/archive/master.zip) manually.

### Run DeSCI on high-speed video
#### `Kobe` video data
2. Test the DeSCI algorithm (for high-speed imaging, that is CACTI on `Kobe` dataset) via
```matlab
test_desci.m
```
and (optionally) video demonstrate the reconstruction results (after running `test_desci.m`) via
```matlab
./figures/fig_desci_video.m
```

### Run DeSCI on hyperspectral images
#### `toy` hyperspectral data
3. Test the DeSCI algorithm (for hyperspectral imaging, that is CASSI on `toy` dataset) via
```matlab
test_desci_cassi.m
```
and (optionally) video demonstrate the reconstruction results (after running `test_desci_cassi.m`) via
```matlab
./figures/fig_desci_cassi.m
```

#### `bird` hyperspectral data
[Optional] Test the DeSCI algorithm (for hyperspectral imaging, that is CASSI on `bird` dataset) via
```matlab
./figures/test_desci_cassi_bird.m
```
and (optionally) video demonstrate the reconstruction results (after running `./figures/test_desci_cassi_bird.m`) via
```matlab
./figures/fig_desci_cassi_bird.m
```

Note that running the `bird` dataset is time- and memory- consuming for the size of the images is $1021\times 703$ with 24 spectral bands. Please ensure that the machine has *at least* 8 physical CPU cores and 32 GB memory to run `./figures/test_desci_cassi_bird.m`, where more CPU cores and memory help for faster reconstruction.

## Structure of directories

| directory  | description  |
| :--------: | :----------- | 
| `algorithms` | MATLAB functions of main algorithms proposed in the paper (original) | 
| `figures`    | MATLAB scripts to reproduce the results and figures in the paper (original) |
| `packages`   | algorithms adapted from the state-of-art algorithms (adapted)|
| `dataset`    | data used for reconstruction (simulated) |
| `results`    | results of reconstruction (after reconstruction) |
| `utils`      | utility functions |

## Platform
The test platform is MATLAB(R) 2017b (and 2018b) operating on Ubuntu 16.04 LTS (x64) with an Intel(R) Core(TM) 18-core processor at 2.60 GHz and 128 GB RAM. It can run on any machine with MATLAB(R) and Parallel Computing Toolbox, operating on Windows(R), Linux, or Mac OS. No GPU is needed to run this code.

It could take hours to run a single measurement depending on the number of frames collapsed into a single measurement and the number of CPU cores of the machine. This is due to the time-consuming iterations of block matching and low-rank approximation via singular value decomposition for weighted nuclear norm minimization of recovered video/spectral patches. The computation issue might be addressed by using generative adversarial networks (GAN) for block matching, and truncated SVD for low-rank estimation. *Notice: Please donot wait for the results immediately after getting the DeSCI code to run.*

## Citation
```
@article{Liu19DeSCI,
   author  = {Liu, Yang and Yuan, Xin and Suo, Jinli and Brady, David J. and Dai, Qionghai},
   title   = {Rank Minimization for Snapshot Compressive Imaging},
   journal = {IEEE Trans. Pattern Anal. Mach. Intell.},
   doi     = {10.1109/TPAMI.2018.2873587},
   year    = {2019},
   volume  = {41},
   number  = {12},
   pages   = {2990 - 3006},
   url     = {https://doi.org/10.1109/TPAMI.2018.2873587},
   type    = {Journal Article}
}
```

## Contact
[Yang Liu, Tsinghua University](mailto:yliu@csail.mit.edu "Yang Liu, Tsinghua University") 

[Xin Yuan, Bell Labs](mailto:xyuan@bell-labs.com "Xin Yuan, Bell labs")  
