# Decompress Snapshot Compressive Imaging (DeSCI)
This repository contains the code for the paper **Rank Minimization for Snapshot Compressive Imaging** (preprint 2018) by [Yang Liu*](https://liuyang12.github.io/), [Xin Yuan*](https://www.bell-labs.com/usr/x.yuan), [Jinli Suo](https://sites.google.com/site/suojinli/), [David J. Brady](https://ece.duke.edu/faculty/david-brady), and [Qionghai Dai](http://media.au.tsinghua.edu.cn/qhdai_new.html) (*Equal contributions).
[[pdf]](https://arxiv.org/pdf/1807.07837.pdf "arXiv preprint")   [[project page]](https://github.com/liuyang12/DeSCI "github repository")

![Video comparison of DeSCI.](/results/video/desci_gmm_gaptv_kobe32.gif?raw=true)

Figure 1. Reconstructed `Kobe` video of DeSCI compared with the state-of-the-art methods, *i.e.*, GMM-TP (TIP'14), MMLE-GMM (TIP'15), MMLE-MFA (TIP'15), and GAP-TV (ICIP'16).

Snapshot compressive imaging (SCI) refers to encoding the three- or higher- dimensional data in a snapshot with distinct masks (or coded aperutre) for each slice of the data. Decompress snapshot compressive imaging (DeSCI) then exploits the nonlocal self-similarity of natural scenes and apply an alternating minimization algorithm to solve the ill-posed problem. 

Performance boost compared with the state-of-the-art methods includes more than 4 dB in terms of PSNR for simulated data and significant improvement for real data, which is addressed in the DeSCI paper. A video comparison of the proposed DeSCI method with the state-of-the-art algorithms is shown in Figure 1. 

This code only contains the simulated high-speed video `Kobe` dataset encoded with shifting masks following the coded aperture shapshot compressive imager (CACTI). It can also be used for snapshot-hyperspectral compressive imaging, such as the coded aperture snapshot spectral imaging (CASSI) system, as shown in the paper. Code for the CASSI data and the CACTI and CASSI data from real systems are available upon request. 

## Usage
0. Requirements are MATLAB(R) with Parallel Computing Toolbox (`parfor` for multi-CPU acceleration).
1. Download this repository via git
```
git clone https://github.com/liuyang12/DeSCI
```
or download the [zip file](https://github.com/liuyang12/DeSCI/archive/master.zip) manually.

2. Test the DeSCI algorithm via
```matlab
test_desci.m
```
and (optionally) video demonstrate the reconstruction results (after running `test_desci.m`) via
```matlab
demo_desci_video.m
```

## Structure of directories

| directory  | discription  |
| :--------: | :----------- | 
| `algorithms` | MATLAB functions of main algorithms proposed in the paper (original) | 
| `packages`   | algorithms adapted from the state-of-art algorithms (adapted)|
| `dataset`    | data used for reconstruction (simulated) |
| `results*`   | results of reconstruction (after reconstruction) |
| `utils`      | utility functions |

## Platform
The test platform is MATLAB(R) 2017b operating on Ubuntu 16.04 LTS (x64) with an Intel(R) Core(TM) 18-core processor at 2.60 GHz and 128 GB RAM. It can run on any machine with MATLAB(R) and Parallel Computing Toolbox, operating on Windows(R), Linux, or Mac OS. No GPU is needed to run this code.

It could take hours to run a single measurement depending on the number of frames collapsed into a single measurement and the number of CPU cores of the machine. This is due to the time-consuming iterations of block matching and low-rank approximation via singular value decomposition for weighted nuclear norm minimization of recovered video/spectral patches. The computation issue might be addressed by using generative adversarial networks (GAN) for block matching, and truncated SVD for low-rank estimation. *Notice: Please donot wait for the results immediately after getting the DeSCI code to run.*

## Citation
```
@article{Liu18DeSCI,
   author = {Liu, Yang and Yuan, Xin and Suo, Jinli and Brady, David J. and Dai, Qionghai},
   title = {Rank Minimization for Snapshot Compressive Imaging},
   journal = {arXiv (preprint)},
   pages = {1807.07837},
   year = {2018},
   type = {Journal Article}
}
```

## Contact
[Yang Liu, Tsinghua University](mailto:y-liu16@mails.tsinghua.edu.cn "Yang Liu, Tsinghua University") 

[Xin Yuan, Bell Labs](mailto:xyuan@bell-labs.com "Xin Yuan, Bell labs")  
