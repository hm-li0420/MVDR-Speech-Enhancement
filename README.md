# MVDR-Speech-Enhancement
这是一个C语言代码实现的MVDR语音增强代码。
## 目录

1. [背景](#背景)
2. [安装](#安装)
3. [使用方法](#使用方法)

## 背景

这是一个C语言实现的MVDR语音增强，整体项目分为两个部分：
1. 基于Capon搜索的声源定位(DOA)实现(参考main.c里的get_angle函数)
2. 将估计的目标声源位置结合MVDR实现语音增强(参考main.c里的beamform函数)

## 使用方法
### 文件路径
include: 头文件路径  
src：c文件路径  
pcm_files：噪声文件例子：45为声源在45度，90为声源在90度  
matlab：MVDR和Capon算法的matlab实现，以及合成点声源的matlab代码 
```bash
│  CMakeLists.txt
│  README.md
│
├─.vscode
│      settings.json
│      tasks.json
│
├─build
│  │  CMakeCache.txt
│  │  cmake_install.cmake
│  │  compile_commands.json
│  │  Makefile
│  │  mic_array.exe
│  ├─.cmake
│  └─CMakeFiles
│
├─include
│      bf.h
│      kiss_fft.h
│      kiss_fftr.h
│      kiss_fft_log.h
│      matrix.h
│      _kiss_fft_guts.h
│
├─matlab
│  ├─main
│  ├─Simulation
│  └─STFT
├─out
│
├─pcm_files
│  ├─45
│  │
│  └─90
│
└─src
        bf.c
        kiss_fftr.c
        main.c
        matrix.c
```
### 运行方法
在build目录找到可执行文件mic_array.exe，运行以下指令：
```bash
$ mic_array mic_0.pcm mic_1.pcm mic2.pcm mic3.pcm out.pcm
```  
mic_0.pcm mic_1.pcm mic2.pcm mic3.pcm分别为四个麦克风的录制音频文件，格式为PCM，out.pcm为输出音频。

## 安装
项目环境为CMakeLists，需要使用者自己搭建好cmake环境
```bash
$ git clone git@github.com:hm-li0420/MVDR-Speech-Enhancement.git
```