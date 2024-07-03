# Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery

Implementation for MICCAI 2024 paper **[Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery](https://)** by [Jiahe Chen](http://), Etsuko Kobayashi, Ichiro Sakuma, and Naoki Tomii.

A novel problem setting and model-based framework for online tissue deformation recovery.

**[\[Paper\]](https://) [\[Sample Dataset\]](https://forms)**

## Demo
https://github.com/bmpelab/trans_window_panoramic_impasto/assets/44491917/c95251c0-c616-4243-a203-20599780d6df


## Problem formulation
![overview](https://github.com/bmpelab/trans_window_panoramic_impasto/assets/44491917/5662d76e-e410-44e6-a462-d027681dfae1)


## Code preparation



<!--We recommend using Miniconda to set up an environment:-->

```bash
cd xxx
```

<!--We managed to test our code on Ubuntu 18.04 with Python 3.6 and CUDA 10.2.-->

## Dataset preparation

You can download the well-prepared dataset from XXX (derived from our ex vivo dataset) and XXX (derived from Hamlyn center dataset). For the EndoNeRF dataset, please download from [here](https://github.com/med-air/EndoNeRF) and process following the instructions below. For more datasets from the Hamlyn Centre, please visit [Hamlyn Centre Laparoscopic / Endoscopic Video Datasets](https://hamlyn.doc.ic.ac.uk/vision/)

The structure of the dataset should looks like:



image and camrea parameters download, ours, endonerf, hamlyn

Note that, in case you capture images by your own system, image rectification may be necessary using the calibrated camera parameters. For image rectification, please follows [OpenCV_stereoRectify](https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html#ga617b1685d4059c6040827800e72ad2b6).

3D reconstruction, 2D optical flow tracking

instrument mask generation

scene flow generation, point map generation

<!--To test our method on your own data, prepare a data directory organized in the following structure:-->

```
+ data1
    |
    |+ depth/           # depth maps
    |+ masks/           # binary tool masks
    |+ images/          # rgb images
    |+ pose_bounds.npy  # camera poses & intrinsics in LLFF format
```

<!--In our experiments, stereo depth maps are obtained by [STTR-Light](https://github.com/mli0603/stereo-transformer/tree/sttr-light) and tool masks are extracted manually. Alternatively, you can use segmentation networks, e.g., [MF-TAPNet](https://github.com/YuemingJin/MF-TAPNet), to extract tool masks. The `pose_bounds.npy` file saves camera poses and intrinsics in [LLFF format](https://github.com/Fyusion/LLFF#using-your-own-poses-without-running-colmap). In our single-viewpoint setting, we set all camera poses to identity matrices to avoid interference of ill-calibrated poses.-->

## Run

Although the data of all frames are prepared, the code iteratively load the data frame by frame; thus, it is a online approach that making use of only the current and previous frames information.

```bash
xxx
```

The reconstructed point clouds will be saved to `logs/{expname}/reconstructed_pcds_{epoch}`. For more options of this reconstruction script, type `python endo_pc_reconstruction.py -h`.

We also build a visualizer to play point cloud animations. To display reconstructed point clouds, type the command as follows.

```bash
python vis_pc.py --pc_dir logs/{expname}/reconstructed_pcds_{epoch}
```

Type `python vis_pc.py -h` for more options of the visualizer.

## Evaluation

We provide ex vivo dataset with ground truth of the deformed surface obtained by 3D scanning. Please download the dataset from [here](http://). Since the structure

```bash
python run_endonerf.py --config configs/{your_config_file}.txt --render_only
```

The rendered images will be saved to `logs/{expname}/renderonly_path_fixidentity_{epoch}/estim/`. Then, you can type the command below to acquire quantitative results:

```bash
python eval_rgb.py --gt_dir /path/to/data/images --mask_dir /path/to/data/gt_masks --img_dir logs/{expname}/renderonly_path_fixidentity_{epoch}/estim/
```

Note that we only evaluate photometric errors due to the difficulties in collecting geometric ground truth. 

## Bibtex

If you find this work helpful, you can cite our paper as follows:

```
@InProceedings{chen_trans-window_2024,
 author = {Chen, Jiahe and Kobayashi, Etsuko and Sakuma, Ichiro and Tomii, Naoki},
 title = {{Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery}},
 booktitle = {Medical {{Image Computing}} and {{Computer Assisted Intervention}} -- {{MICCAI}} 2024},
 pages = {xxx--xxx},
 year = {2024},
}
```

If you use our ex vivo datasets, please cite this work also:

```
xxx
```

If you use the Hamlyn center datasets, please appropriately cite their works following the instructions in [Hamlyn Centre Laparoscopic / Endoscopic Video Datasets](https://hamlyn.doc.ic.ac.uk/vision/).

## Acknowledgement

- readflow code
- natsort
