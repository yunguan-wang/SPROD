"""
This script will extract two intensity features, median and std, as well as six haralick features from image.
It is tested to work with Visium he and if images.
A conda env is available in the project folder. 

# input folder must have the raw tif image and the spot_metadata csv file.
# The spot metadata file MUST have 'X', 'Y', and 'Spor_radius' columns

# example usage:
conda activate /project/shared/xiao_wang/projects/MOCCA/conda_env
python /project/shared/xiao_wang/projects/MOCCA/code/feature_extraction.py \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune \
    if \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune/intermediate \
    
"""
import numpy as np
import pandas as pd
import os
import sys
import argparse
from skimage import io, img_as_float32, morphology, exposure
from skimage.feature import greycomatrix, greycoprops
from itertools import product
from sklearn.preprocessing import minmax_scale
from skimage.color import separate_stains, hdx_from_rgb


def extract_img_features(
    input_path,
    input_type,
    output_path,
    img=None,
    img_meta=None,
    feature_mask_shape="spot",
):
    """
    Extract features from image. Works with IF or HE image from Visium tif files.
    For block feature, a square will be drawn around each spot. Since it is bigger than 
    the spot itself, it is more suitable to extract texture features. 
    For Spot feature, only area in the actual sequencing spot will be uses. 
    It is more suitable to extract intensity features.

    Parameters
    ----------
    input_path : str
        input folder containing all necessary files. 
    input_type : str
        input image type, select from {'if','he'}.
    output_path : str
        output folder path.
    img : None or np.array, optional
        alternative input for image, will override input_path.
    img_meta : None or np.array, optional
        alternative input for image metadata, will override input_path.
    feature_mask_shape : {'spot', 'block'}
        type of feature extracted.
    """
    intensity_fn = os.path.join(
        os.path.abspath(output_path),
        "{}_level_texture_features.csv".format(feature_mask_shape)
    )
    texture_fn = os.path.join(
        os.path.abspath(output_path),
        "{}_level_intensity_features.csv".format(feature_mask_shape)
        )
    if (os.path.exists(intensity_fn)) == (os.path.exists(texture_fn)) == True:
        print('Features are already extracted.')
        return
    if img_meta is None:
        img_meta = pd.read_csv(
            os.path.join(input_path,"Spot_metadata.csv"), index_col=0)
    if img is None:
        img_tif = [x for x in os.listdir(input_path) if "tif" in x][0]
        img_tif = os.path.join(input_path, img_tif)
        if input_type == "if":
            # the indexing is a workaround for the strange Visium if image channels.
            img = io.imread(img_tif)
            img = img_as_float32(img)
            img = (255 * img).astype("uint8")
        else:
            img = io.imread(img_tif)
            # normalize image with color deconv
            print('Normalizing image...')
            img = separate_stains(img, hdx_from_rgb)
            img = minmax_scale(img.reshape(-1, 3)).reshape(img.shape)
            img = np.clip(img, 0, 1)
            img = exposure.equalize_adapthist(img, clip_limit=0.01)
            img = (255 * img).astype("uint8")
    # Hard coded type of Haralick features and Angles for searching for neighboring pixels
    # hard coded number of angles to be 4, meaning horizontal, vertical and two diagonal directions.
    # extracting block shaped features
    if feature_mask_shape == "block":
        tmp = img_meta.sort_values(["Row", "Col"])
        block_y = int(np.median(tmp.Y.values[2:-1] - tmp.Y.values[1:-2]) // 2)
        tmp = img_meta.sort_values(["Col", "Row"])
        block_x = int(np.median(tmp.X.values[2:-1] - tmp.X.values[1:-2]) // 2)
        block_r = min(block_x, block_y)
        block_x = block_y = block_r
    print("Prossessing {}".format(input_path))
    feature_set = [
        "contrast",
        "dissimilarity",
        "homogeneity",
        "ASM",
        "energy",
        "correlation",
    ]
    text_features = []
    intensity_features = []
    for i in range(img_meta.shape[0]):
        if (i + 1) % 100 == 0:
            print("Processing {} spot out of {} spots".format(i + 1, img_meta.shape[0]))
        row = img_meta.iloc[i]
        x, y, r = row[["X", "Y", "Spot_radius"]].astype(int)
        if feature_mask_shape == "spot":
            spot_img = img[x - r : x + r + 1, y - r : y + r + 1]
            spot_mask = morphology.disk(r)
            # only use the spot, not the bbox
            spot_img = np.einsum("ij,ijk->ijk", spot_mask, spot_img)
        else:
            spot_img = img[x - block_x : x + block_x + 1, y - block_y : y + block_y + 1]
            spot_mask = np.ones_like(spot_img[:, :, 0], dtype="bool")

        # extract texture features
        ith_texture_f = []
        for c in range(img.shape[2]):
            glcm = greycomatrix(
                spot_img[:, :, c],
                distances=[1],
                # Angles are arranged in a counter clockwise manner, in radian.
                angles=[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
                levels=256,
                symmetric=True,
                normed=False,
            )
            glcm = glcm[1:, 1:]
            glcm = glcm / np.sum(glcm, axis=(0, 1))
            for feature_name in feature_set:
                ith_texture_f += greycoprops(glcm, feature_name)[0].tolist()
        # The first 6 features are intensity features, and the rest are Haralicks.
        text_features.append(ith_texture_f)

        # extract intensity features
        int_low = 0.2
        int_high = 0.8
        int_step = 0.1
        q_bins = np.arange(int_low, int_high, int_step)
        ith_int_f = []
        for c in range(img.shape[2]):
            for t in q_bins:
                ith_int_f.append(np.quantile(spot_img[:, :, c][spot_mask == True], t))
        intensity_features.append(ith_int_f)

    # Naming the features. f stands for channels, A stands for angles.
    # construct texture feature table
    channels = ["f" + str(i) for i in range(img.shape[2])]
    col_names = product(channels, feature_set, ["A1", "A2", "A3", "A4"])
    col_names = ["_".join(x) for x in col_names]
    text_features = pd.DataFrame(text_features, index=img_meta.index, columns=col_names)
    # construct intensity feature table
    intensity_features = pd.DataFrame(
        intensity_features,
        index=img_meta.index,
        columns=[
            "_".join(x) for x in product(channels, ["{:.1f}".format(x) for x in q_bins])
        ],
    )
    text_features.to_csv(intensity_fn)
    intensity_features.to_csv(texture_fn)
    return text_features, intensity_features

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extract features from image. Works with IF or HE image from Visium tif files.\
            For block feature, a square will be drawn around each spot. Since it is bigger than \
            the spot itself, it is more suitable to extract texture features. \
            For Spot feature, only area in the actual sequencing spot will be uses. \
                It is more suitable to extract intensity features.",
    )

    parser.add_argument(
        "input_path", type=str, help="Input folder containing all necessary files."
    )
    parser.add_argument(
        "input_type", type=str, help="Input image type, select from {'if','he'}."
    )
    parser.add_argument(
        "--output_path",
        "-o",
        type=str,
        default="./intermediate",
        help="Output folder path.",
    )
    parser.add_argument(
        "--feature_mask_shape",
        "-m",
        type=str,
        default="spot",
        help="Type of feature extracted. {'spot', 'block'}",
    )

    args = parser.parse_args()
    # Todo: decide if checking metadata function should be added, or force the user to provide the correct format.
    input_path = os.path.abspath(args.input_path)
    input_type = args.input_type
    output_path = os.path.abspath(args.output_path)
    feature_mask_shape = args.feature_mask_shape

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    print("Processing {}".format(input_path))
    output_fn = output_path + "/{}_level_texture_features.csv".format(
        feature_mask_shape
    )
    if os.path.exists(output_fn):
        print("Image feature already extracted. Stop operation.")
    else:
        _ = extract_img_features(
            input_path, input_type, output_path,
            feature_mask_shape=feature_mask_shape)

# Debugging codes
# f_t, f_i = extract_img_features(
#     '','HE','',img = center_patch_sep, img_meta = center_patch_meta,
#     feature_mask_shape='spot')

# # valid_cols = [x for x in f_block.columns if 'homo' not in x]
# # valid_cols = [x for x in valid_cols if 'corr' not in x]
# # f_pca_b = PCA(n_components=10).fit_transform(scale(f_block[valid_cols]))
# # clusters_b = [str(x) for x in k_means(f_pca_b,5)[1]]

# f_pca_b = PCA(n_components=10).fit_transform(scale(f_i))
# clusters_b = [str(x) for x in k_means(f_pca_b,3)[1]]

# _ = sns.mpl.pyplot.figure(figsize = (16,16))
# ax = sns.scatterplot(
#     y = center_patch_meta.X, x = center_patch_meta.Y,
#     hue=clusters_b, hue_order = sorted(set(clusters_b)), markers = ['h'],
#     linewidth=0, alpha=0.25, s=250)
# ax.set_facecolor('grey')
# io.imshow(center_patch[:,:,0], alpha=0.9)
# io.imshow(center_patch_sep[:,:,0],cmap='gray')