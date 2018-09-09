# a2

## If part 4 is not running (merge conflicts) look at branch: minhtong to run
## Part 1: Custom Billboards (Image Warping and Homogeneous Matrix)

#### Running the Code
```
make
./a2 part1 image_for_billboard.png
```

#### Part 1.1: Lincoln Transformation
We had created a forward transformation function (apply_transformation) to test the given transformation matrix on the Lincoln memorial image. Afterwards we created a function to apply inverse warping (apply_transformation_inverse_warping) to create a better image that would match the example output shown in the assignment. This function created a result image of the same size as the input filled with 0s. The transformation matrix passed in is inverted. Then each location in the black result image is transformed by the inverted transformation matrix to find the original pixel location and the RGB values. 

##### 1.1 Results
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/lincoln_warped.png "lincoln")

The result could have been improved with the use of bilinear interpolation.

#### Part 1.2: Computing Transformation Matrix
The get_transformation function was created to compute the transformation function when given 4 corresponding coordinate pairs. This function assumes that the 9th position in the matrix is 1. The function solves a system of equations for 8 unknowns (9th is assumed as 1) and the results are returned as a 3x3 transformation matrix.
The computed transformation matrix is also printed to the console.
Resource used:  http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/EPSRC_SSAZ/node11.html

##### 1.2 Results
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2book_result.png "book")

#### Part 1.3: Billboards
For this section the overlay_images function was created. It takes the two images and their corresponding points. Then it uses the get_transformation function to find the transformation matrix required to transform the image to the billboard position. Next it loops through the billboard image and when it comes to the location of where the transformation would appear (a position on the billboard), it finds the corresponding pixel value in the original image to be placed on the billboard.

##### 1.3 Results
images/part1/nlss.jpg reference:  https://www.reddit.com/r/NLSSCircleJerk/comments/7uu4fh/i_replaced_harambe_with_hafu_on_the_overlay/

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2synthetic_billboard1.png "b1")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2synthetic_billboard2.png "b2")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2synthetic_billboard3.png "b3")

#### Issues during development of Part 1
1) It took some time to get an understanding of the CImg library.
2) Forgot to divide out by w (3rd dimension value) when transforming.

## Part 2: Blending

#### Running the Code
In order to run this part's code enter the following commands in the project directory
```
make
./a2 part2 image1.png image2.png mask.png
```
with image1, image2, and mask all being the same size images.

#### Overview
In order to make all of our output images clear and human viewable, they have been normalized between 0 and 255.

The program is implemented in the same form as the updated instructions on Piazza.

A few helper functions have also been implemented to make the code more readable, including

- upscale, for the gaussian and laplacian pyramids

- half_mask, for automatically generating masks split down the middle for two images

- half_blend, for atomatically blending two images down the middle

- half_blend_panel, for creating a display of two images blended down the middle


- gaussian_pyramid, returning a list containing each level of the gaussian pyramid

- laplacian_pyramid, returning a list containing each level of the laplacian pyramid

- blended_pyramid, returning a list containing the non-combined blended pyramid


The basic function for blending to images is

- blend_images, which takes both images, a mask, a filter, and the number of pyramid levels

which returns the final blended image.

The blend_images function only requires the mask to be the same size as the smallest image in order to function properly, if the two images are different sizes, the larger one is resized to the smaller one's size and centered around the middle of the image.

The blending working on the apple and orange example can be seen below, along with the intermediate pyramids.

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/apple_dis.jpg "apple")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/blended_dis.png "blended")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/orange_dis.jpg "orange")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/mask_dis.jpg "mask")

Apple pyramids

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau1_000000.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau1_000001.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau1_000002.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau1_000003.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau1_000004.png "gaussian pyramid")


![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap1_000000.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap1_000001.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap1_000002.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap1_000003.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap1_000004.png "laplacian pyramid")


Orange pyramids

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau2_000000.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/part2/gau2_000001.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau2_000002.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau2_000003.png "gaussian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/gau2_000004.png "gaussian pyramid")

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap2_000000.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap2_000001.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap2_000002.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap2_000003.png "laplacian pyramid")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/lap2_000004.png "laplacian pyramid")

Another interesting application of this blending program is the prototyping of future wardrobe options, as can be seen below.

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/tru_dis.png "trump")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/blendedtrump_good_dis.png "blended")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/mel_dis.png "melania")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/trumpmask_good_dis.png "mask")

The mask made a huge difference in the quality of the result, as can be seen with the same images blended with a more simple mask.

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/tru_dis.png "trump")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/blendedtrump_half_dis.png "blended")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/mel_dis.png "melania")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/images/part2/trumpmask_half_dis.png "mask")

## Part 3

#### Running the Code
Run make to compile the program and then input the command
```
./a2 part3 eiffel_18.jpg eiffel_19.jpg
```
with eiffel_18.jpg being source image and eiffel_19.png being destination image between which the matching points needs to be computed. It outputs two images one with SIFT matching pairs "part3_sift_matches.jpg" and another with matching pairs after applying RANSAC "part3_ransac_matches.jpg"

#### Overview
In this task we are trying to find the matching feature points between two given images so as to identify if both the images belong to same location/place/area. First we apply SIFT on both the images and retrieve feature points (128-dimensional descriptor) in each image. For each of the descriptor in first image we find the ratio of the Euclidean distances between closest descriptor and next closest descriptor from the second image. The range of these ratio values will be between 0 and 1, values closer to 0 depicting similar points and values closer to 1 depicting least similarity. We tried a set of values to check the empirical threshold which maximizes the probability of selecting correct matches at the same time eliminating the noisy matches. We decided this threshold to be **_0.65_** and named it as **_sift_threshold_**, considering any threshold more than that though increases the probability of finding correct matches also increases noisy matches.

After filtering out the SIFT matching descriptor-pairs, we used the RANSAC algorithm to check the geometric arrangement of these descriptor-pairs. The implementation of this algorithm goes with selecting 4 random descriptor-pairs with which the projective transformation (homography) is computed. The first points in a descriptor-pair is transfered and the transformed point is decided to be an inliner if it is in a certain Euclidean distance threshold from the original second point in the pair. This threshold is decided to be **_3_** and named as **_ransac_threshold_**. This procedure is repeated multiple times so as to find the best projective transformation that results in maximum number of inliners. Number of **_iterations_** are fixed as **_2000_** if number of SIFT matching pairs is less than 20, and if it greater than 20 number of iterations is computed as min(num_pairs\*100, 4000).

#### Qulaitative Tests

##### Test_1:
```
./a2 part3 eiffel_18.jpg eiffel_19.jpg
```
Number of Matches = 37  
Drawing Lines for SIFT matching pairs and writing to part3_sift_matches.jpg  
Drawing Lines for matching pairs after RANSAC and writing to part3_ransac_matches.jpg  

**part3_sift_matches_eiffel**
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/part3_sift_matches_eiffel.jpg "part3_sift_matches_eiffel")

**part3_ransac_matches_eiffel**
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/part3_ransac_matches_eiffel.jpg "part3_ransac_matches_eiffel")

##### Test_2:
```
./a2 part3 bigben_2.jpg bigben_3.jpg
```
Number of Matches = 26  
Drawing Lines for SIFT matching pairs and writing to part3_sift_matches.jpg  
Drawing Lines for matching pairs after RANSAC and writing to part3_ransac_matches.jpg  

**part3_sift_matches_bigben**
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/part3_sift_matches_bigben.jpg "part3_sift_matches_bigben")

**part3_ransac_matches_bigben**
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/part3_ransac_matches_bigben.jpg "part3_ransac_matches_bigben")

##### Test_3:
```
./a2 part3 eiffel_18.jpg bigben_3.jpg
```
Number of Matches = 0  
No SIFT matches found..!!

## Part 4: Panorama
#### Running the Code
In order to run this part's code enter the following commands in the project directory
```make
./a2 part4 images/part4/Hill_1.jpg  images/part4/Hill_2.jpg  images/part4/Hill_3.jpg
```
#### Assumptions made
While running the code we assume that input images are given in the order such that the first image has a guaranteed overlap with the second image and the second image has a guaranteed overlap with the third and so on

#### Design Decision
We transform the first image into the prespective of the next image until the very last. For example if we have 3 images the the first image prespective is chaged to second. Then stitch these two together and then the resultant image is transformed to the prespective of the third image. So every panorama is with the prespective of the last image.

#### Design Decision Effects
Since we decided to change the prespective to the very last image, the first images gets transformed multiple times and the left images appear larger in prespective. This can be seen in the results below.



#### Functions and Structs Used
- d_pairs struct to get matching pairs 
- offset struct for finding the offset in pixels
- Sift Matching pairs to find the matching sif features based on SSD
- ransac to apply RANSAC algorithm and get the best 4 pairs for the transformation
- padImage to add padding to images so that they can be properly blended and used in apply_transformation 

#### Changes to apply_transformation from part 1
We had to change the apply transformation form part 1 because we ended up getting images that once we transformed had parts outside the images. We manually calculated the min_x, min_y, max_x, max_y by looping through every transformed pixel and then used these values to offset the image properly to get it inside the frame.
######Note: 
We later realized that we can use the corners of the original image to get the min_x and min_y values rather that looping through every pixel. But we it was to late to make this change to the code

#### Issues Faced
1) Offsetting images 
2) Weird color appearing in the resultant panorama
3) Image Blurring with generated mask
4) Transforming image using inverse warping caused the some part of image to be cut out



### Results

#### Book Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/book_panorama.png)

#### Hill Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/mountain_panorama.png)

#### TV Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/tv_panorama.png)

#### Roof Top(3)  Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/roof_3image_panorama.png)

#### Roof Top(4)  Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/roof_4images_panorama.png)

#### Yosemite Panorama
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A2/yosimite_panorama.png)











