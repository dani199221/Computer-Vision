## Part 1

### 1.1
#### Running the Code
Run make to compile the program and then input the command
```
./watermark 1.1 input_image.png out_image.png
```
with input_image as the image to create the spectogram from, with out_image being the name of the resulting png file.

#### Overview
In order to make the output png clear, we also normalize the magnitude of the image to have values between 0 and 255.

You can see a sample spectogram of Yoshi below.

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/yoshi.png "yoshi")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/yoshi-spec.png "yoshi spectogram")
### 1.2
#### Running the Code
Run make to compile the program and then input the command
```
./watermark 1.2 input_image.png out_image.png
```
#### Overview
In order to remove the noise we set the values of the imaginary and complex values that were out of the ordinary to zero,
these corresponded the diagonal noise in the image.  They were the only outliers with above average values at such high frequencies,
so a simple conditional singled them out to be zeroed.

The results can be seen below.

![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/noise1.png "noisy")
![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/cleaned.png "cleaned")

### 1.3
#### Running the Code
Run make to compile the program and then input the command

To add watermark:
```
./watermark 1.3 Lenna.png marked_Lenna.png add 333
```

To check watermark:
```
./watermark 1.3 marked_Lenna.png marked_Lenna.png check 333
```

To test the robustness:
```
./watermark 1.3 Lenna.png marked_Lenna.png quant_test 333
```


#### Overview
We are injecting the watermark by generating a random binary vector by seeding a given N, which makes surte that the results are reproducable and also to fetch the same binary vector to check an elready embedded watermark. To embed the water mark first the image is converted to its fourier domain, and then the binary vector is used to update the selected bins in the real portion of DFT.

The bins selected are equidistant from each other across a circle of radius *r* and every bin has a corresponding bin selected which is symmetric about the origin. So, for a fixed *l* (length of the binary vector) we have updated *2l* bins using the formula R'(x,y) = R(x,y) + alpha\*|R(x,y)|\*v<sub>i</sub>, for a selected *alpha*. As the updated bin values are a linear combination of the original bin values and the embedded binary vector, we decided to use either [Pearson Correlation Coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) or [Point-Biserial Coefficient](https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient) (which is better when calculating correlation between a continuous vector and a dichotomous vector). Both these metrics gave almost similar results, so we went ahead with Pearson Correlation Coefficient. Below are the optimum values used for the parameters to handle this task.

- *r* = 8
- *l* = 32
- *alpha* = 0.9
- *t* = 0.3

#### Parameter Selection and Issues:
- As *r* values increases the bins are selected from a sparse area and is very difficult to detect the modified values after embedding. At the same time if *r* is very small it limited the number of bins (*l*) that we select because the points on these circles are finite and the maximum *l* that can be selected reduces with *r*. So we searched on a grid of values for different images and fixed the *r* value to be 8.
- Similarly small values for *l* made it difficult to find the embedded vectors which increased False negatives. Increasing *l* increased the probability of finding the embedded vectors, but a higher *l* had a lot of False Positives.
- Selecting *alpha* turned out to be very interesting, because having higher values gave very good results both in terms of accuracy and reducing False Positives, but higher values of *alpha* (eg. 2, 3 etc.) introduced noticeable distortions onto the image. So we have fixed *alpha* value to be 0.9 as a tradeoff between the accuracy and image distortion.
- Correlation coefficient threshold is the easiest one to select as a lot of insights were drawn with all the earlier results while selecting remaining parameters. Having *t* = 0.4 along with higher alpha gave tremendous results, but as we had to lower *alpha* value we tried a lower treshold and fixed to 0.3. Much lower values of the threshold increased both False Negatives and False Positives.

#### Qualitative Test:

**Before Embedding watermark**
![](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/Lenna.png) 

**After Embedding watermark**
![](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/marked_Lenna.png)

```
./watermark 1.3 Lenna.png marked_Lenna.png add 333

./watermark 1.3 marked_Lenna.png marked_Lenna.png check 333
```

In: marked_Lenna.png  Out: Lenna.png  
image rows: 512  
image columns: 512  
N value: 333  
1.3  
Watermark detected..!  


**Before Embedding watermark**
![](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/cleaned.png) 

**After Embedding watermark**
![](https://github.com/dani199221/Computer-Vision/blob/master/A1/part1/marked_cleaned.png)

```
./watermark 1.3 cleaned.png marked_cleaned.png add 333

./watermark 1.3 marked_cleaned.png marked_cleaned.png check 333
```

In: marked_cleaned.png  Out: marked_cleaned.png  
image rows: 512  
image columns: 512  
N value: 333  
1.3  
Watermark detected..!  


#### Quantitative Test:
To check the robustness of this water mark detection we embedded a water mark corresponding to a given number N, and checked whether check_image function can detect this watermark. Later we took 100 random N's and cheked whether check_image function detects any of these 100 water marks, if detected they are counted as False Positives.

```
./watermark 1.3 Lenna.png marked_Lenna.png quant_test 333
```
Watermark Survived..!
False Positves = 13

```
./watermark 1.3 cleaned.png marked_cleaned.png quant_test 333
```
Watermark Survived..!
False Positives = 16

## Part2:

### HOW TO RUN: 
```
make 
./detect image.png
```

### DESIGN DECISIONS:

We are using template matching to detect squares in our image. We cropped out ICs from the original images to use them as templates for detection. We found out that the ic with pins on all 4 sides can act as pretty good detectors for IC in our images so we opted to use this as our primary template for detection.  

1) Convolution General:
    We have four nested loop that loop through the whole image and the kernel and convolve the image and save the result in the output. While convolving the edges we considered the indexes out of boundary to be zero. To cater for even kernels we wrote a function to pad the kernel with zeros to make it odd length so that we have a center for convolution. The function assumes that the kernel is flipped. We have also added an additional function to flip the kernel before it can be sent to the convolution_general function. 

2) Convolution Separable:
    We loop through the image twice, once with the row filter and the result is convolved with the column filter. We handle the boundary conditions the same as in the convolution_general function. To handle the even filters we pad with zeroes to get to an odd length for both row and column filters   


3) Sobel_filter
   Implements the sobel 3x3 filter. The _gx bool variable is used to swap between separable(_gx = True) and general convolution(_gx = False). 
   https://www.wikiwand.com/en/Sobel_operator. The wiki page was consulted to see the sobel filter values.

4) Edge_detect
   We smooth the input_image with a gaussian filter and pass it to the edge_detect function which is then convolved with the sobel filter to detect edges. We then threshold the image and our template used for template matching and pass it to the template_matching function which looks for all the matching parts in the image and returns the result. We also thresholded our image and filter such that everything above 225 is set to 255 and everything below it is set to 0.

5) Template_matching
   The function takes as input the edge detected image and the edge detected filter and looks for the filter throught the complete image and tries to find all similar  occurrances of the template in the image. Template matching is done by finding sum of absoulute differences in pixel intensities in our template and a sliding window on our image.We threshold these results at 1.39587e+06. This value was found by testing each gievn image with our template and finding the max of all the minimum  differences that detected some IC in the image. We then sum all the values from our template and divide the absolute difference with this and subtracted from 1 to get the confidence level of the detected IC. The detected edges are stored in a global variable result which is them used im main function for plotting squares.
    To get the best window for the ic, we compare the confidencce level of other windows around it and chose the one with the highest confidence. 

    Input:  Edge detected image (SDoubleplane) , Edge detected filter (Sdoubleplane)
    Output: none
    https://www.wikiwand.com/en/Template_matching. The wiki page was used to understand template matching.
### Accuracy 
Since we are using a constant size window, the window wont ever change size to perfectly match the IC. If there is a small Ic in the image, it might not be detected at all and if the IC in the image is very large then it might detect quadrants of the IC as separate ICs. However if the size of the IC in the image is close to our template the match is generally good. There were cases where the algorithm ended up detecting ports as IC. The code is relatively fast on smaller images however it ends up taking plenty of time on images of larger size. It took around 10mins for ic_6.png to be completed. 

Python script results:
    <p>ic_1.png Mean average precision for ic = 0.750000</p>
    ![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part2/detected_ic1.png "ic1")
    <p>ic_2.png Mean average precision for ic = 1.000000</p>
    ![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part2/detected_ic2.png "ic2")
    <p>ic_3.png Mean average precision for ic = 0.000000</p>
    ![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part2/detected_ic3.png "ic3")
    <p>ic_4.png Mean average precision for ic = 0.713333</p> 
    ![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part2/detected_ic4.png "ic4")
    <p>ic_5.png Mean average precision for ic = 0.000000</p>
    ![alt text](https://github.com/dani199221/Computer-Vision/blob/master/A1/part2/detected_ic5.png "ic5")

### How to Improve 
1) The size of the template can be varied in accordance with the image or we could run through the image with multiple templates detect different size ICs. 
2) Faster Convolution using fourier transform instead of naive approach can be used to speed up the code.

Other Approaches Cosidered. We also thought of transforming parts of our image using a sliding window and the kernel to the fourier domain and comparaing the distances to find minumums which would be matches.

