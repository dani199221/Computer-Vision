//
// detect.cpp : Detect integrated circuits in printed circuit board (PCB) images.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// PUT YOUR NAMES HERE
//
//

#include "SImage.h"
#include "SImageIO.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};


//Template matching to find the IC in the image
//https://www.wikiwand.com/en/Template_matching
// input---> image, filter
// output---> None
// A global variable is maintained to  detect the box 
// 3354631 total for ic_1_copy2
vector <DetectedBox> boxes;
vector<DetectedBox>result;

void template_matching(const SDoublePlane &image, SDoublePlane &filter){
    
    double min_sad  = 1.39587e+06; 
    for(int x = 0;x<image.rows()-filter.rows();x++){
        for(int y = 0; y<image.cols()-filter.cols(); y++){
            double sad = 0;
            for(int i=0; i<filter.rows(); i++){
                for(int j=0; j< filter.cols(); j++){
                    double p_image = image[x+i][y+j];
                    double p_filter = filter[i][j];

                    sad += abs(p_image -p_filter);
                }
            }

            if (min_sad > sad){
                
                DetectedBox d;
                d.row = x, d.col = y;  
                d.width = filter.cols(), d.height = filter.rows();
                d.confidence = 1 -sad/3354631.0;
                bool will_push = true;
                for(vector<DetectedBox>::iterator it = boxes.begin(); it < boxes.end(); it++){
                    if (abs(it->row -d.row) < filter.rows() +20 && abs(it->col - d.col) < filter.cols() +20){
                       will_push = false; 
                       if (it->confidence < d.confidence){
                           //boxes.erase(it);
                           it->row = d.row, it->col = d.col, it->confidence = d.confidence;
                           //will_push = true;
                       }
                       break;
                    }
                }
                if (will_push){
                    boxes.push_back(d);
                }
                    
            }
        }
    }
    for(unsigned i = 0; i < boxes.size(); i++){
        bool will_push = true;
        for(unsigned j = 0; j < boxes.size(); j++){
            if (i != j){
                if(abs(boxes[i].row - boxes[j].row) < filter.rows() + 20 && abs(boxes[i].col - boxes[j].col) < filter.cols()+ 20){
                   if(boxes[i].confidence < boxes[j].confidence){
                       will_push = false;
                       break;
                   }
                }
            }
        }
        if (will_push){
            result.push_back(boxes[i]);
        }
    }
                                        
                    

}

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &ics)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &ics, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

void print_kernel(const SDoublePlane &filter){
    for(int i=0;i<filter.rows();i++){
        for(int j=0;j<filter.cols();j++){
            cout<< filter[i][j] << "      ";
        }
        cout << endl;
    }
}

SDoublePlane flip_kernel(const SDoublePlane &input){
    SDoublePlane output(input.rows(), input.cols());
    int row = 0, col = 0;
    for(int i = input.rows() - 1; i >= 0; i--){
        col = 0;
        for(int j = input.cols() - 1; j >= 0; j--){
            output[row][col] = input[i][j];
            col++;
        }
        row++;
    }             
    return output;
}

// Check if the index exists for the convolution
bool indexExists(const SDoublePlane &input, const int row, const int col){
  // cout << row << ", " << col << "\n"; 
  return row < input.rows() && row >= 0 && col < input.cols() && col >= 0;
}



// Takes a 1d even filter and extends it to be odd
// fills extra space with 0
// dimension: false = row vector, true = col vector
SDoublePlane makeOdd1D(const SDoublePlane &filter1D, bool dimension){
  SDoublePlane result;

  // col vector
  if(filter1D.rows() % 2 == 0 && dimension){
    result = SDoublePlane(filter1D.rows() + 1, 1);
  }

  // row vector
  if(filter1D.cols() % 2 == 0 && !dimension){
    result = SDoublePlane(1, filter1D.cols() + 1);
  }

  // already odd
  if(filter1D.rows() % 2 == 1 && filter1D.cols() % 2 == 1){
    result = SDoublePlane(filter1D.rows(), filter1D.cols());
  }

  // copy values
  for(int i = 0; i < filter1D.rows(); i++){
    for(int j = 0; j < filter1D.cols(); j++){
      result[i][j] = filter1D[i][j];
    }
  }

  return result;
}

// Convolve an image with a separable convolution kernel
// ASSUMES FILTERS HAVE ALREADY BEEN REFLECTED
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{



  // get odd length row and col vectors
  SDoublePlane odd_row_filter = makeOdd1D(row_filter, true);
  SDoublePlane odd_col_filter = makeOdd1D(col_filter, true);

  SDoublePlane intermediate(input.rows(), input.cols());
  SDoublePlane output(input.rows(), input.cols());
  // Convolution code here

  // find centers of row and col filters
  int rowCenterIndex = odd_row_filter.cols() / 2;
  int colCenterIndex = col_filter.rows() / 2;

  // apply row filter
  for(int row = 0; row < input.rows(); row++){
    for(int col = 0; col < input.cols(); col++){
      // cout << row << ", " << col << "      " << input.rows() << ", " << input.cols() <<"\n"; 
      // add center
      double rowVal = input[row][col] * odd_row_filter[0][rowCenterIndex];

      // loop through rest of row kernel
      for(int i = 1; i <= rowCenterIndex; i++){
        // add upper
        if(indexExists(input, row, col - i)){
          rowVal += input[row][col - i] * odd_row_filter[0][rowCenterIndex - i];
          // cout << row_filter[0][rowCenterIndex - i] << " " << rowCenterIndex - i << " (" << row << " ," << col << ")\n";
        }

        // add lower
        if(indexExists(input, row, col + i)){
          rowVal += input[row][col + i] * odd_row_filter[0][rowCenterIndex + i];
        }
      }

      intermediate[row][col] = rowVal;
    }
  }

  // apply col filter
  for(int row = 0; row < intermediate.rows(); row++){
    for(int col = 0; col < intermediate.cols(); col++){
      // add center
      double colVal = intermediate[row][col] * col_filter[colCenterIndex][0];

      // loop through rest of col kernel
      for(int i = 1; i <= colCenterIndex; i++){
        // add left
        if(indexExists(intermediate, row - i, col)){
          colVal += intermediate[row - i][col] * col_filter[colCenterIndex - i][0];
        }

        // add right
        if(indexExists(intermediate, row + i, col)){
          colVal += intermediate[row + i][col] * col_filter[colCenterIndex + i][0];
        }
      }

      output[row][col] += colVal;
     
    }
  }
    
  return output;
}

//check if the kernel is odd if not, then add row or col to make it odd
SDoublePlane make_odd(const SDoublePlane &input){
    int rows = input.rows(), col = input.cols();
    if (input.rows() % 2 == 0){
        rows++;
    }
    if (input.cols() % 2 == 0){
        col++;
    }

    
    if(rows > input.rows() || col > input.cols()){

        SDoublePlane res(rows, col);
        for(int i = 0; i< input.rows();i++){
            for (int j=0; j< input.cols(); j++){
                res[i][j] = input[i][j];
            }
        }
        return res;

    }else{
        return input;
    }
}
   
// Convolve an image with a  convolution kernel
// pseudo code for convolution
//https://www.wikiwand.com/en/Kernel_(image_processing)
//Assuming the filter is flipped before calling the function
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter){
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane filter2 = make_odd(filter);

    int max = INT_MIN;
    int min = INT_MAX;

    for(int i = 0; i < input.rows(); i++){     // iterate over all the rows of the image
        for (int j = 0; j< input.cols(); j++){ // iterate over all the cols of the image
            for (int k =-filter2.rows()/2; k <= filter2.rows()/2 ; k++){         // iterate over all the rows of the filter 
                for (int l = -filter2.cols()/2; l <=filter2.cols()/2; l++){      // iterate over all the cols of the filter
                    if (indexExists(input, i+k, j+l))
                        output[i][j] += input[i+k][j+l] * filter2[k+filter2.rows()/2][l+filter2.cols()/2]; //convolution process
                }
            }
            if (output[i][j] > max){
                max = output[i][j];
            }else if(output[i][j]< min){
                min = output[i][j];
            }
        }
    }
    //bring the values betweem 0-255
    for(int i = 0 ;i< input.rows(); i++){
        for(int j = 0; j< input.cols(); j++){
            output[i][j] = ((output[i][j] - min) * 255)/ (max - min);
        }
    }
  
    return output;
}


// Apply a sobel operator to an image, returns the result
//https://www.wikiwand.com/en/Sobel_operator 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane output_x(input.rows(), input.cols());
    SDoublePlane output_y(input.rows(), input.cols());
     
	if (_gx){	
		// Implement a sobel gradient estimation filter with 1-d filters
		
		SDoublePlane row_f(1,3); 
		row_f[0][0] = 1,row_f[0][1] = 2,row_f[0][2] = 1;
		SDoublePlane col_f(3,1); 
		col_f[0][0] = -1,col_f[1][0] = 0,col_f[2][0] = 1;
		output_y = convolve_separable(input, row_f, col_f);

		row_f[0][0] = 1,row_f[0][1] = 0,row_f[0][2] = -1;
		col_f[0][0] = 1,col_f[1][0] = 2,col_f[2][0] = 1;
		output_x = convolve_separable(input, row_f, col_f);
	
	}else{	
		//Implement a sobel gradient with 2-d filters

		SDoublePlane sobi(3,3);
		sobi[0][0] = 1,sobi[0][1] = 0,sobi[0][2] = -1; 
		sobi[1][0] = 2,sobi[1][1] = 0,sobi[1][2] = -2; 
		sobi[2][0] = 1,sobi[2][1] = 0,sobi[2][2] = -1; 
		output_y = convolve_general(input, sobi);

		sobi[0][0] = 1,sobi[0][1] = 2,sobi[0][2] = 1; 
		sobi[1][0] = 0,sobi[1][1] = 0,sobi[1][2] = 0; 
		sobi[2][0] = -1,sobi[2][1] =-2,sobi[2][2] =-1; 
		output_x = convolve_general(input, sobi);


    }
    //get average of both horizontal and vertical edges in the final image
    for(int i=0; i<output.rows(); i++){
        for (int j=0; j<output.cols(); j++){
            output[i][j] = sqrt( pow(output_x[i][j],2) + pow(output_y[i][j],2));
        }
    }

	return output;

}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
    SDoublePlane input2(input.rows(), input.cols());
    SDoublePlane ic= SImageIO::read_png_file("template1.png");

    input2 = sobel_gradient_filter(input, true);

    
    //threshold the image
    for(int i = 0; i < input2.rows(); i++){
        for (int j = 0; j< input2.cols(); j++){
            if (input2[i][j] > thresh){
                input2[i][j] = 255;
            }else{
                input2[i][j] = 0;
            }
        }
    }
    
    //template matching to find the ic in the image
    template_matching(input2, ic); 
  
  return input2;
}


//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
 
    SDoublePlane rowFilter(1, 3);
    rowFilter[0][0] = .25;
    rowFilter[0][1] = .5;
    rowFilter[0][2] = .25;

    SDoublePlane colFilter(3, 1);
    colFilter[0][0] = .25;
    colFilter[1][0] = .5;
    colFilter[2][0] = .25;

  //convloving with a gaussian filter to get rid of the pcb trace lines 
  input_image = convolve_separable(input_image, rowFilter, colFilter);
  // finding edges in the image
  SDoublePlane output_image = find_edges(input_image, 225);

  //vector to hold the boxes to be drawn
  vector<DetectedBox> ics = result;
  vector<DetectedBox> empty_ics;

   for(vector<DetectedBox>::iterator it = ics.begin(); it != ics.end(); it++)
        cout<< it->confidence<< endl;
     write_detection_txt("detected.txt", ics); 
    write_detection_image("detected.png", ics, input_image);
    write_detection_image("edges.png", empty_ics, output_image);

}
