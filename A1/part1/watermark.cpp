//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// PUT YOUR NAMES HERE
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>
#include <cmath>

using namespace std;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input);

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
int check_image(const SDoublePlane &input, int N);

// Pad Image if row count and column count is not equal and are not powers of 2
SDoublePlane padImage(const SDoublePlane &input);

//Test the robustness of watermarking with 100 random numbers
void testRobustness(const SDoublePlane &input, int N);

double min_plane(const SDoublePlane&);
double max_plane(const SDoublePlane&);
double average(const SDoublePlane&);
void NormalizeToGray(SDoublePlane&);


SDoublePlane remove_interference(const SDoublePlane &input)
{
  SDoublePlane new_image = SDoublePlane(input.rows(), input.cols());
  SDoublePlane real = SDoublePlane(input.rows(), input.cols());
  SDoublePlane imag = SDoublePlane(input.rows(), input.cols());
  fft(input, real, imag);

  SDoublePlane s1 = fft_magnitude(real, imag);
  double avs = average(s1);
  double maxs = max_plane(s1);

  SImageIO::write_png_file("uncleaned_spec.png", s1, s1, s1);

  for (int r=0; r < real.rows(); r++){
    for (int c=0; c < real.cols(); c++){
      if ( (abs(r - (real.rows()/2)) < 110) && (abs(r - (real.rows()/2)) > 95) && (s1[r][c] > ( (avs + ( (maxs - avs)/ 3 ) ) ) ) ){
        cout << "Above average pixel at Row: " << r << " Col: " << c << endl;
        imag[r][c] = 0;
        real[r][c] = 0;
      }
    }
  }
  ifft(real, imag, new_image);
  return new_image;
}

double t = 0.3;
int l = 32;
int r = 8;
double alpha = 0.9;

int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }

    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

    SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    cout << "image rows: " << input_image.rows() << endl;
    cout << "image columns: " << input_image.cols() << endl;

    SImageIO::write_png_file(outputFile.c_str(), input_image, input_image, input_image);

    // Check the dimensions and pad the image if required
    if(input_image.rows() != input_image.cols() || (input_image.rows()  & (input_image.rows()-1)))
    	input_image = padImage(input_image);

    if(part == "1.1")
    {
	SDoublePlane fft_real(input_image.rows(), input_image.cols());
    	SDoublePlane fft_imag(input_image.rows(), input_image.cols());
    	fft(input_image, fft_real, fft_imag);
    	SDoublePlane fft_mgnt = fft_magnitude(fft_real, fft_imag);
    	SImageIO::write_png_file(outputFile.c_str(), fft_mgnt, fft_mgnt, fft_mgnt);
    }

    else if(part == "1.2")
    {
      SDoublePlane removed = remove_interference(input_image);
      SImageIO::write_png_file(outputFile.c_str(), removed, removed, removed);
    }
    else if(part == "1.3")
      {
	if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	string op(argv[4]);
	int N = atoi(argv[5]);
  	cout << "N value: " << N << endl;
	if(op == "add")
	{
	  SDoublePlane marked_image = mark_image(input_image, N);
	  SImageIO::write_png_file(outputFile.c_str(), marked_image, marked_image, marked_image);
	}
	else if(op == "check")
	{
	  int result = check_image(input_image, N);
	  cout << "===============\n" << "1.3" << endl;
	  if(result) {
		  cout << "Watermark detected..!" << endl;
	  } else {
		  cout << "Watermark not found..!" << endl;
	  }
	  cout << "===============" << endl;
	}
	else if(op == "quant_test")
	{
	  testRobustness(input_image, N);
	}
	else
	  throw string("Bad operation!");
      }
    else
      throw string("Bad part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}


double average(const SDoublePlane& s)
{
  double sum = 0;
  int count = 0;
  for (int r=0; r < s.rows(); r++){
    for (int c=0; c < s.cols(); c++){
       sum += s[r][c];
       count++;
    }
  }
  return sum/count;
}

double max_plane(const SDoublePlane& s)
{
  double max = -100000000;
  for (int r=0; r < s.rows(); r++){
    for (int c=0; c < s.cols(); c++){
      if (s[r][c] > max)
      {
        max = s[r][c];
      }
    }
  }
  return max;
}

double min_plane(const SDoublePlane& s)
{
  double min = 10000000000;
  for (int r=0; r < s.rows(); r++){
    for (int c=0; c < s.cols(); c++){
      if (s[r][c] < min)
      {
        min = s[r][c];
      }
    }
  }
  return min;
}

int getNearestDimensionToPad(const SDoublePlane &input)
{
	int nearest_2k;
	if((input.rows()  & (input.rows()-1)) || (input.cols()  & (input.cols()-1))) {
		int max_dim = max(input.rows(), input.cols());
		cout << "Max Dim = " << max_dim << endl;
		nearest_2k = pow(2, ceil(log2(max_dim)));
	}
	cout << "Input Rows = " << input.rows() << endl;
	cout << "Input Cols = " << input.cols() << endl;
	cout << "Nearest 2k  = " << nearest_2k << endl;
	return nearest_2k;
}

SDoublePlane padImage(const SDoublePlane &input)
{
	int nearest_2k = getNearestDimensionToPad(input);
	int extra_rows = nearest_2k - input.rows();
	int extra_cols = nearest_2k - input.cols();
	int left_cols, top_rows;
	if(extra_rows > 0){
		top_rows = extra_rows/2;
	}
	if(extra_cols > 0){
		left_cols = extra_cols/2;
	}
	SDoublePlane paddedImage(nearest_2k, nearest_2k);
	double **data = input.row_pointers();
	for(int i = 0; i < input.rows(); i++) {
		memcpy(&paddedImage[top_rows+i][left_cols], data[i], input.cols()*sizeof(double));
	}
	return paddedImage;
}

SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag)
{
	SDoublePlane fft_mgnt(fft_real.rows(), fft_real.cols());
	for(int i = 0; i < fft_real.rows(); i++) {
		for(int j = 0; j < fft_real.cols(); j++) {
			fft_mgnt[i][j] = log(sqrt(pow(fft_real[i][j], 2) + pow(fft_imag[i][j], 2)));
		}
	}
  NormalizeToGray(fft_mgnt);
	return fft_mgnt;
}

void generateRandomBinaryVector(int N, int l, double* vec)
{
	srand(N);
	for(int i = 0; i < l; i++) {
		vec[i] = rand()%2;
	}
}

void moveCircleInSteps(int &x, int &y, char &move, char &dir, int step, int r) {
	if(move == 'x') {
		if(dir == '+') {
			if(x + step <= r) {
				x = x + step;
			} else {
				step = step - abs(x - r);
				x = r;
				move = 'y';
				dir = '-';
				moveCircleInSteps(x, y, move, dir, step, r);
			}
		} else {
			if(abs(x - step) <= r) {
				x = x - step;
			} else {
				step = step - abs(x + r);
				x = -r;
				move = 'y';
				dir = '+';
				moveCircleInSteps(x, y, move, dir, step, r);
			}
		}
	} else {
		if(dir == '+') {
			if(y + step <= r) {
				y = y + step;
			} else {
				step = step - abs(y - r);
				y = r;
				move = 'x';
				dir = '+';
				moveCircleInSteps(x, y, move, dir, step, r);
			}
		} else {
			if(abs(y - step) <= r) {
				y = y - step;
			} else {
				step = step - abs(y + r);
				y = -r;
				move = 'x';
				dir = '-';
				moveCircleInSteps(x, y, move, dir, step, r);
			}
		}
	}
}

void getEquidistantPoints(int r, int l, int* x_coord, int* y_coord)
{
	char move = 'y';
	char dir = '-';	//Start moving clock-wise
	int x = r;	//Start at this position
	int y = 0;	//Start at this position
	int step = 8*r/(2*l);	//Steps to be taken to get next point
	//cout << "Number of points required: " << 2*l << endl;
	//cout << "Step Size required: " << step << endl;
	x_coord[0] = x;
	y_coord[0] = y;
	for(int i = 1; i<2*l; i++) {
		moveCircleInSteps(x, y, move, dir, step, r);
		x_coord[i] = x;
		y_coord[i] = y;
	}
}

struct Point
{
	int x;
	int y;
};

Point getCenterPoint(const SDoublePlane &input) {
	int k = input.rows()/2;
	double maxval = 0;
	Point p;
	for(int i = k-1; i <= k; i++) {
		for(int j = k-1; j <= k; j++) {
			if(input[i][j] > maxval) {
				maxval = input[i][j];
				p.x = i;
				p.y = j;
			}
		}
	}
	return p;
}

SDoublePlane mark_image(const SDoublePlane &input, int N)
{
	SDoublePlane fft_real(input.rows(), input.cols());
	SDoublePlane fft_imag(input.rows(), input.cols());
	fft(input, fft_real, fft_imag);

	double v[l];
	generateRandomBinaryVector(N, l, v);
	int x_coord[2*l];
	int y_coord[2*l];
	getEquidistantPoints(r, l, x_coord, y_coord);
	
	Point p = getCenterPoint(fft_real);
	for(int i = 0; i < l; i++) {
		fft_real[p.x + x_coord[i]][p.y + y_coord[i]] = fft_real[p.x + x_coord[i]][p.y + y_coord[i]] + alpha*fabs(fft_real[p.x + x_coord[i]][p.y + y_coord[i]])*v[i];
		fft_real[p.x + x_coord[i+l]][p.y + y_coord[i+l]] = fft_real[p.x + x_coord[i+l]][p.y + y_coord[i+l]] + alpha*fabs(fft_real[p.x + x_coord[i+l]][p.y + y_coord[i+l]])*v[i];
	}
	SDoublePlane output(input.rows(), input.cols());
	ifft(fft_real, fft_imag, output);
	return output;
}

double mean(double* vec, int length)
{
	double sum = 0;
	for(int i = 0; i < length; i++) {
		sum += vec[i];
	}
	return sum/length;
}

double sd(double* vec, int length)
{
	double sum = 0;
	double mu = mean(vec, length);
	for(int i = 0; i < length; i++) {
		sum += pow((vec[i] - mu), 2);
	}
	return sqrt(sum/(length-1));
}

double sumOfProducts(double* v1, double* v2, int length)
{
	double sum = 0;
	for(int i = 0; i < length; i++) {
		sum += v1[i]*v2[i];
	}
	return sum;
}

double pearsonsCoefficient(double* v1, double* v2, int length)
{
	double mu1 = mean(v1, length);
	double mu2 = mean(v2, length);
	double sd1 = sd(v1, length);
	double sd2 = sd(v2, length);
	double sop = sumOfProducts(v1, v2, length);
	double pearson_cor_coeff = (sop - length*mu1*mu2)/((length-1)*sd1*sd2);
	//cout << "Pearson Correlation Coefficient = " << pearson_cor_coeff << endl;
	return pearson_cor_coeff;
}

int check_image(const SDoublePlane &input, int N)
{

	int result = 0;
	SDoublePlane fft_real(input.rows(), input.cols());
	SDoublePlane fft_imag(input.rows(), input.cols());
	fft(input, fft_real, fft_imag);
	double v[l];
	generateRandomBinaryVector(N, l, v);
	int x_coord[2*l];
	int y_coord[2*l];
	getEquidistantPoints(r, l, x_coord, y_coord);
	Point p = getCenterPoint(fft_real);
	double bins[l];
	for(int i = 0; i < l; i++) {
		bins[i] = fft_real[p.x + x_coord[i]][p.y + y_coord[i]];
	}
  double pearson_cor_coeff = pearsonsCoefficient(v, bins, l);
	if(fabs(pearson_cor_coeff) > t)
  {
		result = 1;
	}
	return result;
}

double dist(double x, double y, int size)
{
  if (size == 0)
  {
    return 0;
  }
  return abs(x - size/2) + abs(y - size/2);
}

void NormalizeToGray(SDoublePlane& plane)
{
  double min = min_plane(plane);
  double max = max_plane(plane);
  for (int r=0; r < plane.rows(); r++){
    for (int c=0; c < plane.cols(); c++){
      plane[r][c] = (255 * (plane[r][c] - min))/(max - min);
    }
  }
}

void testRobustness(const SDoublePlane &input, int N)
{
	SDoublePlane marked_imag = mark_image(input, N);
	int result_orig = check_image(input, N);
	int result_marked = check_image(marked_imag, N);
	int false_pos = 0;
	for(int n = 207; n < 307; n++) {
		 if(check_image(marked_imag, N*n)) {
			 false_pos += 1;
		 }
	}
	if(result_marked) cout << "Watermark Survived..!" << endl;
	else cout << "Watermark not found..!" << endl;
	cout << "False Positives = " << false_pos << endl;

}

