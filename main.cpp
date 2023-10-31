#include <cmath>
#include <tuple>
#include <iostream>
#include <string>

/*
  BIG IDEA:
  (a) Computer (x,y,z) point based on parameters
  (b) Compute associated luminance value for (x,y,z) coordinate based on light vector (xl, yl)
  (c) Project (x,y,z) coordinate along w/ luminance onto x',y' plane using zBuffer
  (d) Draw out coordinate using ASCII
*/

struct pPoint {
  double* origPoint;
  double* lumVec;
  double* projPoint; 
}; 

int screenWidth = 50;
int screenHeight = 50;
double K = 7;
double R1 = 20; 
double R2 = 5;
double z0 = 10;
double theta_inc = 0.07;
double phi_inc = 0.07;

void debugPlot();
void populateFrameBuffer(std::tuple<char, double>** FrameBuffer);
void clearFrameBuffer(std::tuple<char, double>** FrameBuffer);
double* projectPoint(double* origPoint, double K);
double* computePoint(double R1, double R2, double theta, double phi, double z0);
double** makeXRotMat(double phi);
double** makeYRotMat(double phi);
double** makeZRotMat(double phi);
double* squareMatVec(double** Mat, double* Vec, int dim);




int main(int argc, char *argv[])
{
  std::tuple<char, double>** FrameBuffer = new std::tuple<char, double>*[screenHeight];
  for (int i = 0; i < screenHeight; i++) FrameBuffer[i] = new std::tuple<char, double>[screenWidth];
  clearFrameBuffer(FrameBuffer);
  populateFrameBuffer(FrameBuffer);
  for (int i = 0; i < screenHeight; i++)
  {
    for (int j = 0; j < screenWidth; j++)
    {
      std::cout << std::get<char>(FrameBuffer[i][j]);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}

void debugPlot()
{
  double theta_curr = 0;
  double rotations[8] = {0, M_PI_4, M_PI_2, 3*M_PI_4, M_PI, M_PI+M_PI_4, M_PI+M_PI_2, M_PI+3*M_PI_4};
  for(int i = 0; i < 8; i++)
  {
    while(theta_curr < 2*M_PI)
    {
      double* point = computePoint(R1, R2, theta_curr, rotations[i], z0);
      std::cout << '(' << point[0] << ',' << point[1] << ',' << point[2] << ')' << std::endl;
      theta_curr += theta_inc;
    }
    theta_curr = 0;
  }
}

void populateFrameBuffer(std::tuple<char, double>** FrameBuffer)
{
  double theta_curr = 0; 
  double phi_curr = 0;
  while(phi_curr < 2*M_PI)
  {
    theta_curr = 0;
    while (theta_curr < 2*M_PI)
    {
      double* orig_point = computePoint(R1, R2, theta_curr, phi_curr, z0);
      double* proj_point = projectPoint(orig_point, K);
      double* p_point = computePoint(R1, R2, theta_curr, M_PI, z0);
      //std::cout << '(' << orig_point[0] << ',' << orig_point[1] << ',' << orig_point[2] << ')' << std::endl;
      //std::cout << '(' << proj_point[0] << ',' << proj_point[1] << ',' << '0' << ')' << std::endl;
      int x_cord = (screenWidth/2) + ((int) proj_point[0]);
      int y_cord = (screenHeight/2) + ((int) proj_point[1]);
      if(x_cord < 0 || x_cord >= screenWidth || y_cord < 0 || y_cord >= screenHeight)
      {
        theta_curr += theta_inc;
        continue;
      }
      std::tuple<char,double> curr_buff = FrameBuffer[x_cord][y_cord];
      double z_buff = std::get<double>(curr_buff);
      if (z_buff == 0 || (proj_point[2] < z_buff))
      {
        FrameBuffer[x_cord][y_cord] = std::make_tuple('@', proj_point[2]);
      }
      theta_curr += theta_inc;
    }
    phi_curr += phi_inc;
  }
}

void clearFrameBuffer(std::tuple<char, double>** FrameBuffer)
{
  for (int i = 0; i < screenHeight; i++)
  {
    for (int j = 0; j < screenWidth; j++)
    {
      FrameBuffer[i][j] = std::make_tuple(' ', 0);
    }
  }
}

double* projectPoint(double* origPoint, double K)
{
  double z_inv = 1/origPoint[2];
  double x_prime = K*z_inv*origPoint[0];
  double y_prime = K*z_inv*origPoint[1];
  double* p_point = new double[3];
  p_point[0] = x_prime;
  p_point[1] = y_prime;
  p_point[2] = origPoint[2];
  return p_point;
}

double* computePoint(double R1, double R2, double theta, double phi, double z0)
{
  double *projPoint, *origPoint;
  double** rotMat = makeZRotMat(phi);
  origPoint = new double[3];
  origPoint[0] = R1 + R2*cos(theta);
  origPoint[1] = 0;
  origPoint[2] = z0+R2*sin(theta);
  projPoint = squareMatVec(rotMat, origPoint, 3);
  delete origPoint;
  delete rotMat[0];
  delete rotMat[1];
  delete rotMat[2];
  delete rotMat;
  return projPoint;
}

double** makeXRotMat(double phi)
{
  double** rotMat = new double*[3];
  for (int i = 0; i < 3; i++) rotMat[i] = new double[3];
  rotMat[0][0] = 1;
  rotMat[0][1] = 0;
  rotMat[0][2] = 0;
  rotMat[1][0] = 0;
  rotMat[1][1] = cos(phi);
  rotMat[1][2] = -sin(phi);
  rotMat[2][0] = 0;
  rotMat[2][1] = sin(phi); 
  rotMat[2][2] = cos(phi);
  return rotMat;

}

double** makeYRotMat(double phi)
{
  double** rotMat = new double*[3];
  for (int i = 0; i < 3; i++) rotMat[i] = new double[3];
  rotMat[0][0] = cos(phi);
  rotMat[0][1] = 0;
  rotMat[0][2] = sin(phi);
  rotMat[1][0] = 0;
  rotMat[1][1] = 1;
  rotMat[1][2] = 0;
  rotMat[2][0] = -sin(phi);
  rotMat[2][1] = 0; 
  rotMat[2][2] = cos(phi);
  return rotMat;
}

double** makeZRotMat(double phi)
{
  double** rotMat = new double*[3];
  for (int i = 0; i < 3; i++) rotMat[i] = new double[3];
  rotMat[0][0] = cos(phi);
  rotMat[0][1] = -sin(phi);
  rotMat[0][2] = 0;
  rotMat[1][0] = sin(phi);
  rotMat[1][1] = cos(phi);
  rotMat[1][2] = 0;
  rotMat[2][0] = 0;
  rotMat[2][1] = 0; 
  rotMat[2][2] = 1;
  return rotMat;

}

// Perform a MatVec operation on dim*dim Mat and appropriate dimensioned vector
double* squareMatVec(double** Mat, double* Vec, int dim)
{
  double* newVec = new double[dim];
  for (int i = 0; i < dim; i++)
  {
    newVec[i] = 0;
    for (int j = 0; j < dim; j++)
    {
      newVec[i] += Mat[i][j] * Vec[j];
    }
  }
  return newVec;
}