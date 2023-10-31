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
void projectPoint(pPoint* pointStruct, double K);
void computePoint(double R1, double R2, double theta, double phi, double z0, pPoint* pointStruct);
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
  pPoint* pointStruct = new pPoint;
  double theta_curr = 0;
  double rotations[8] = {0, M_PI_4, M_PI_2, 3*M_PI_4, M_PI, M_PI+M_PI_4, M_PI+M_PI_2, M_PI+3*M_PI_4};
  for(int i = 0; i < 8; i++)
  {
    while(theta_curr < 2*M_PI)
    {
      computePoint(R1, R2, theta_curr, rotations[i], z0, pointStruct);
      double* point = pointStruct->origPoint;
      std::cout << '(' << point[0] << ',' << point[1] << ',' << point[2] << ')' << std::endl;
      theta_curr += theta_inc;
    }
    theta_curr = 0;
  }
}

void populateFrameBuffer(std::tuple<char, double>** FrameBuffer)
{
  pPoint* pointStruct = new pPoint;
  double theta_curr = 0; 
  double phi_curr = 0;
  while(phi_curr < 2*M_PI)
  {
    theta_curr = 0;
    while (theta_curr < 2*M_PI)
    {
      computePoint(R1, R2, theta_curr, phi_curr, z0, pointStruct);
      projectPoint(pointStruct, K);
      //std::cout << '(' << orig_point[0] << ',' << orig_point[1] << ',' << orig_point[2] << ')' << std::endl;
      //std::cout << '(' << proj_point[0] << ',' << proj_point[1] << ',' << '0' << ')' << std::endl;
      int x_cord = (screenWidth/2) + ((int) ((pointStruct->projPoint)[0]));
      int y_cord = (screenHeight/2) + ((int) ((pointStruct->projPoint)[1]));
      if(x_cord < 0 || x_cord >= screenWidth || y_cord < 0 || y_cord >= screenHeight)
      {
        theta_curr += theta_inc;
        continue;
      }
      std::tuple<char,double> curr_buff = FrameBuffer[x_cord][y_cord];
      double z_buff = std::get<double>(curr_buff);
      if (z_buff == 0 || (((pointStruct->projPoint)[2]) < z_buff))
      {
        FrameBuffer[x_cord][y_cord] = std::make_tuple('@', ((pointStruct->projPoint)[2]));
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

void projectPoint(pPoint* pointStruct, double K)
{
  double z_inv = 1/((pointStruct->origPoint)[2]);
  double x_prime = K*z_inv*((pointStruct->origPoint)[0]);
  double y_prime = K*z_inv*((pointStruct->origPoint)[1]);
  double* pPoint = new double[3];
  pPoint[0] = x_prime;
  pPoint[1] = y_prime;
  pPoint[2] = ((pointStruct->origPoint)[2]);
  pointStruct->projPoint = pPoint;
}

void computePoint(double R1, double R2, double theta, double phi, double z0, pPoint* pointStruct)
{
  double *projPoint, *origPoint, *lumPointOrig, *lumPointProj;
  double** rotMat = makeZRotMat(phi);
  origPoint = new double[3];
  origPoint[0] = R1 + R2*cos(theta);
  origPoint[1] = 0;
  origPoint[2] = z0+R2*sin(theta);
  lumPointOrig = new double[3];
  lumPointOrig[0] = cos(theta);
  lumPointOrig[1] = 0;
  lumPointOrig[2] = sin(theta);
  projPoint = squareMatVec(rotMat, origPoint, 3);
  lumPointProj = squareMatVec(rotMat, lumPointOrig, 3);
  delete origPoint;
  delete rotMat[0];
  delete rotMat[1];
  delete rotMat[2];
  delete rotMat;
  delete lumPointOrig;
  pointStruct->origPoint = projPoint;
  pointStruct->lumVec = lumPointProj;
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