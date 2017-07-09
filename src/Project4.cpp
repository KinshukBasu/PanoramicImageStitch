#include <iostream>
#include <cstdlib>
#include <vector>
#include "ifs.h"
#include "flip.h"
#include "math.h"

using namespace std;
int cmax_result, rmax_result;
IFSIMG Leftimg, Rightimg,Result,ifs_Ldx, ifs_Ldy,ifs_Rdx,ifs_Rdy;
float **leftimg, **rightimg, **finalimg;
float **Ldx, **Ldy, **Rdx, **Rdy;
float **H;
float **Hinv;


int checkValid(int x,int y,int i,int j)
{
	int maxrow, maxcol;

	maxcol = ifsdimen(Leftimg,0);			
	maxrow = ifsdimen(Leftimg,1);

	if(y+i<0 || y+i>=maxrow)
		return(0);
	else if(x+j<0 || x+j>=maxcol)
		return(0);

	return(1);

}

void buildTemplate(int L_or_R, int imgx,int imgy, float imgtemplate[][5])
{
	//float imgtemplate[5][5];
	int i,j;

	for(i=-2;i<=2;i++)
	{
		for(j=-2;j<=2;j++)
		{
			if(checkValid(imgx,imgy,i,j))
			{
				if(L_or_R==0)
				{
					imgtemplate[2+i][2+j] = sqrt(Ldx[imgy+i][imgx+j]*Ldx[imgy+i][imgx+j] + Ldy[imgy+i][imgx+j]*Ldy[imgy+i][imgx+j]);
				}
				else if(L_or_R==1)
				{
					imgtemplate[2+i][2+j] = sqrt(Rdx[imgy+i][imgx+j]*Rdx[imgy+i][imgx+j] + Rdy[imgy+i][imgx+j]*Rdy[imgy+i][imgx+j]);
				}
			}
			else
			{
				imgtemplate[2+i][2+j]=0;
			}
		}
	}

	return;
}


float diffMeasure(float a[][5], float b[][5])
{
	int i,j;
	float diff=0.0;

	for(i=0;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			diff += sqrt((a[i][j] - b[i][j])*(a[i][j] - b[i][j]));
		}
	}
	return(diff);
}


std::vector<int> detectCorner (IFSIMG img)
{
	std::vector<int> results;
	int *len;
  	len = ifssiz(img);
	int maxcol = len[1];
	int maxrow = len[2];

	cout << "Size of the image is: "<< endl << "Row: " << maxrow << " Col: " << maxcol << endl;

	IFSIMG Ix, Iy, IxIy;
	Ix = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	Iy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	IxIy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);


	fldx(img, Ix);
	fldy(img, Iy);
	fldx(Iy, IxIy);


	IFSIMG Ixx, Iyy, detect;
	Ixx = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	Iyy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);

	fldx(Ix, Ixx);
	fldy(Iy, Iyy);


	ifsfree(Ix, IFS_FR_ALL);
	ifsfree(Iy, IFS_FR_ALL);

	IFSIMG Sxx, Syy, Sxy;
	Sxx = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	Syy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	Sxy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);

	flDoG(Ixx, Sxx, 1.0, 2, 2);
	flDoG(Iyy, Syy, 1.0, 2, 2);
	flDoG(IxIy, Sxy, 1.0, 2, 2);

	ifsfree(Ixx, IFS_FR_ALL);
	ifsfree(Iyy, IFS_FR_ALL);
	ifsfree(IxIy, IFS_FR_ALL);


	cout << "Detecting Corner Points" << endl;
	detect = ifscreate((char *) "float", len, IFS_CR_ALL, 0);

	float detH, R, k = 0.04, threshold = 10^6;
		for (int row = 0; row < maxrow; row++)
		{
			for (int col = 0; col < maxcol; col++)
			{
				detH = ((ifsfgp(Sxx, row, col)*ifsfgp(Syy, row, col)) - (ifsfgp(Sxy, row, col)*ifsfgp(Sxy, row, col)));
				R = detH - (k * (ifsfgp(Sxx, row, col)+ifsfgp(Syy, row, col)) * (ifsfgp(Sxx, row, col)+ifsfgp(Syy, row, col)));

				if (R > threshold)
				ifsfpp(detect, row, col, R);
			}
		}


		float value = 0;
		for (int row = 0; row < maxrow; row++)
		{
			for (int col = 0; col < maxcol; col++)
			{
				value = ifsfgp(detect, row, col);
				if (value < 50)
					ifsfpp(detect, row, col, 0);
				else
				{
					ifsfpp(detect, row, col, 255);
					results.push_back(col);
					results.push_back(row);
				}
			}
		}

		ifsfree(Sxy, IFS_FR_ALL);
		ifsfree(Sxx, IFS_FR_ALL);
		ifsfree(Syy, IFS_FR_ALL);

		return results;
}


void findH(std::vector<int>& corres)
{
	int N,i,j,p;
	float det;


	N = corres.size()/4;
	//cout<<"\nno of correspondences = "<<N<<endl;

	
	float **Amatrix;
	Amatrix = matrix(1,(2*N),1,8);
	for(i=1;i<=N;i++)
	{
		p = i-1;

		Amatrix[(2*i)-1][1] = 0;
		Amatrix[(2*i)-1][2] = 0;
		Amatrix[(2*i)-1][3] = 0;
		Amatrix[(2*i)-1][4] = -corres[4*p];
		Amatrix[(2*i)-1][5] = -corres[(4*p)+1];
		Amatrix[(2*i)-1][6] = -1;
		Amatrix[(2*i)-1][7] = corres[(4*p)]*corres[(4*p)+3];
		Amatrix[(2*i)-1][8] = corres[(4*p)+1]*corres[(4*p)+3];

		Amatrix[(2*i)][1] = corres[4*p];
		Amatrix[(2*i)][2] = corres[(4*p)+1];
		Amatrix[(2*i)][3] = 1;
		Amatrix[(2*i)][4] = 0;
		Amatrix[(2*i)][5] = 0;
		Amatrix[(2*i)][6] = 0;
		Amatrix[(2*i)][7] = -corres[(4*p)]*corres[(4*p)+2];
		Amatrix[(2*i)][8] = -corres[(4*p)+1]*corres[(4*p)+2];

	}


	//double *dvector = new double[2*N];			//Remember to deallocate
	float **D;
	D = matrix(1,(2*N),1,1);


	for(i=1;i<=N;i++)
	{
		p = i-1;
		D[(2*i)-1][1]= -corres[(4*p)+3];
		D[(2*i)][1] = corres[(4*p)+2];
	}



	float **Atrans;
	Atrans = matrix(1,8,1,(2*N));
	transpose(Amatrix,(2*N),8,Atrans);


	float **At_andA;
	At_andA = matrix(1,8,1,8);
	ifsmatmult(Atrans, Amatrix, At_andA, 8,(2*N),(2*N),8);

	double **temp1;
	temp1 = dmatrix(1,8,1,8);
	for(i=1;i<=8;i++)			//Assign values of (float)At_andA to (double)temp1
	{
		for(j=1;j<=8;j++)
		{
			temp1[i][j] = double(At_andA[i][j]);
		}
	}

	double **temp2;
	temp2 = dmatrix(1,8,1,8);
	det = ifsinverse(temp1, temp2, 8);

	float **AtA_inv;
	AtA_inv = matrix(1,8,1,8);
	for(i=1;i<=8;i++)			//Assign values of (double)temp2 to (float)AtA_inv
	{
		for(j=1;j<=8;j++)
		{
			AtA_inv[i][j] = float(temp2[i][j]);
		}
	}


	float **A_triple;
	A_triple = matrix(1,8,1,(2*N));
	ifsmatmult(AtA_inv, Atrans, A_triple, 8,8,8,(2*N));

	float **h;
	h = matrix(1,8,1,1);								
	ifsmatmult(A_triple, D, h, 8,(2*N),(2*N),1);



	H = matrix(1,3,1,3);
	for(i=1;i<=3;i++)			//Convert h vector into H matrix
	{
		for(j=1;j<=3;j++)
		{
			if(i*j != 9)
			{
				H[i][j] = h[(3*(i-1)+j)][1];
			}
			else
			{
				H[3][3] = 1;
			}
		}
	}

}

void mapback()
{
	int i,j,maxcol,maxrow,det;
	int xdash,ydash,xold,yold;
	maxcol = ifsdimen(Leftimg,0);
	maxrow = ifsdimen(Leftimg,1);
	finalimg = fifs2farr(Result);

	int origin_y = 100;
	int origin_x = maxcol;

	double **temp1;
	temp1 = dmatrix(1,3,1,3);
	for(i=1;i<=3;i++)			//Assign values of (float)H to (double)temp1
	{
		for(j=1;j<=3;j++)
		{
			temp1[i][j] = double(H[i][j]);
		}
	}
	double **temp2;
	temp2 = dmatrix(1,3,1,3);


	Hinv = matrix(1,3,1,3);

	det = ifsinverse(temp1, temp2, 3);

	for(i=1;i<=3;i++)			//Assign values of (double)temp2 to (float)Hinv
	{
		for(j=1;j<=3;j++)
		{
			Hinv[i][j] = float(temp2[i][j]);
		}
	}

	cout<<"DET(Hinv) ="<<det<<endl;
	

	
	for(i=0;i<maxrow;i++)		//Write the right side image to output image
	{
		for(j=0;j<maxcol;j++)
		{
			finalimg[origin_y+i][origin_x+j] = rightimg[i][j];
		}
	}

	float **wvector;
	wvector = matrix(1,3,1,1);
	float **Xdash;
	Xdash = matrix(1,3,1,1);


	for(i=0;i<rmax_result;i++)
	{
		for(j=0;j<cmax_result;j++)
		{
			xdash = j-origin_x;
			ydash = i-origin_y;

			Xdash[1][1] = xdash;
			Xdash[2][1] = ydash;
			Xdash[3][1] = 1;

			ifsmatmult(Hinv, Xdash, wvector, 3,3,3,1);

			xold = int(wvector[1][1]/wvector[3][1]);
			yold = int(wvector[2][1]/wvector[3][1]);

			if( (xold>=0 && xold<maxcol) && (yold>=0 && yold<maxrow) )
			{
				finalimg[i][j] = leftimg[yold][xold];
			}


		}
	}
	ifspot(Result, (char *)"stitched.ifs");

}


std::vector<int> testdata()
{
	std::vector<int> corres;
	int N,i;

	int data[] = {436,100,291,92,449,48,304,44,458,65,313,60,455,209,307,199,503,63,354,62,520,63,367,62,505,96,355,93,520,96,367,93,
				506,116,355,112,521,116,368,112,507,144,356,138,523,143,369,138,510,164,357,157,525,164,370,157,510,193,357,183,
				526,193,371,184,594,60,428,66,608,60,441,66,595,93,430,96,611,93,442,96,598,114,431,115,614,114,445,115,
				601,142,433,140,616,141,445,139};


	N = sizeof(data)/4;
	cout<<"\nN= "<<N<<endl;

	for(i=0;i<N;i++)
	{
		corres.push_back(data[i]);
	}
	

	return(corres);

}

/*
std::vector<int> matchPoints(std::vector<int>& leftpoints, std::vector<int>& rightpoints, int N_left, int N_right)
{
	
	int i,j,lx,ly,rx,ry,rx_match,ry_match;
	float diff,diff_threshold,mindiff;
	std::vector<int> results;
	float left_template[5][5];
	float right_template[5][5];

	diff_threshold = 50.0;

	for(i=0;i<N_left;i++)
	{
		lx = leftpoints[2*i];
		ly = leftpoints[(2*i)+1];

		//Build a descriptor for the neighbourhood of this point
		buildTemplate(0,lx,ly, left_template);

		mindiff = 99999;
		rx_match =0;
		ry_match =0;

		for(j=0;j<N_right;j++)
		{
			rx = rightpoints[2*j];
			ry = rightpoints[(2*j)+1];

			if(abs(ry-ly)>=25 || abs(abs(rx-lx)-280)>30)				//Purely Experimental
				continue;

			buildTemplate(1,rx,ry, right_template);


			diff = diffMeasure(left_template, right_template);

			if(diff<mindiff)
			{
				mindiff = diff;
				rx_match = rx;
				ry_match = ry;
			}

		}
		if(mindiff<diff_threshold)
		{
			results.push_back(lx);
			results.push_back(ly);
			results.push_back(rx_match);
			results.push_back(ry_match);
		}

	}

	return(results);


}
*/


std::vector<int> matchPoints(std::vector<int>& leftpoints, std::vector<int>& rightpoints, int N_left, int N_right)
{
	
	int i,j,lx,ly,rx,ry,rx_match,ry_match,rx_orig,ry_orig;
	float diff,diff_threshold,mindiff;
	std::vector<int> results;
	float left_template[5][5];
	float right_template[5][5];

	diff_threshold = 50.0;

	for(i=0;i<N_left;i++)
	{
		lx = leftpoints[2*i];
		ly = leftpoints[(2*i)+1];

		//Build a descriptor for the neighbourhood of this point
		buildTemplate(0,lx,ly, left_template);

		mindiff = 99999;
		rx_match =0;
		ry_match =0;

		for(j=0;j<N_right;j++)
		{
			rx_orig = rightpoints[2*j];
			ry_orig = rightpoints[(2*j)+1];

			if(abs(ry_orig-ly)>=25 || abs(abs(rx_orig-lx)-280)>30)				//Purely Experimental
				continue;

			for(rx=rx_orig-10;rx<=rx_orig+10;rx++)
			{
				for(ry=ry_orig-10;ry<=ry_orig+10;ry++)
				{
					buildTemplate(1,rx,ry, right_template);

					diff = diffMeasure(left_template, right_template);

					if(diff<mindiff)
					{
						mindiff = diff;
						rx_match = rx;
						ry_match = ry;
					}
				}
			}

		}

		if(mindiff<diff_threshold)
		{
			results.push_back(lx);
			results.push_back(ly);
			results.push_back(rx_match);
			results.push_back(ry_match);
		}

	}

	return(results);


}


void imageGradients()
{
	int temp;

	temp = flDoG(Leftimg,ifs_Ldx,1,1,0);
	temp = flDoG(Leftimg,ifs_Ldy,1,1,1);

	temp = flDoG(Rightimg,ifs_Rdx,1,1,0);
	temp = flDoG(Rightimg,ifs_Rdy,1,1,1);

	Ldx = fifs2farr(ifs_Ldx);
	Ldy = fifs2farr(ifs_Ldy);

	Rdx = fifs2farr(ifs_Rdx);
	Rdy = fifs2farr(ifs_Rdy);
	return;

}


void RANSAC(std::vector<int>& corres)
{
	int i,j,r,N,duplicate,index,test,ctr,iterations,testpts,maxctr;
	int randompts[12];
	std::vector<int> new_corres;

	testpts = 12;

	float **leftpt, **rightpt, **temp, **final;
	leftpt = matrix(1,3,1,1);
	rightpt = matrix(1,3,1,1);
	temp = matrix(1,3,1,1);
	final = matrix(1,1,1,1);

	N = corres.size()/4;
	cout<<"\nno of correspondences = "<<N<<endl;

	iterations=0;
	do
	{
		maxctr=0;
		for(i=0;i<testpts;i++)		//Choosing 10 random points
		{	
			duplicate=0;
			r = rand()%N;
			for(j=0;j<i;j++)
			{
				if(randompts[j]==r)
					duplicate=1;
			}

			if(!duplicate)
			{
				randompts[i]=r;
				;
			}
			else
			{
				i-=1;
			}

		}

		
		for(i=0;i<testpts;i++)
		{
			index = randompts[i];

			new_corres.push_back(corres[(4*index)]);
			new_corres.push_back(corres[(4*index)+1]);
			new_corres.push_back(corres[(4*index)+2]);
			new_corres.push_back(corres[(4*index)+3]);
		}

		findH(new_corres);

		ctr=0;
		for(i=0;i<N;i++)
		{
			leftpt[1][1] = corres[(4*i)];
			leftpt[2][1] = corres[(4*i)+1];
			leftpt[3][1] = 1.0;

			rightpt[1][1] = corres[(4*i)+2];
			rightpt[2][1] = corres[(4*i)+3];
			rightpt[3][1] = 1.0;

			test = ifsmatmult(H,leftpt,temp,3,3,3,1);
			//test = ifsmatmult(leftpt, temp,final,1,3,3,1);

			
			if(temp[3][1]==0)
				continue;

			temp[1][1] = temp[1][1]/temp[3][1];
			temp[2][1] = temp[2][1]/temp[3][1];
			temp[3][1] = temp[3][1]/temp[3][1];

			test=0;
			for(j=1;j<=2;j++)
			{
				test+= sqrt((temp[j][1] - rightpt[j][1])*(temp[j][1] - rightpt[j][1]) );
			}
			
			//cout<<final[1][1]<<"\t";

			if(test<7)
			{
				ctr++;
			}

		}
		if(ctr>maxctr)
			maxctr=ctr;
		
		iterations+=1;

	}while( ctr<= (0.5*N) && iterations<=100);

	cout<<"Ransac iterations = "<<iterations<<endl;
	cout<<"Inliers = "<<ctr<<endl;
	cout<<"Max = "<<maxctr<<endl;

	return;
}


int main()
{

	int maxcol,maxrow, N_left, N_right,i,N;
	std::vector<int> corres,test;
	std::vector<int> leftpoints, rightpoints;


	//Import both images
	Leftimg = ifspin((char *)"leftpic.ifs");
	leftimg = fifs2farr(Leftimg);
	Rightimg = ifspin((char *)"rightpic.ifs");
	rightimg = fifs2farr(Rightimg);

	//Get Sizes
	
	maxcol = ifsdimen(Leftimg,0);			
	maxrow = ifsdimen(Leftimg,1);

	cmax_result = 2*maxcol;
	rmax_result = maxrow+100;
	
	int len[3] = {2,cmax_result,rmax_result};
	Result = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	ifs_Ldx = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	ifs_Ldy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	ifs_Rdx = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	ifs_Rdy = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
	//finalimg = fifs2farr(Result);

	imageGradients();

	
	leftpoints = detectCorner(Leftimg);	
	rightpoints = detectCorner(Rightimg);
	//Assume they're both vectors, with index 2n representing x cood and 2n+1 representing y cood of nth interest point
	
	N_left  = leftpoints.size()/2;
	N_right = rightpoints.size()/2;
	

	//Write a function to find correspondences and return a vector with only corresponding points in it,
	// in the format (x1,y1,x2,y2)??
	
	cout<<"N_left = "<<N_left<<endl;
	cout<<"N_right = "<<N_right<<endl;

	corres = matchPoints(leftpoints, rightpoints,N_left, N_right);
	
	N = corres.size()/4;
	cout<<"N = "<<N<<endl;

	for(i=0;i<N;i++)
	{
		cout<<corres[(4*i)]<<" "<<corres[(4*i)+1]<<"\t"<<corres[(4*i)+2]<<" "<<corres[(4*i)+3]<<endl;
	}
	
	//corres = testdata();


	//findH(corres);

	RANSAC(corres);
	mapback();


}


