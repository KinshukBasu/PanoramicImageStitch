#include <iostream>
#include "math.h"
#include "ifs.h"
#include "flip.h"

using namespace std;

IFSIMG detectCorner(IFSIMG img, int len[]);

int main() {
IFSIMG leftPic, rightPic;

cout << "Reading images from disk" << endl;
leftPic = ifspin((char *) "leftpic.ifs");
rightPic = ifspin((char *) "rightpic.ifs");

cout << "Checking the size of the image" <<endl;
int *len = ifssiz(leftPic);


IFSIMG cornerLeft, cornerRight;

cornerLeft = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
cornerLeft = detectCorner(leftPic, len);
cout << "Writing the contents of cornerLeft to the disk" << endl;
ifspot(cornerLeft,(char *) "CornerLeft.ifs");

// overlapping points on the left image
for (int row = 0; row < len[2]; row++) {
	for (int col = 0; col < len[1]; col++) {
		if (ifsfgp(cornerLeft, row, col) == 255)
			ifsfpp(leftPic, row, col, 255);
	}
}

cout << "Writing the contents of leftPic to the disk" << endl;
ifspot(leftPic,(char *) "overlapLeft.ifs");


cornerRight = ifscreate((char *) "float", len, IFS_CR_ALL, 0);
cornerRight = detectCorner(rightPic, len);
cout << "Writing the contents of cornerRight to the disk" << endl;
ifspot(cornerRight,(char *) "CornerRight.ifs");

// overlapping points on the right image
for (int row = 0; row < len[2]; row++) {
	for (int col = 0; col < len[1]; col++) {
		if (ifsfgp(cornerRight, row, col) == 255)
			ifsfpp(rightPic, row, col, 255);
	}
}

cout << "Writing the contents of rightPic to the disk" << endl;
ifspot(rightPic,(char *) "overlapRight.ifs");

// Freeing len
free(len);

// Freeing Images
ifsfree(leftPic, IFS_FR_ALL);
ifsfree(rightPic, IFS_FR_ALL);
ifsfree(cornerLeft, IFS_FR_ALL);
ifsfree(cornerRight, IFS_FR_ALL);

return 0;
}


IFSIMG detectCorner (IFSIMG img, int len[]) {

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

	/*
	cout << "Writing the contents of Sxx to the disk" << endl;
	ifspot(Sxx,(char *) "Sxx.ifs");
	cout << "Writing the contents of Syy to the disk" << endl;
	ifspot(Syy,(char *) "Syy.ifs");
	cout << "Writing the contents of Sxy to the disk" << endl;
	ifspot(Sxy,(char *) "Sxy.ifs");
	*/

	ifsfree(Ixx, IFS_FR_ALL);
	ifsfree(Iyy, IFS_FR_ALL);
	ifsfree(IxIy, IFS_FR_ALL);


	cout << "Detecting Corner Points" << endl;
	detect = ifscreate((char *) "float", len, IFS_CR_ALL, 0);

	float detH, R, k = 0.04, threshold = 10^6;
		for (int row = 0; row < maxrow; row++) {
			for (int col = 0; col < maxcol; col++) {
				detH = ((ifsfgp(Sxx, row, col)*ifsfgp(Syy, row, col)) - (ifsfgp(Sxy, row, col)*ifsfgp(Sxy, row, col)));
				R = detH - (k * (ifsfgp(Sxx, row, col)+ifsfgp(Syy, row, col)) * (ifsfgp(Sxx, row, col)+ifsfgp(Syy, row, col)));

				if (R > threshold)
				ifsfpp(detect, row, col, R);
			}
		}

/*
		for(int row = 0; row < maxrow; row++) {
			for (int col = 0; col < maxcol; col++) {
				if((ifsfgp(detect,row,col-1) == 0 && ifsfgp(detect,row,col+1) == 0) || (ifsfgp(detect,row-1,col) == 0 && (ifsfgp(detect,row+1,col) == 0 )) || (ifsfgp(detect,row+1,col-1) == 0 && ifsfgp(detect,row-1,col+1) == 0) || (ifsfgp(detect,row-1,col-1) == 0 && (ifsfgp(detect,row+1,col+1) == 0 )))
				ifsfpp(detect,row,col,0);
			}
		}
*/

		float value = 0;
		for (int row = 0; row < maxrow; row++) {
			for (int col = 0; col < maxcol; col++) {
					value = ifsfgp(detect, row, col);
					if (value < 50)
						ifsfpp(detect, row, col, 0);
					else
						ifsfpp(detect, row, col, 255);
		}
	}

		ifsfree(Sxy, IFS_FR_ALL);
		ifsfree(Sxx, IFS_FR_ALL);
		ifsfree(Syy, IFS_FR_ALL);

		return detect;
}
