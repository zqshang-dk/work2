#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXSTRPATH  256      /* max length of file path */
#define MAXOBS      50		 /* max number of observations */

/* decode angle observations in degree, minutes and seconds --------------------
*  args  ：  char*   p      I   strings that contains angle obs in deg, min and sec
*            int     n      I   number of decoded value
*            double* v      O   the array where decoded values are stored
*  return    int     i      :    the number of decoded value
*-----------------------------------------------------------------------------*/
static int decodef(char *p, int n, double *v)
{
	int i;
	const char *delim = " ";

	/* initialize the arrays v */
	for (i = 0; i < n; i++) v[i] = 0;

	/* split strings and converted to 'double' */
	for (i = 0, p = strtok(p, delim); p&&i < n; p = strtok(NULL, delim)) {
		v[i++] = atof(p);
	}

	return i;
}

/* point positioning functions -------------------------------------------------
*  args  ：  double* obs    I   obs data{X, Y, Z, pseudorange};
*            int     n      I   number of obs
*            double* x      O   unknown parameters
*  return    :    0:failed; 1:sucess
*-----------------------------------------------------------------------------*/
static int pntpos(const double obs[][4], const int n, double* x) 
{
	int stat = 1;

	return stat;
}

/* main functions ------------------------------------------------------------*/
int main(int argc, char *argv[]) {
	FILE *fp,*fop;
	int i = 0, j,epochnum,gpsw,obsnum=0;
	double obs[MAXOBS][4] = { 0 }, gpss, x[4] = { 0 };
	char filename[MAXSTRPATH] = "E:/STUDY/Sophomore1/最优估计/第二次上机实习/work2/CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt";  /* 数据文件路径(绝对路径) */
	char outfile[MAXSTRPATH] = "outpos.txt"; /* 数据输出文件(与源码同路径)*/
	char buff[MAXSTRPATH],satid[4];

	if (!(fp = fopen(filename, "r")) || !(fop = fopen(outfile, "w+"))) {
		printf("error file directory, please check again\n");
		return 0;
	}

	/* 逐历元读取数据 */
	while (fgets(buff, sizeof(buff), fp)) {
		j = 0;
		if (strlen(buff) < 20) continue;

		if (strchr(buff, '#')) {
			sscanf(buff + 2, "%d %d %lf %d", &epochnum, &gpsw, &gpss, &obsnum);
		}
		if (obsnum < 4) continue;

		for (i = 0; i < obsnum; i++) {
			fgets(buff, sizeof(buff), fp);

			/* 卫星prn */
			strncpy(satid, buff, 3);
			satid[3] = '\0';

			/* 逐行解析卫星坐标和伪距观测值 */
			decodef(buff + 3, 4, obs[j++]);
		}

		/* 定位解算函数 */
		pntpos(obs, obsnum, x);

		fprintf(fop, "%2d %10.3lf %14.3lf %14.3lf %14.3lf\n", gpsw, gpss, x[0], x[1], x[2]);
	}

	fclose(fp); fclose(fop);
	return 0;
}