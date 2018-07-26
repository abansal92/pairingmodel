#include <stdio.h>  	
#include <stdlib.h>		
#include <math.h>


/*
                   _ooOoo_
                  o8888888o
                  88" . "88
                  (| -_- |)
                  O\  =  /O
               ____/`---'\____
             .'  \\|     |//  `.
            /  \\|||  :  |||//  \
           /  _||||| -:- |||||-  \
           |   | \\\  -  /// |   |
           | \_|  ''\---/''  |   |
           \  .-\__  `-`  ___/-. /
         ___`. .'  /--.--\  `. . __
      ."" '<  `.___\_<|>_/___.'  >'"".
     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
     \  \ `-.   \_ __\ /__ _/   .-` /  /
======`-.____`-.___\_____/___.-`____.-'======
                   `=---='
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           ·ð×æ±£ÓÓ       ÓÀÎÞBUG
*/


#define MaxDim 500
#define MaxDegenerate 700
#define Maxjm  20
#define Maxnumber 30

#define sgn(x) ( (x) % 2 ? -1 : 1 )

int particle_number, j_number;
int j[Maxjm], l_number[Maxjm], s_number[Maxjm], N_number[Maxjm];
double M_uplimit;

int temp_n[MaxDim][Maxnumber];
double temp_Jz[MaxDim][Maxnumber];
int temp_count;

//int ni_m[500][50];

//int part_N[Maxjm][500][Maxnumber];
//double part_Jz[Maxjm][500][Maxnumber];
//int part_particles_num;
int MaxtixDim;


int degenerate;
int rangeforj[Maxjm];   //j
double totalM[MaxDim];  //m
double parameterG;
double energyspace;

double eps = 10e-9;

static double Tmatrix[50][50][50][50];

void ParticlesDistribution();
void ProduceMbasis();
double Overlap(int *ket, int *bar);
double OverlapBit(int *ket, int *bar, int nump);
void OnebodyOperator();
double TwobodyOperator(int *ket, int *bar, int nump);
double FockMatrix(int p, int q);
double HbarN(int ap, int bp, int ip, int jp);
void InitTmatix();

void iteration();


void main()
{
	int ket[10], bar[10];
	int ap, bp, ip, jp;
	parameterG = -0.25;
	energyspace = 1.;
	ParticlesDistribution();
	InitTmatix();

	//printf("Tmaxtix  %f \n", Tmatrix[4][5][0][1] + HbarN(4, 5, 0, 1) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)) );
	//printf("Tmaxtix  %f \n", Tmatrix[4][5][1][0] + HbarN(4, 5, 1, 0) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)));


	iteration();

	//ProduceMbasis();
	//printf("produce M-scheme basis ok!!!\n");
	//OnebodyOperator();
	//printf("produce One-body Operator elements ok!!!\n");
	// TEST
	//ket
	//int phase;
	//phase = 0;//5621
/*	ket[0] = 0;
	ket[1] = 0;
	ket[2] = 0;
	ket[3] = 0;
	ket[4] = 1;
	ket[5] = 1;
	ket[6] = 0;
	ket[7] = 0;
	//bar
	bar[0] = 1;
	bar[1] = 1;
	bar[2] = 0;
	bar[3] = 0;
	bar[4] = 0;
	bar[5] = 0;
	bar[6] = 0;
	bar[7] = 0;
	*/
	/*	for (ap = particle_number; ap < degenerate; ap++)
		{
			for (bp = particle_number; bp < degenerate; bp++)
			{
				for (ip = 0; ip < particle_number; ip++)
				{
					for (jp = 0; jp < particle_number; jp++)
					{
						if ( HbarN(ap,bp,ip,jp) != 0  )
						printf("Hn %d %d %d %d  %lf \n",ap+1,bp+1,ip+1,jp+1 ,HbarN(ap, bp, ip, jp));
					}
				}
			}
		}*/

		/*for (ap = particle_number; ap < degenerate; ap++)
		{
			for (bp = particle_number; bp < degenerate; bp++)
			{
				for (ip = 0; ip < particle_number; ip++)
				{
					for (jp = 0; jp < particle_number; jp++)
					{
						if ( (Tmatrix[ap][bp][ip][jp] + Tmatrix[ap][bp][jp][ip])!= 0)
							printf("T  %d %d %d %d  %lf  %lf  %lf \n", ap + 1, bp + 1, ip + 1, jp + 1, (Tmatrix[ap][bp][ip][jp]+ Tmatrix[ap][bp][jp][ip]), Tmatrix[ap][bp][ip][jp], Tmatrix[ap][bp][jp][ip]);
					}
				}
			}
		}*/


	//phase++;
	//printf("Fock  %lf \n", FockMatrix(4,4));
	//printf("Fock  %lf \n", FockMatrix(5, 5));
	//printf("Fock  %lf \n", FockMatrix(0, 0));
	//printf("Fock  %lf \n", FockMatrix(1, 1));
	//printf("V  %lf \n", TwobodyOperator(ket, bar,2));
	//printf("Test overlap  %lf \n", OverlapBit(ket, ket, 4));
	//printf("Tmaxtix  %f \n", Tmatrix[4][5][0][1]);
	//printf("Tmaxtix  %f \n", Tmatrix[4][5][1][0]);

	//printf("Tmaxtix  %f \n", Tmatrix[4][5][0][1] + HbarN(4, 5, 0, 1) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)) );
	//printf("Tmaxtix  %f \n", Tmatrix[4][5][1][0] + HbarN(4, 5, 1, 0) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)));

	//printf("Tmaxtix  %f \n", Tmatrix[6][4][1][0]);	
	//printf("Tmaxtix  %f \n", Tmatrix[4][5][0][1]);
	//printf("Hn   %lf \n", HbarN(4, 5, 2, 3));
	//printf("Hn   %lf \n", HbarN(4, 5, 3, 2));
	//printf("Hn   %lf \n", HbarN(5, 5, 1, 1));
	//printf("Hn   %lf \n", HbarN(4, 5, 1, 0));
}

  
void ParticlesDistribution()
{
	FILE *fp, *wfp, *wfp2;
	char cmd[300];
	int i,i1;
	int N, NJ_sum, num;
	int l, N_J[MaxDim], l1, tempnum;
	int particleDistribution[MaxDim];
	char *feedbackp;
	double tempM;


	num = 0;
	fp = fopen("input.dat", "r");
	wfp = fopen("basis/ParticlesDistribution.dat", "w");
	wfp2 = fopen("basis/ParticlesDistribution2.dat", "w");
	feedbackp = fgets(cmd, 300, fp);
	if (sscanf(cmd, "n=%d  jm=%d", &particle_number, &j_number) != 2)
	{
		printf("error reading input files\n");
		exit(4);
	}
	if (j_number > Maxjm)
	{
		printf("too many shells\n");
		exit(4);
	}
	for (i = 0; i < j_number; i++)
	{
		if (fscanf(fp, "%d %d %d %d", &N_number[i], &j[i], &l_number[i], &s_number[i]) != 4)
		{
			printf("input file err!!!\n");
		}
		//printf("%d ", j[i]);
	}
	//printf("%d  %d\n", j_number, particle_number);

	degenerate = 0;
	rangeforj[0] = 0;
	for (i = 0; i < j_number; i++)                         /* folding */
	{
		for (i1 = 0; i1 < j[i] + 1; i1++)
		{
			rangeforj[degenerate+i1] = i;
			totalM[degenerate + i1] = -1.*j[i] * 1.0 / 2.0 + i1;
			//printf("%f  ", totalM[degenerate + i1]);
		}

		degenerate += j[i]+1;
	}
	//printf("\n");

	NJ_sum = 1;
	for (i = 0; i < degenerate; i++)                         /* folding */
	{

		N_J[i] = 2;
		//printf("  %d  ,",N_J[i]);
		NJ_sum *= N_J[i];
		//printf("  %d  ,",NJ_sum);
	}

	for (l = 0; l <= NJ_sum; l++)
	{	                                   /* unfolding it */
		l1 = 1;
		for (i = 0; i < degenerate; i++)
		{
			particleDistribution[i] = (l / l1) % N_J[i];
			l1 *= N_J[i];
			//printf(" %d ", particleDistribution[i]);
		}
		//printf("\n");
		tempnum = 0;
		for (i = 0; i < degenerate; i++)
		{
			tempnum += particleDistribution[i];
		}
		//printf("%d \n", tempnum);
		if (tempnum != particle_number)
		{
			//printf(" out ! \n");
			continue;
		}
		/*********** write files ***********/
		fprintf(wfp, "%d", num);
		fprintf(wfp2, "%d", num);
		for (i = 0; i < degenerate; i++)
		{
			fprintf(wfp, "	%d", particleDistribution[i]);
			if (particleDistribution[i] == 1)
			{
				fprintf(wfp2, "	%d", i+1);
			}
		}
		num++;
		tempM = 0;
		for (i = 0; i < degenerate; i++)
		{
			//rangeforj[i]
			if (particleDistribution[i] == 1)
			{
				tempM += totalM[i];
			}
		}
		if (tempM > M_uplimit)
		{
			M_uplimit = tempM;
		}
		fprintf(wfp, "	%10.1f", tempM);
		fprintf(wfp, "\n");
		fprintf(wfp2, "	%10.1f", tempM);
		fprintf(wfp2, "\n");
	}

	fclose(fp);
	fclose(wfp);
	fclose(wfp2);
	MaxtixDim = num;

}

void ProduceMbasis()
{
	FILE *fp, *mfp;
	int i, j, num, num_PD,tempcount;
	double M_max_alpha[Maxjm], M_we_cal;
	static int ni_m[500][30];
	int count;
	char f1[200];


	fp = fopen("basis/ParticlesDistribution2.dat", "r");
	i = j_number;
	num = 1;
	//printf("ok!!\n");
	while (i == particle_number)
	{
		//printf("index:%d  ", num);
		i = 0;
		if (fscanf(fp, "%d", &tempcount) != 1)
		{
			break;
		}

		for (i = 0; i < particle_number; i++)
		{
			if (fscanf(fp, "%d", &ni_m[num][i]) != 1)
			{
				break;
			}
			//printf(" %d", ni_m[num][i]);
		}
		//printf("\n");
		if (fscanf(fp, "%lf", &M_max_alpha[num]) != 1)
		{
			break;
		}
		num++;
		if (num > MaxDim)
		{
			printf("too many distributions, terminated!!!\n");
			fclose(fp);
			return;
		}
	}
	fclose(fp);
	//printf("****************\n");


	count = 0;
	
	for (M_we_cal = M_uplimit; M_we_cal >= -M_uplimit; M_we_cal--)
	{

		sprintf(f1, "basis/M%.1f.dat", M_we_cal);
		mfp = fopen(f1, "w");
		//fprintf(mfp, "	M	num\n");
		for (num_PD = 1; num_PD < num; num_PD++) //alpha
		{
			//temp_count = 0;
			//printf("PD=%d", num_PD);
			if (M_max_alpha[num_PD] != M_we_cal)
			{
				continue;
			}
			
			for (j = 0; j < particle_number; j++)
			{
				fprintf(mfp, "%d	", ni_m[num_PD][j]);
				printf(" %d", ni_m[num_PD][j]);
			}
			printf("   %d\n",num_PD);
			fprintf(mfp, "%.1lf\n", M_we_cal);

		}
		fclose(mfp);
	}
	return;
}

double Overlap(int *ket, int *bar)
{
	int i, j;
	int tempket[MaxDim], tempbar[MaxDim];
	for ( i = 0; i <particle_number ; i++)
	{
		tempket[ ket[i] ] += 1;
		tempbar[ bar[i] ] += 1;
		if (tempket[ket[i]] > 1  || tempbar[bar[i]] > 1)
		{
			return 0.;
		}
	}
	for ( i = 0; i < degenerate; i++)
	{
		if ( ket[i] != bar[i] )
		{
			return 0.;
		}
	}
	return 1.;
}

double OverlapBit(int *ket, int *bar, int nump)
{
	int i, j;
	int tempket[MaxDim], tempbar[MaxDim];
	int templ, tempr;
	/*for (i = 0; i <particle_number; i++)
	{
		tempket[ket[i]] += 1;
		tempbar[bar[i]] += 1;
		if (tempket[ket[i]] > 1 || tempbar[bar[i]] > 1)
		{
			return;
		}
	}*/
	templ = 0;
	tempr = 0;
	for (i = 0; i < degenerate; i++)
	{
		templ += ket[i];
		tempr += bar[i];
		if (ket[i] != bar[i] || ket[i] > 1 || bar[i] >1 )
		{
			return 0.;
		}
	}
	if (templ != nump || tempr != nump)
	{
		//printf("number of particles in overlap err%d %d!!\n", templ, tempr);
		return 0.;
	}
	return 1.;
}

void OnebodyOperator()
{
	int i, i1;
	FILE *fp;
	double M_max_alpha[MaxDegenerate];
	int ni[MaxDim][200];
	int num, num_PD;
	double M_cal;
	char cmd[100];
	double Matrixele[MaxDegenerate][MaxDegenerate];
	int temp_i, temp_j;

	/*******  defined one-body opertor  *******/
	int j1, j2, m1, m2;
	double OBO[MaxDegenerate][MaxDegenerate];
	if (degenerate > MaxDegenerate)
	{
		printf("too great dim!!\n");
		return;
	}
	//printf("%d\n",degenerate);
	for (m1 = 0; m1 < degenerate; m1++)
	{
		for (m2 = 0; m2 < degenerate; m2++)
		{
			OBO[m1][m2] = 0;
			if (m1 != m2)
			{
				continue;
			}
			OBO[m1][m2] = 1;
		}
	}
	/******************************************/
	//read basis
	num = 0;
	sprintf(cmd, "basis/ParticlesDistribution2.dat");
	fp = fopen(cmd, "r");
	i = j_number;
	//printf("ok!!\n");
	while (i == particle_number)
	{
		i = 0;
		if (fscanf(fp, "%d", &num_PD) != 1)
		{
		break;
		}
		for (i = 0; i < particle_number; i++)
		{
			if (fscanf(fp, "%d", &ni[num][i]) != 1)
			{
				break;
			}
			//printf(" %d", ni[num][i]);
		}
		//printf("\n");
		if (fscanf(fp, "%lf", &M_max_alpha[num]) != 1)
		{
			break;
		}
		//printf("   M=%lf\n", M_max_alpha[num]);
		num++;
		if (num > MaxDim)
		{
			printf("too many distributions, terminated!!!\n");
			fclose(fp);
			return;
		}
	}
	fclose(fp);

	//printf("%d\n", num);

	/******************************************/
	//cal matrix elements
	for (M_cal = M_uplimit; M_cal >= -1. * M_uplimit; M_cal--)
	{
		//printf("%f  %d\n", M_cal,num);

		sprintf(cmd,"elements/OneBody_%.1f.dat", M_cal);
		fp = fopen(cmd, "w");
		for ( m1 = 0; m1 < num; m1++)
		{
			temp_i = 0;
			temp_j = 1;
			for (m2 = 0; m2 < num; m2++)
			{
				if (M_max_alpha[m1] != M_max_alpha[m2] || M_max_alpha[m1] != M_cal )
				{
					continue;
				}
				temp_i++;
				Matrixele[m1][m2] = OBO[m1][m2]*Overlap(ni[m1],ni[m2]);
				fprintf(fp, "%.5lf	", Matrixele[m1][m2]);
			}//m1
			if ( temp_i != 0)
			{
				fprintf(fp, "\n");
				//temp_j++;
			}

		}
		fclose(fp);
	}//M_cal
}

double TwobodyOperator(int *ket, int *bar, int nump)
{
	int i, i1,k,kp,he1;
	FILE *fp;
	int ni[MaxDim][200];
	int num, num_PD;
	double M_cal;
	char cmd[100];
	int temp_i, temp_j;
	int tempket[100], tempbar[100];
	int position[MaxDegenerate], positioncut;
	double val;
	int phase;
	/*******  defined two-body opertor  *******/
	// H = -G \sum (a*k+) (a*k-) (ak-) (ak+)  
	int j1, j2, m1, m2;
	if (degenerate > MaxDegenerate)
	{
		printf("too great dim!!\n");
		return 0.;
	}
	//printf("%d\n",degenerate);
	positioncut = 0;
	for (i = 0; i < degenerate; i++)
	{
		tempbar[i] = bar[i];
		//printf("%d ", bar[i]);
		if ( bar[i] == 1 )
		{
			position[positioncut] = i;
			positioncut++;
		}
		else if (bar[i] > 1)
		{
			printf("basis err!!!\n");
		}
	}

	/******************************************/
	//cal matrix elements
	val = 0;
	kp = degenerate / 2;
	//printf("dim %d   %d\n",degenerate,kp);
	for ( k = 0; k < kp; k++)
	{
		for (he1 = 0; he1 < kp; he1++)
		{
			//printf(" %d %d \n", k,he1);
			for (i = 0; i < degenerate; i++)
			{
				tempbar[i] = bar[i];
				//printf(" %d ", tempbar[i]);
			}
			tempbar[2 * k]--;
			tempbar[2 * k + 1]--;
			if (tempbar[2 * k] < 0 || tempbar[2 * k + 1] < 0)
			{
				//printf("%d  %d  %d \n", tempbar[2 * k], tempbar[2 * k + 1], k);
				continue;
			}
			tempbar[2 * he1 + 1]++;
			tempbar[2 * he1]++;
			if (tempbar[2 * kp] >= 2  || tempbar[2 * kp + 1] >= 2 )
			{
				continue;
			}
			//printf("       %f   %f  \n", val,  OverlapBit(ket, tempbar, nump));
			val +=  OverlapBit(ket, tempbar, nump);
		}
	}
	val = -1. * parameterG * val;
	return val;
}

double FockMatrix(int p, int q)
{
	double val;
	int i,j;
	int ket[20], bar[20];
	int phase;
	double singleparticleenergy[100];

	for ( i = 0; i < degenerate; i++)
	{
		singleparticleenergy[i] = energyspace * N_number[i/2];
		//printf("  %d  %f \n", i, singleparticleenergy[i]);
	}

	val = 0.;
    if ( p == q )
	{
		val += singleparticleenergy[p];
	}
	//printf("%d   %f   %f\n", p,val,singleparticleenergy[p]);
	for (i = 0; i < particle_number; i++)
	{
		if ( i == p || i ==q )
		{
			continue;
		}
		for (j = 0; j < degenerate; j++)
		{
			ket[j] = 0;
			bar[j] = 0;
		}
		phase = 0;
		ket[p] = 1;
		ket[i] = 1;
		if (i> p)
		{
			phase++;
		}
		bar[q] = 1;
		bar[i] = 1;
		if (i> q)
		{
			phase++;
		}
		val += sgn(phase) * TwobodyOperator(ket,bar,2);
		//printf("TTT  %f  \n", TwobodyOperator(ket, bar, 2));
		/*if ( TwobodyOperator(ket, bar, 2) != 0. )
		{
			printf("%d  \n",i);
		}*/
	}
	
	return val;

}

double HbarN(int ap, int bp, int ip, int jp)
{
	double val, tempval;
	int tempk[50], tempb[50];
	int he1,he2,he3,he4,he5;
	int phase;

	val = 0;

	// 1st term
	for ( he1 = 0; he1 < degenerate; he1++)
	{
		tempb[he1] = 0;
		tempk[he1] = 0;
	}
	phase = 0;
	tempk[ap] = 1;
	tempk[bp] = 1;
	if (bp > ap)
	{
		phase++;
	}
	tempb[ip] = 1;
	tempb[jp] = 1;
	if (jp > ip)
	{
		phase++;
	}
	val += sgn(phase) * TwobodyOperator(tempk,tempb,2);
	//printf("test  %lf\n", sgn(phase) * TwobodyOperator(tempk, tempb, 2));

	//2nd term
	//val = 0;
	for (he1 = particle_number; he1 < degenerate; he1++)  // sum c
	{
		val += FockMatrix(bp, he1)*Tmatrix[ap][he1][ip][jp] - FockMatrix(ap, he1)*Tmatrix[bp][he1][ip][jp];
		//printf("test  %lf   %f  %f\n", val, FockMatrix(bp, he1),Tmatrix[ap][he1][ip][jp]);
	}
	//printf("test  %lf\n", val);

	//3th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		val += FockMatrix(he1, ip)*Tmatrix[ap][bp][jp][he1] - FockMatrix(he1,jp)*Tmatrix[ap][bp][ip][he1];
	}
	//printf("test  %lf\n", val);

	//4th term
	//val = 0;
	for (he1 = particle_number; he1 < degenerate; he1++)  // sum c
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum d
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[ap] = 1;
			tempk[bp] = 1;
			if (bp > ap)
			{
				phase++;
			}
			tempb[he1] = 1;
			tempb[he2] = 1;
			if (he2 > he1)
			{
				phase++;
			}
			val += 0.5 * sgn(phase) * TwobodyOperator(tempk,tempb,2) * Tmatrix[he1][he2][ip][jp];
		}
	}
	//printf("test  %lf\n", val);


	//5th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = 0; he2 < particle_number; he2++)  // sum l
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[he1] = 1;
			tempk[he2] = 1;
			if (he2 > he1)
			{
				phase++;
			}
			tempb[ip] = 1;
			tempb[jp] = 1;
			if (jp > ip)
			{
				phase++;
			}
			val += 0.5 * sgn(phase) * TwobodyOperator(tempk,tempb,2) * Tmatrix[ap][bp][he1][he2];
		}
	}
	//printf("test  %lf\n", val);

	//6th term
	tempval = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[he1] = 1;
			tempk[bp] = 1;
			if (bp > he1)
			{
				phase++;
			}
			tempb[he2] = 1;
			tempb[jp] = 1;
			if (jp > he2)
			{
				phase++;
			}
			tempval += sgn(phase) * TwobodyOperator(tempk,tempb,2) * Tmatrix[ap][he2][ip][he1];
		}
	}
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[he1] = 1;
			tempk[ap] = 1;
			if (ap > he1)
			{
				phase++;
			}
			tempb[he2] = 1;
			tempb[jp] = 1;
			if (jp > he2)
			{
				phase++;
			}
			tempval += -1.0 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * Tmatrix[bp][he2][ip][he1];
		}
	}
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[he1] = 1;
			tempk[bp] = 1;
			if (bp > he1)
			{
				phase++;
			}
			tempb[he2] = 1;
			tempb[ip] = 1;
			if (ip > he2)
			{
				phase++;
			}
			tempval += -1.0 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * Tmatrix[ap][he2][jp][he1];
		}
	}
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he3 = 0; he3 < degenerate; he3++)
			{
				tempb[he3] = 0;
				tempk[he3] = 0;
			}
			phase = 0;
			tempk[he1] = 1;
			tempk[bp] = 1;
			if (bp > he1)
			{
				phase++;
			}
			tempb[he2] = 1;
			tempb[jp] = 1;
			if (jp > he2)
			{
				phase++;
			}
			tempval += 0.5 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * Tmatrix[ap][he2][ip][he1];
		}
	}
	//val += tempval;

	//7th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he4 = 0; he4 < particle_number; he4++)  // sum l
			{
				for (he5 = particle_number; he5 < degenerate; he5++)  // sum d
				{
					for (he3 = 0; he3 < degenerate; he3++)
					{
						tempb[he3] = 0;
						tempk[he3] = 0;
					}
					phase = 0;
					tempk[he1] = 1;
					tempk[he4] = 1;
					if (he4 > he1)
					{
						phase++;
					}
					tempb[he2] = 1;
					tempb[he5] = 1;
					if (he5 > he2)
					{
						phase++;
					}
					tempval = Tmatrix[ap][he2][ip][he1] * Tmatrix[he5][bp][he4][jp] - Tmatrix[ap][he2][jp][he1] * Tmatrix[he5][bp][he4][ip];
					tempval += Tmatrix[bp][he2][jp][he1] * Tmatrix[he5][ap][he4][ip] - Tmatrix[bp][he2][ip][he1] * Tmatrix[he5][ap][he4][jp];
					val += 0.5 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * tempval;
				}
			}
		}
	}
	//printf("test  %lf\n", val);

	//8th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he4 = 0; he4 < particle_number; he4++)  // sum l
			{
				for (he5 = particle_number; he5 < degenerate; he5++)  // sum d
				{
					for (he3 = 0; he3 < degenerate; he3++)
					{
						tempb[he3] = 0;
						tempk[he3] = 0;
					}
					phase = 0;
					tempk[he1] = 1;
					tempk[he4] = 1;
					if (he4 > he1)
					{
						phase++;
					}
					tempb[he2] = 1;
					tempb[he5] = 1;
					if (he5 > he2)
					{
						phase++;
					}
					tempval = Tmatrix[he2][he5][ip][he1] * Tmatrix[ap][bp][he4][jp] - Tmatrix[he2][he5][jp][he1] * Tmatrix[ap][bp][he4][ip];
					val += 0.5 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * tempval;
				}
			}
		}
	}
	//printf("test  %lf\n", val);

	//9th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he4 = 0; he4 < particle_number; he4++)  // sum l
			{
				for (he5 = particle_number; he5 < degenerate; he5++)  // sum d
				{
					for (he3 = 0; he3 < degenerate; he3++)
					{
						tempb[he3] = 0;
						tempk[he3] = 0;
					}
					phase = 0;
					tempk[he1] = 1;
					tempk[he4] = 1;
					if (he4 > he1)
					{
						phase++;
					}
					tempb[he2] = 1;
					tempb[he5] = 1;
					if (he5 > he2)
					{
						phase++;
					}
					tempval = Tmatrix[ap][he2][he1][he4] * Tmatrix[he5][bp][ip][jp] - Tmatrix[bp][he2][he1][he4] * Tmatrix[he5][ap][ip][jp];
					val += 0.5 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * tempval;
				}
			}
		}
	}
	//printf("test  %lf\n", val);

	//10th term
	//val = 0;
	for (he1 = 0; he1 < particle_number; he1++)  // sum k
	{
		for (he2 = particle_number; he2 < degenerate; he2++)  // sum c
		{
			for (he4 = 0; he4 < particle_number; he4++)  // sum l
			{
				for (he5 = particle_number; he5 < degenerate; he5++)  // sum d
				{
					for (he3 = 0; he3 < degenerate; he3++)
					{
						tempb[he3] = 0;
						tempk[he3] = 0;
					}
					phase = 0;
					tempk[he1] = 1;
					tempk[he4] = 1;
					if (he4 > he1)
					{
						phase++;
					}
					tempb[he2] = 1;
					tempb[he5] = 1;
					if (he5 > he2)
					{
						phase++;
					}
					tempval = Tmatrix[he2][he5][ip][jp] * Tmatrix[ap][bp][he1][he4];
					val += 0.25 * sgn(phase) * TwobodyOperator(tempk, tempb, 2) * tempval;
				}
			}
		}
	}
	//printf("test  %lf\n", val);


	return val;
}

void InitTmatix()
{
	int ket[20], bar[20], temp;
	int he1, he2, he3, he4;
	int phase;
	//zero
	for (he1 = 0; he1 < degenerate; he1++)   //i
	{
		for (he2 = 0; he2 < degenerate; he2++)  //j
		{
			for (he3 = 0; he3 < degenerate; he3++)  //a
			{
				for (he4 = 0; he4 < degenerate; he4++) //b
				{
					Tmatrix[he1][he2][he3][he4] = 0.;
				}
			}
		}
	}

	for ( he1 = 0; he1 < particle_number; he1++)   //i
	{
		for (he2 = 0; he2 < particle_number; he2++)  //j
		{
			for ( he3 = particle_number; he3 < degenerate; he3++)  //a
			{
				for ( he4 = particle_number; he4 < degenerate; he4++) //b
				{
					for ( temp = 0; temp < degenerate; temp++)
					{
						ket[temp] = 0;
						bar[temp] = 0;
					}
					phase = 0;
					ket[he3] = 1;
					ket[he4] = 1;
					if (he4 > he3)
					{
						phase++;
					}
					bar[he1] = 1;
					bar[he2] = 1;
					if (he2 > he1)
					{
						phase++;
					}
					Tmatrix[he3][he4][he1][he2] = -1. * sgn(phase) * (TwobodyOperator(ket, bar, 2) / (FockMatrix(he1, he1) + FockMatrix(he2, he2) - FockMatrix(he3, he3) - FockMatrix(he4, he4)));
					//printf("  %lf  ", Tmatrix[he3][he4][he1][he2]);
				}
			}
		}
	}
}

void iteration()
{
	
	int ap, bp, ip, jp;
	int h1, h2, h3, h4;
	double Ec, Elast;
	double temp11;
	int ket[20], bar[20], temp,phase, times;
	double tempTmatrix[10][10][10][10];

	times = 0;
	Ec = 1000;
	Elast = 0;
	for ( ap = particle_number; ap < degenerate; ap++)
	{
		for (bp = particle_number; bp < degenerate; bp++)
		{
			for ( ip = 0; ip < particle_number; ip++)
			{
				for ( jp = 0; jp < particle_number; jp++)
				{
					for (temp = 0; temp < degenerate; temp++)
					{
						ket[temp] = 0;
						bar[temp] = 0;
					}
					phase = 0;
					ket[ip] = 1;
					ket[jp] = 1;
					if (jp > ip)
					{
						phase++;
					}
					bar[ap] = 1;
					bar[bp] = 1;
					if (bp > ap)
					{
						phase++;
					}		
					if (Tmatrix[ap][bp][ip][jp] == 0)
					{
						continue;
					}
					Elast += 0.25 *sgn(phase)* TwobodyOperator(ket, bar, 2) * Tmatrix[ap][bp][ip][jp];
				}
			}
		}
	}
	while ( fabs(Elast - Ec) > eps )
	{
		Elast = Ec;
		//printf("Tmaxtix  %f \n", Tmatrix[4][5][0][1] + HbarN(4, 5, 0, 1) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)));
		//printf("Tmaxtix  %f \n", Tmatrix[4][5][1][0] + HbarN(4, 5, 1, 0) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(1, 1) - FockMatrix(0, 0)));

		for (h1 = particle_number; h1 < degenerate; h1++)
		{
			for (h2 = particle_number; h2 < degenerate; h2++)
			{
				for (h3 = 0; h3 < particle_number; h3++)
				{
					for (h4 = 0; h4 < particle_number; h4++)
					{
						tempTmatrix[h1][h2][h3][h4] = 0;
						//Tmatrix[ap][bp][ip][jp] += 1;
						//printf("new Tmaxtix  %f \n", Tmatrix[4][5][0][1] + HbarN(4, 5, 0, 1) / (FockMatrix(4, 4) + FockMatrix(5, 5) - FockMatrix(0, 0) - FockMatrix(1, 1)));
						//temp11 = Tmatrix[h1][h2][h3][h4] + HbarN(h1, h2, h3, h4) / (FockMatrix(h1, h1) + FockMatrix(h2, h2) - FockMatrix(h3, h3) - FockMatrix(h4, h4));
						tempTmatrix[h1][h2][h3][h4] = Tmatrix[h1][h2][h3][h4] - HbarN(h1, h2, h3, h4) / (FockMatrix(h1, h1) + FockMatrix(h2, h2) - FockMatrix(h3, h3) - FockMatrix(h4, h4));
						//printf("Tmaxtix   %d %d %d %d %f %f \n", h1, h2, h3, h4, Tmatrix[h1][h2][h3][h4],temp11);
					}
				}
			}
		}
		for (ap = particle_number; ap < degenerate; ap++)
		{
		for (bp = particle_number; bp < degenerate; bp++)
		{
		for (ip = 0; ip < particle_number; ip++)
		{
			for (jp = 0; jp < particle_number; jp++)
			{
				if (tempTmatrix[ap][bp][ip][jp] != 0)
				{
					Tmatrix[ap][bp][ip][jp] = tempTmatrix[ap][bp][ip][jp];
					//Tmatrix[ap][bp][ip][jp] += HbarN(ap, bp, ip, jp) / (FockMatrix(ap, ap) + FockMatrix(bp, bp) - FockMatrix(ip, ip) - FockMatrix(jp, jp));
					//printf("T %d %d %d %d  %lf \n", ap + 1, bp + 1, ip + 1, jp + 1, Tmatrix[ap][bp][ip][jp]);
					//printf("Tmaxtix  %f \n", HbarN(ap, bp, ip, jp) / (FockMatrix(ap, ap) + FockMatrix(bp, bp) - FockMatrix(ip, ip) - FockMatrix(jp, jp)));
					//printf("Tmaxtix  %f \n", HbarN(ap, bp, jp, ip) / (FockMatrix(ap, ap) + FockMatrix(bp, bp) - FockMatrix(ip, ip) - FockMatrix(jp, jp)));
					//printf("Tmaxtix   %d %d %d %d    %f \n", ap + 1, bp + 1, ip + 1, jp + 1, Tmatrix[ap][bp][ip][jp]);
				}
			}
		}
		}
		}
		Ec = 0;
		for (ap = particle_number; ap < degenerate; ap++)
		{
			for (bp = particle_number; bp < degenerate; bp++)
			{
				for (ip = 0; ip < particle_number; ip++)
				{
					for (jp = 0; jp < particle_number; jp++)
					{
						for (temp = 0; temp < degenerate; temp++)
						{
							ket[temp] = 0;
							bar[temp] = 0;
						}
						phase = 0;
						ket[ip] = 1;
						ket[jp] = 1;
						if (jp > ip)
						{
							phase++;
						}
						bar[ap] = 1;
						bar[bp] = 1;
						if (bp > ap)
						{
							phase++;
						}
						if (Tmatrix[ap][bp][ip][jp] == 0)
						{
							continue;
						}
						Ec += 0.25 * sgn(phase) * TwobodyOperator(ket, bar, 2) * Tmatrix[ap][bp][ip][jp];
					}
				}
			}
		}
		printf("%d   %f   %f     %f\n", times, fabs(Elast - Ec),Ec,Elast);
		times++;
	}
	printf("Correlation energy: %10.10lf\n",Ec);

}
