#include "udf.h"
#include "dynamesh_tools.h"
#include "string.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#define pi 3.14159265
#define ID 5			/*ID of foil */
#define frequency 1
#define density 997		/* kg/m2 */

/*change*/
#define Nn 800			 // number of nodes horizontal on the body (important to ensure during case creation
real ds = 0.00011875;	/*based on generated mesh*/ //Lsensitive
real len = 0.095;		//Lsensitive

//real area = 0.00032;	/* m2*/
//real mass = 0.31904;	/*kg*/


int posID;
int i,j;
real r;
real J[3];

// functions 
real Outys;
real Outws;

real w0 = 0.64;
real w1 = 0.16;
real ws1, ws2;

real cgx;

//swimmer motion variables

real dx, dy, dys;
real dx_temp, dx_sum;
real dy_temp;
real y, y2;
real xr, yr;
real y1_temp, y2_temp;
real alpha;

real temp_x1, temp_y1, temp_x2, temp_y2;

real cgx_temp, cgy_temp;

//Centre of gravity variables
real A, SumA, SumAx, SumAy;

//kinematic variables
real c1 = -0.01; 
real c2 = 0.25; 

//global motion
real acc[3], vel[3], pos[3];
real ang_acc[3], ang_vel[3], angle[3];

real f[3], m[3], cg[3];

/*Runge Kutta*/
real k1, k2, k3, k4;


double ys(double InX)
{

	//Lsensitive
	if (CURRENT_TIME < 1 / frequency)
	{
		return Outys = 0.0125 * (InX / len + 0.03125) / 1.03125 * sin(2 * pi *(InX / len - CURRENT_TIME * frequency))*0.5*(1 - cos(pi*CURRENT_TIME*frequency));
		//return Outys = c1 * (((InX / len) + c2) / (1 + c2) * sin(2 * pi *(InX / len - CURRENT_TIME * frequency))*0.5*(1 - cos(pi*CURRENT_TIME*frequency)));
	}
	else
	{
		return Outys = 0.0125 * (InX / len + 0.03125) / 1.03125 * sin(2 * pi *(InX / len - CURRENT_TIME * frequency));
		//return Outys = c1 * (((InX / len) + c2) / (1 + c2) * sin(2 * pi *(InX / len - CURRENT_TIME * frequency)));
	}
}
double ws(double InS)
{

	//Lsensitive
	if (InS < 0.004)
	{
		Outws = sqrt(2*0.004 * ds- ds*ds);			/*Change*/
	}
	else
	{
		Outws = 0.004*(len-ds)/(len-0.004);
	}
	return Outws;
}

DEFINE_INIT(init, d)
{

	/*** Calculate initial centre of gravity ***/
	for (i = 0; i < Nn; i++)
	{
		//Calculate area of element between node i and i+1
		ws1 = ws(i*ds);
		ws2 = ws((i + 1)*ds);

		A = (ws1 + ws2) / 2 * ds;
		SumA = SumA + A;				//Result checked: correct!

		dx_sum = 0.0; // loop to calculate dx
		for (j = 0; j < i; j++)
		{
			y1_temp = ys(j*ds);				//ys(s)
			y2_temp = ys((j + 1)*ds);		//ys(s+ds)
			dy_temp = y2_temp - y1_temp;	//dy = y2-y1

			dx = sqrt(ds*ds - dy_temp * dy_temp);

			dx_sum = dx_sum + dx;
		}

		alpha = asin((ys((i + 1)*ds) - ys(i*ds)) / ds);

		cgx = dx_sum + cos(alpha)*(ws1 + 2 * ws2) / (3 * (ws2 + ws1)) *ds;

		SumAx = SumAx + A * cgx; //(cgx * cos(angle[2]) - ys(cgx) * sin(angle[2]));				//	(x * cos(angle[2]) - y * sin(angle[2]);		
		SumAy = SumAy + A * ys(cgx); //(cgx * sin(angle[2]) + ys(cgx) * cos(angle[2]));			//  (x * sin(angle[2]) + y * cos(angle[2]);		//define xi and yi based on body orientation towards the swimmer

		//Calculate the mass moment of Inertia

		r = sqrt(	(cg[0] - pos[0] - dx_sum)*	(cg[0] - pos[0] - dx_sum) +	(cg[1] - pos[1] - ys(dx_sum) )	*	(cg[1] - pos[1] - ys(dx_sum))); // ERROR


		J[2] = J[2] + A * density* r*r; // Area* mass *r^2
	}


	cg[0] = SumAx / SumA + pos[0];
	cg[1] = SumAy / SumA + pos[1];
	cg[2] = 0.0;

#if PARALLEL
	if (I_AM_NODE_ZERO_P)
	{
		PRF_CSEND_REAL(node_host, cg, 3, myid);
		PRF_CSEND_REAL(node_host, J, 3, myid);
		PRF_CSEND_REAL(node_host, pos, 3, myid);
	}
#if RP_HOST

		PRF_CRECV_REAL(node_zero, cg, 3, node_host);
		PRF_CRECV_REAL(node_zero, J, 3, node_host);
		PRF_CRECV_REAL(node_zero, pos, 3, node_host);

#endif
#endif

}
DEFINE_EXECUTE_AT_END(execute_at_end)
{
#if !RP_NODE
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;
	FILE *fp6;
#endif

#if !RP_HOST

	Domain *d = Get_Domain(1);
	Thread *tf = Lookup_Thread(d, ID);

	//Calculate center of mass cgx and cgy coordinate
	SumA = 0.0;
	SumAx = 0.0;
	SumAy = 0.0;
	J[2] = 0.0;

	for (i = 0; i < Nn; i++)
	{
		//Calculate area of element between node i and i+1
		ws1 = ws(i*ds);
		ws2 = ws((i+1)*ds);
		
		A = (ws1 + ws2)/2 * ds;
		SumA = SumA + A;				//Result checked: correct!

		dx_sum = 0.0; // loop to calculate dx
		for (j = 0; j < i; j++)
		{
			y1_temp = ys(j*ds);				//ys(s)
			y2_temp = ys((j + 1)*ds);		//ys(s+ds)
			dy_temp = y2_temp - y1_temp;	//dy = y2-y1

			dx = sqrt(ds*ds - dy_temp * dy_temp);

			dx_sum = dx_sum + dx;
		}

		alpha = asin((ys((i+1)*ds)-ys(i*ds)) / ds);

		cgx = dx_sum + cos(alpha)*(ws1+2*ws2) / (3 * (ws2 + ws1)) *ds;

		SumAx = SumAx + A * cgx;		//  (cgx * cos(angle[2]) - ys(cgx) * sin(angle[2]));				//	(x * cos(angle[2]) - y * sin(angle[2]);		
		SumAy = SumAy + A * ys(cgx);	// (cgx * sin(angle[2]) + ys(cgx) * cos(angle[2]));			//  (x * sin(angle[2]) + y * cos(angle[2]);		//define xi and yi based on body orientation towards the swimmer

		//Calculate the mass moment of Inertia

		r = sqrt(((cg[0] - pos[0]) - dx_sum)*((cg[0] - pos[0]) - dx_sum) + ((cg[1] - pos[1]) - ys(dx_sum))*((cg[1] - pos[1]) - ys(dx_sum)));
		J[2] = J[2] + A * density* r*r; // Area* mass *r^2

	}


	Message("Area %f radius %f MoI %f \n", SumA, r, J[2]);

	N3V_S(f, =, 0.0);
	N3V_S(m, =, 0.0);


	cg[0] = SumAx / SumA + pos[0];
	cg[1] = SumAy / SumA + pos[1];
	cg[2] = 0.0;
	
	//calculate moment of Inertia
	Compute_Force_And_Moment(d, tf, cg, f, m, FALSE);							/*True if called on honst otherwise FALSE*/



	N3V_S(acc, =, 0.0);															/*reset acc to zero before calculating*/

	/*angular motion*/

	ang_acc[0] = 0.0;
	ang_acc[1] = 0.0;
	ang_acc[2] = m[2] / J[2];

	ang_vel[0] = 0.0;	
	ang_vel[1] = 0.0;
	ang_vel[2] = ang_vel[2] + ang_acc[2] * CURRENT_TIMESTEP; /*cacluate omega_z acceleration*/

	angle[0] = 0.0;
	angle[1] = 0.0;
	angle[2] = angle[2] + ang_vel[2] * CURRENT_TIMESTEP;


	/*linear motion*/
	acc[0] =  f[0] / (density * SumA);												/*cacluate x acceleration*/
	acc[1] =  f[1] / (density * SumA);												/*cacluate y acceleration*/
	acc[2] =  f[2] / (density * SumA);												/*cacluate z acceleration*/

	vel[0] = vel[0] + acc[0] * CURRENT_TIMESTEP;
	vel[1] = vel[1] + acc[1] * CURRENT_TIMESTEP;
	vel[2] = vel[2] + acc[2] * CURRENT_TIMESTEP;
	 
	pos[0] = pos[0] + vel[0] * CURRENT_TIMESTEP;
	pos[1] = pos[1] + vel[1] * CURRENT_TIMESTEP;
	pos[2] = pos[2] + vel[2] * CURRENT_TIMESTEP;

	//Message("Theta %f \n", angle[2]);
	 

#endif
#if PARALLEL
	if (I_AM_NODE_ZERO_P)
	{
		PRF_CSEND_REAL(node_host, cg, 3, myid);
		PRF_CSEND_REAL(node_host, f, 3, myid);
		PRF_CSEND_REAL(node_host, m, 3, myid);
		PRF_CSEND_REAL(node_host, J, 3, myid);
		PRF_CSEND_REAL(node_host, acc, 3, myid);
		PRF_CSEND_REAL(node_host, vel, 3, myid);
		PRF_CSEND_REAL(node_host, pos, 3, myid);
		PRF_CSEND_REAL(node_host, ang_acc, 3, myid);
		PRF_CSEND_REAL(node_host, ang_vel, 3, myid);
		PRF_CSEND_REAL(node_host, angle, 3, myid);


	}
	#if RP_HOST
		PRF_CRECV_REAL(node_zero, cg, 3, node_host);
		PRF_CRECV_REAL(node_zero, f, 3, node_host);
		PRF_CRECV_REAL(node_zero, m, 3, node_host);
		PRF_CRECV_REAL(node_zero, J, 3, node_host);
		PRF_CRECV_REAL(node_zero, acc, 3, node_host);
		PRF_CRECV_REAL(node_zero, vel, 3, node_host);
		PRF_CRECV_REAL(node_zero, pos, 3, node_host);
		PRF_CRECV_REAL(node_zero, ang_acc, 3, node_host);
		PRF_CRECV_REAL(node_zero, ang_vel, 3, node_host);
		PRF_CRECV_REAL(node_zero, angle, 3, node_host);

	#endif
#endif

#if !RP_NODE

	fp = fopen("0linear.txt", "a");
	fprintf(fp, "%f\t ", CURRENT_TIME);
	fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t", acc[0], acc[1], vel[0] , vel[1], pos[0], pos[1] );
	fprintf(fp, "\n");
	fclose(fp);

	fp1 = fopen("1CG.txt", "a");
	fprintf(fp1, "%f\t ", CURRENT_TIME);
	fprintf(fp1, "%f\t%f\t%f\t%f\t%f\t%f\t", cg[0], cg[1], cg[2], J[2], ys(0.0), ys(len));
	fprintf(fp1, "\n");
	fclose(fp1);

	fp2 = fopen("2momentum.txt", "a");
	fprintf(fp2, "%f\t ", CURRENT_TIME);
	fprintf(fp2, "%f\t%f\t%f\t", m[0], m[1], m[2]);
	fprintf(fp2, "\n");
	fclose(fp2);

	fp3 = fopen("3force.txt", "a");
	fprintf(fp3, "%f\t ", CURRENT_TIME);
	fprintf(fp3, "%f\t%f\t%f\t", f[0], f[1], f[2]);
	fprintf(fp3, "\n");
	fclose(fp3);

	fp4 = fopen("4angular.txt", "a");
	fprintf(fp4, "%f\t ", CURRENT_TIME);
	fprintf(fp4, "%f\t%f\t%f\t", ang_acc[2], ang_vel[2], angle[2]);
	fprintf(fp4, "\n");
	fclose(fp4);

	fp5 = fopen("5Analysis.txt", "a");
	fprintf(fp5, "%f\t ", CURRENT_TIME*frequency);
	fprintf(fp5, "%f\t%f\t%f\t%f\t%f\t", -vel[0] / len, -vel[1] / len, -(vel[0] * cos(-angle[2]) - vel[1] * sin(-angle[2])) / len, -(vel[0] * sin(-angle[2]) + vel[1] * cos(-angle[2])) / len, ys(len));
	fprintf(fp5, "\n");
	fclose(fp5);

#endif

}
DEFINE_GRID_MOTION(foil, domain, dt, time, dtime)
{
	Thread *tf = DT_THREAD(dt);
	face_t f;
	Node *v;
	int n;

	SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));
	begin_f_loop(f, tf)
	{
		if PRINCIPAL_FACE_P(f, tf)
		{

			f_node_loop(f, tf, n)
			{
				v = F_NODE(f, tf, n);							// obtain global face node number
				if (NODE_POS_NEED_UPDATE(v))
				{
					NODE_POS_UPDATED(v);

					if (N_TIME < 1)
					{
						N_UDMI(v, 0) = NODE_Y(v) - pos[1];				// ws - defined by created case
						N_UDMI(v, 1) = (NODE_X(v) - pos[0]) / len * Nn;	// Position ID //200 nodes along the ex achse in created case

						//Message("N_UDMI(v, 1)!!!!!!!!!!!! %f \n", N_UDMI(v, 1));
					}

					posID = round(N_UDMI(v, 1));			//round() important otherwise some points are placed wrong

					//Message("posID %i PosID_raw %f X_cord %f \n", posID, N_UDMI(v, 1), NODE_X(v));

						// calculate dy
					dys = ys(posID*ds);

					//calculate dx
					dx_sum = 0.0;
					for (i = 0; i < posID; i++)
					{
						y1_temp = ys(i*ds);				//ys(s)
						y2_temp = ys((i + 1)*ds);		//ys(s+ds)
						dy_temp = y2_temp - y1_temp;

						dx = sqrt(ds*ds - dy_temp * dy_temp);
						dx_sum = dx_sum + dx;
					}

					temp_x1 = dx_sum + pos[0];
					temp_y1 = N_UDMI(v, 0) + dys + pos[1];		//initial y + dys function increment + global position y

					NODE_Y(v) = temp_y1;

					//calculate alpha
					y= ys(posID*ds);				//ys(s)
					y2 = ys((posID + 1)*ds);		//ys(s+ds)
					dy = y2 - y;

					alpha = asin(dy / ds);

					xr = dx_sum + pos[0];
					yr = ys(posID*ds) + pos[1];

					temp_x2 = xr + (temp_x1 - xr)*cos(alpha) - (temp_y1 - yr)*sin(alpha);
					temp_y2 = yr + (temp_x1 - xr)*sin(alpha) + (temp_y1 - yr)*cos(alpha);

					NODE_X(v) = cg[0] + (temp_x2 - cg[0]) * cos(angle[2]) - (temp_y2 - cg[1]) * sin(angle[2]);
					NODE_Y(v) = cg[1] + (temp_x2 - cg[0]) * sin(angle[2]) + (temp_y2 - cg[1]) * cos(angle[2]);


				}
			}
		}
	}
	end_f_loop(f, tf);
	
}