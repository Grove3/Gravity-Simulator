/*  For linux:  gcc gravity.c -lglut -lGL –lGLU -fopenmp -lm -o gravity */
/*  For windows:  -lfreeglut -lopengl32 -lglu32 */


#include <GL/gl.h>				/* Header File For OpenGL Library */
#include <GL/glu.h>				/* Header File For The GLU Library */
#include <GL/freeglut.h>		/* Header File For The GLUT Library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>


/* Define the Window dimensions */
#define WIN_WIDTH	800	/* WINDOW Width in pixels */
#define WIN_HEIGHT	800	/* WINDOW Height in pixels */
#define G 1 /* Universal Gravaty */

/* Define the key to exit the system */
#define ESC	27

/* function prototypes */
void myDraw(void);
void myKey(unsigned char key,int x, int y);

// Simulation parameters
#define VECTOR_DIM 3
#define MAX_BODIES 10000

typedef struct {
	double vector[VECTOR_DIM];
} Vector;

typedef struct {
	double mass;
	Vector position;
	Vector velocity;
	Vector accel;
	double radius;
} Object;


// Global variables
Object bodies[MAX_BODIES];
clock_t Sim_time;
int no_bodies = 32;
int no_collisions = 0;
float ti = 0.f;
float zoom = 0.6;
float xzoom = 0;
float no_threads = 1;
float co_rest = 0.95;
float x_axis = 0;
float y_axis = 0;
double deltaTime = 0.00005;
double sml = 0.3;
char trail = 0;



// Vector function prototypes
void VectorSet(Vector *a, double vec[]);
Vector VectorAdd(Vector a, Vector b);
Vector VectorSubtract(Vector a, Vector b);
double VectorMag(Vector a);
Vector VectorScalar(double scalar, Vector a);
void VectorPrint(Vector a);
Vector VectorUnit(Vector a);
double VectorDot(Vector a, Vector b);



// Gravity function prototypes
void initialise(void);
void simulateGravity(int data);
void collidingBodies(Object *b1, Object *b2);




/* main program.*/
int main(int argc, char *argv[])
{
	int i;

	glutInit(&argc, argv); 				// initialise gl utility system
	glutInitWindowSize(WIN_WIDTH, WIN_HEIGHT); // define window size
	glutCreateWindow ("Gravity Simulator"); 	// create window with title

	initialise();
	glutTimerFunc(deltaTime * 10, simulateGravity, 0); /* Registering Simulator Calc */

	/* Register the callback functions to the current window */
	glutDisplayFunc(myDraw);  		// register function for drawing event
	glutKeyboardFunc(myKey);		// register function for keyboard event

	/* Draw window and loop forever. Any keypress will call mykey() */
	glutMainLoop();

	return (0);
}

void myDraw(void)
{

    int i, j;
    int trails[3] = {0};


    glViewport(0, 0, WIN_WIDTH, WIN_HEIGHT);
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset

    gluPerspective(45.0f, 1, 0.10f, 100.0f);
    //glOrtho( -1, 1, -1, 1, -10, 100 );

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    /* clear the screen */
    glClear(GL_COLOR_BUFFER_BIT);

    /* Position the camera */
    gluLookAt( 0, 0, -2, 0, 0, 0, 0, 1, 0);
    glRotatef(x_axis, 1, 0, 0);
    glRotatef(y_axis, 0, 1, 0);
    glScalef(zoom,zoom,zoom);


	glLineWidth(1);
	glBegin(GL_LINE_LOOP);
		glColor3f(1, 1, 1);
		glVertex3f(0.7, 0.7, 0.7);
		glVertex3f(-0.7, 0.7, 0.7);
		glVertex3f(-0.7, -0.7, 0.7);
		glVertex3f(0.7, -0.7, 0.7);

	glEnd();

	glLineWidth(1);
	glBegin(GL_LINE_LOOP);
		glColor3f(1, 1, 1);
		glVertex3f(0.7, 0.7, -0.7);
		glVertex3f(-0.7, 0.7, -0.7);
		glVertex3f(-0.7, -0.7, -0.7);
		glVertex3f(0.7, -0.7, -0.7);

	glEnd();

	glLineWidth(1);
	glBegin(GL_LINES);
		glColor3f(1, 1, 1);
		glVertex3f(0.7, 0.7, 0.7);
		glVertex3f(0.7, 0.7, -0.7);

		glVertex3f(-0.7, 0.7, 0.7);
		glVertex3f(-0.7, 0.7, -0.7);

		glVertex3f(-0.7, -0.7, 0.7);
		glVertex3f(-0.7, -0.7, -0.7);

		glVertex3f(0.7, -0.7, 0.7);
		glVertex3f(0.7, -0.7, -0.7);

	glEnd();

	glLineWidth(1);
	glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex3f(0.7, -0.7, -0.7);
		glVertex3f(0.1, -0.7, -0.7);

		glColor3f(0, 1, 0);
		glVertex3f(0.7, -0.7, -0.7);
		glVertex3f(0.7, -0.1, -0.7);

		glColor3f(0, 0, 1);
		glVertex3f(0.7, -0.7, -0.7);
		glVertex3f(0.7, -0.7, 0.1);

		glEnd();


    for( i = 0; i < no_bodies; i++)
    {
        /* Draws white spheres */
        glPushMatrix();
        if ( i % 2 == 0)
            {
                glColor3f( (float)i/no_bodies, 0, 0.5 );
            }
        else if ( i % 3 == 0)
            {
                glColor3f( 0, (float)i/no_bodies, 0);
            }
        else
            {
                glColor3f( 0, 1, (float)i/no_bodies );
            }
        glTranslatef(bodies[i].position.vector[0],bodies[i].position.vector[1],bodies[i].position.vector[2]);
        glRotatef(i, 1, 0, 0);
        glRotatef(i, 0, 1, 0);
        glRotatef(i, 0, 0, 1);
        glutSolidSphere( bodies[i].radius * 1.5, 20, 12);
        glPopMatrix();

        if ( trail == 1)
        {
            /*  Draws trajectories lines */
            glBegin(GL_LINES);
            if ( i % 2 == 0)
            {
                glColor3f( (float)i/no_bodies, 0, 0.5 );
            }
            else if ( i % 3 == 0)
            {
                glColor3f( 0, (float)i/no_bodies, 0);
            }
            else
            {
                glColor3f( 0, 1, (float)i/no_bodies );
            }

            glVertex3f( bodies[i].position.vector[0],bodies[i].position.vector[1],bodies[i].position.vector[2] );
            glVertex3f( bodies[i].position.vector[0] + 0.01 * bodies[i].velocity.vector[0],bodies[i].position.vector[1] + 0.01 * bodies[i].velocity.vector[1], bodies[i].position.vector[2] + + 0.01 * bodies[i].velocity.vector[2] );

		}

        glEnd();

    }

    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();

    glOrtho( -1, 1, -1, 1, -1, 1 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	//Text for Window

	char str1[100], str2[100], str3[100], str4[100], str5[100], str6[100], str7[100];
	ti = ((double)( clock() - Sim_time ))/CLOCKS_PER_SEC;

	glColor3f(1, 1, 1);
	glRasterPos2f(-0.95, 0.95);
	sprintf(str1, "Trace %s (Press 1 to turn %s)", trail == 0 ? "OFF" : "ON", trail == 1 ? "OFF" : "ON");
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str1);

	glColor3f(1, 1, 1);
	glRasterPos2f(-0.95, 0.91);
	sprintf(str7, "Rotate around X-axis (q - Q), Rotate around Y-axis (w - W)");
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str7);

	glRasterPos2f(-0.95, -0.79);
	sprintf(str2, "Simulation (r) =  %.4fs, (theta = 0.00, phi = 0.00), Zoom (+z/-x) = %.2f", ti, (zoom + 0.4));
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str2);

	glRasterPos2f(-0.95, -0.83);
	sprintf(str3, "Number of bodies (+b/-v): %i", no_bodies);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str3);

	glRasterPos2f(-0.95, -0.87);
	sprintf(str4, "Number of threads (+y/-t): %.1f", no_threads);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str4);

	glRasterPos2f(-0.95, -0.91);
	sprintf(str5, "Number of collisions: %i", no_collisions);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str5);

	glRasterPos2f(-0.95, -0.95);
	sprintf(str6, "Coefficient of restitution (+k/-l): %.2f", co_rest);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12, str6);

	glFlush();


}


void myKey(unsigned char key, int x, int y)
{
	int win;
	switch (key)
	{
		case 'z':				//zoom in
			zoom += 0.1;
			break;
		case 'x':				//zoom out
			zoom -= 0.1;
			break;
		case 'q':				//Rotate X - axis
			x_axis += 1;
			break;
		case 'Q':				//Rotate X - axis
			x_axis -= 1;
			break;
		case 'w':				//Rotate Y - axis
			y_axis += 1;
			break;
		case 'W':				//Rotate Y - axis
			y_axis -= 1;
			break;
		case 'b':				//increase number of bodies
			no_bodies *= 2;
			initialise();
			break;
		case 'v':				//decrease number of bodies
			no_bodies /= 2;
			if (no_bodies < 2)
				no_bodies = 2;
			initialise();
			break;
		case 'y':				//increase number of threads
			no_threads ++;
			break;
		case 't':				//decrease number of threads
			no_threads --;
			if(no_threads < 1)
				no_threads = 1;
			break;
		case 'k':				//Coefficient of restitution
			co_rest += 0.05;
			if(co_rest > 1)
				co_rest = 1;
			break;
		case 'l':				//Coefficient of restitution
			co_rest -= 0.05;
			if(co_rest < 0)
				co_rest = 0;
			break;
		case 'r':				//Rest simulation
			initialise();
			break;
		case '1':				//Turn on trails
			trail = !trail;
			break;
	}

	/* check if keypressed was the ESC key then close window and leave program */
	if (key == ESC) {
		win=glutGetWindow() ; 		/* get identifier for current window */
		glutDestroyWindow(win) ; 	/* close the current window */
		exit(0) ; 					/* leave the program */
	}
}

// gravity functions
/*	Initialising the positions, velocity, and acceleration of objects. Also reset the simulation time */
void initialise(void)
{
	int i, j;


	for(i = 0; i < no_bodies; i++)
		{
			bodies[i].mass = 50;
			bodies[i].radius = bodies[i].mass * 0.0002;
			for(j = 0; j < VECTOR_DIM; j++)
			{
			  bodies[i].position.vector[j] = (float)rand() / ((float)RAND_MAX/2) - 1;
			  bodies[i].velocity.vector[j] = 0;
			  bodies[i].accel.vector[j] = 0;
			}
		}
		Sim_time = clock();
		no_collisions = 0;

}

/*	Computing forces, acceleration, velocity, positions and collisions of all objects in a time step.*/
void simulateGravity(int data)
{
	int k;
    double zeroVector[VECTOR_DIM] = {0};
    Vector newPos[no_bodies];

    omp_set_num_threads( no_threads ); /* Setting number of threads to create */
    #pragma omp parallel /* Starting parallel execution */
    {
        int i, j, id;
        Vector object_a = {0};
        Vector object_ba;
        double object_ba_mag;
        Vector object_ba_unit;
        double force_mag;
        Vector object_b;
        Vector accel;
        Vector vel_object_a2;
        Vector delta_object_a2;
        Vector new_object_a2;


        no_threads = omp_get_num_threads(); /* Number of threads created */
        id = omp_get_thread_num(); /* Current thread id */
	for(i = id; i < no_bodies; i = i + no_threads) /* Distributing bodies for each thread */
        {
            VectorSet( &object_a, zeroVector);
            for(j = 0; j < no_bodies; j++)
            {
                if( i == j ) continue;
                {
					object_ba = VectorSubtract(bodies[j].position, bodies[i].position);

					object_ba_mag = VectorMag(object_ba);

					object_ba_unit = VectorUnit(object_ba);

					force_mag = ((bodies[j].mass * bodies[i].mass * G) / (object_ba_mag * object_ba_mag) + (sml * sml));

					object_b = VectorScalar(force_mag, object_ba_unit);

                	object_a = VectorAdd(object_a, object_b); /* Calculate gravitational force and add to previously calculated force due to other bodies */
                }
            }

			accel = VectorScalar(1.0/bodies[i].mass, object_a);

			vel_object_a2 = VectorAdd(bodies[i].velocity, VectorScalar(deltaTime, accel));

			delta_object_a2 = VectorAdd(VectorScalar(deltaTime, bodies[i].velocity), VectorScalar(0.5 * deltaTime * deltaTime, accel));

			new_object_a2 = VectorAdd(bodies[i].position, delta_object_a2);

			bodies[i].velocity = vel_object_a2;

			bodies[i].accel = accel;

            newPos[i] = new_object_a2; /* Calculate displacement vector */
        }
    }

    /* Apply new positions */
    for( k = 0; k < no_bodies; k++)
    {
        bodies[k].position = newPos[k];
    }

	omp_set_num_threads( no_threads ); /* Setting number of threads to create */
    #pragma omp parallel /* Starting parallel execution */
    {

        int i, j, id;
        double distance;

        no_threads = omp_get_num_threads();
        id = omp_get_thread_num();

        for( i = id; i < no_bodies; i = i + no_threads)
        {
            for( j = 0; j < no_bodies; j++)
            {
                if( i == j ) continue;

                distance = VectorMag( VectorSubtract( bodies[i].position, bodies[j].position ));

                if( distance <= (bodies[i].radius + bodies[j].radius))
                {
                    collidingBodies( &bodies[i], &bodies[j]);
                    no_collisions++;
                }
            }
        }

    }

    glutPostRedisplay();
    glutTimerFunc(10, simulateGravity, 0);

}

/*	Simulate perfectly elastic collisions. 
	Compute and update the new velocities of two objects right after the collision
	Includes coefficient of restitution*/
void collidingBodies(Object *a, Object *b)
{
    Vector x12;
    Vector v12;
    Vector newV1;

    Vector x21;
    Vector v21;
    Vector newV2;

    x12 = VectorSubtract( a->position, b->position);
    v12 = VectorSubtract( a->velocity, b->velocity);

    newV1 = VectorSubtract( a->velocity, VectorScalar( ( 2 * b->mass / ( a->mass + b->mass) * ( VectorDot( v12, x12 ) / ( VectorMag(x12) * VectorMag(x12)) )), x12 ) );

    x21 = VectorSubtract( b->position, a->position);
    v21 = VectorSubtract( b->velocity, a->velocity);

    newV2 = VectorSubtract( b->velocity, VectorScalar( ( 2 * a->mass / ( a->mass + b->mass) * ( VectorDot( v21, x21 ) / ( VectorMag(x21) * VectorMag(x21)) )), x21 ) );

    a->velocity = VectorScalar(co_rest, newV1);
    b->velocity = VectorScalar(co_rest, newV2);

    a->position = VectorAdd( a->position, VectorScalar(deltaTime, a->velocity));
    b->position = VectorAdd( b->position, VectorScalar(deltaTime, b->velocity));

}

// Vector functions
void VectorSet(Vector *a, double vec[])
{
	int i;

	for (i = 0; i < VECTOR_DIM; i++)
		a->vector[i] = vec[i];
}

/* Add vectors and save the result in c */
Vector VectorAdd(Vector a, Vector b)
{
	Vector c;
	int i;

	for (i = 0; i < VECTOR_DIM; i++)
		c.vector[i] = a.vector[i] + b.vector[i];

	return (c);
}

/* Subtract vectors and save the result in c */
Vector VectorSubtract(Vector a, Vector b)
{
	Vector c;
	int i;

	for (i = 0; i < VECTOR_DIM; i++)
		c.vector[i] = a.vector[i] - b.vector[i];

	return (c);
}

/* Compute the unit vector of a and save the result in b */
Vector VectorUnit(Vector a)
{
	Vector b;
	int i;
    double mag;

    mag = VectorMag(a);

	for (i = 0; i < VECTOR_DIM; i++)
		b.vector[i] = a.vector[i] / mag;

	return (b);
}

/* Computes the magnitude of vector a and return it */
double VectorMag(Vector a)
{
	double mag = 0;
	int i;

	for (i = 0; i < VECTOR_DIM; i++)
		mag += a.vector[i] * a.vector[i];

	mag = sqrt(mag);

	return (mag);
}

/* Multiplys vector a by a scalar number and save the result in b */
Vector VectorScalar(double scalar, Vector a)
{
	Vector b;
	int i;

    for (i = 0; i < VECTOR_DIM; i++)
        b.vector[i] = a.vector[i] * scalar;

	return (b);
}

void VectorPrint(Vector a)
{
	int i;

	printf("[ ");

	for (i = 0; i < VECTOR_DIM; i++)
		printf("%lg, ", a.vector[i]);

	printf(" ] ");
}

/* 	Performs a dot product between vectors a and b and save the result in dot */
double VectorDot(Vector a, Vector b)
{
	double dot = 0;
	int i;

    for (i = 0; i < VECTOR_DIM; i++)
        dot += a.vector[i] * b.vector[i];

	return (dot);
}


