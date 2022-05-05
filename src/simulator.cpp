#include "FIELD_2D.h"
#include "Watercolor2D.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif

#include "util/ReadPPM.h"
#include "util/FFMPEG_MOVIE.h"

using namespace std;

// resolution of the field
int x_res = 200;
int y_res = 200;

// the field being drawn and manipulated
Eigen::ArrayXXf field(x_res, y_res);

// the simulation object
/* Watercolor2D simulator(x_res, y_res); */
Watercolor2D* simulator = new Watercolor2D(x_res, y_res);

// the resolution of the OpenGL window -- independent of the field resolution
int x_screen_res = 850;
int y_screen_res = 850;

// Text for the title bar of the window
string windowLabel("Computer Generated Watercolor");

FFMPEG_MOVIE movie;

// mouse tracking variables
int x_mouse = -1;
int y_mouse = -1;
int mouse_button = -1;
int mouse_state = -1;
int mouse_modifiers = -1;

// current grid cell the mouse is pointing at
int x_field = -1;
int y_field = -1;

// animate the current runEverytime()?
bool animate = true;

// draw the grid over the field?
bool draw_grid = false;

// print out what the mouse is pointing at?
bool draw_values = true;

// currently capturing frames for a movie?
bool capture_movie = false;

// the current viewer eye position
float eye_center[] = {0.5, 0.5, 1};

// current zoom level into the field
float zoom = 1.0;

// which reaction-diffusion model to simulate
bool use_gs = true;

// forward declare the caching function here so that we can
// put it at the bottom of the file
void runOnce();

// forward declare the timestepping function here so that we can
// put it at the bottom of the file
void runEverytime();

///////////////////////////////////////////////////////////////////////
// Figure out which field element is being pointed at, set x_field and
// y_field to them
///////////////////////////////////////////////////////////////////////
void refreshMouseFieldIndex(int x, int y)
{
  // make the lower left the origin
  y = y_screen_res - y;

  float xNorm = (float)x / x_screen_res;
  float yNorm = (float)y / y_screen_res;

  float half_zoom = 0.5 * zoom;
  float x_world_min = eye_center[0] - half_zoom;
  float x_world_max = eye_center[0] + half_zoom;

  // get the bounds of the field in screen coordinates
  //
  // if non-square textures are ever supported, change the 0.0 and 1.0 below
  float xMin = (0.0 - x_world_min) / (x_world_max - x_world_min);
  float xMax = (1.0 - x_world_min) / (x_world_max - x_world_min);

  float y_world_min = eye_center[1] - half_zoom;
  float y_world_max = eye_center[1] + half_zoom;

  float yMin = (0.0 - y_world_min) / (y_world_max - y_world_min);
  float yMax = (1.0 - y_world_min) / (y_world_max - y_world_min);

  float x_scale = 1.0;
  float y_scale = 1.0;

  if (x_res < y_res)
    x_scale = (float)y_res / x_res;
  if (x_res > y_res)
    y_scale = (float)x_res / y_res;

  // index into the field after normalizing according to screen
  // coordinates
  x_field = x_scale * x_res * ((xNorm - xMin) / (xMax - xMin));
  y_field = y_scale * y_res * ((yNorm - yMin) / (yMax - yMin));

  // clamp to something inside the field
  x_field = (x_field < 0) ? 0 : x_field;
  x_field = (x_field >= x_res) ? x_res - 1 : x_field;
  y_field = (y_field < 0) ? 0 : y_field;
  y_field = (y_field >= y_res) ? y_res - 1 : y_field;
}

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
// dump the field contents to a GL texture for drawing
///////////////////////////////////////////////////////////////////////
void updateTexture(Eigen::ArrayXXf &texture)
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, texture.rows(), texture.cols(), 0, GL_LUMINANCE, GL_FLOAT, texture.data());

  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glEnable(GL_TEXTURE_2D);
}

///////////////////////////////////////////////////////////////////////
// draw a grid over everything
///////////////////////////////////////////////////////////////////////
void drawGrid()
{
  glColor4f(0.1, 0.1, 0.1, 1.0);

  float dx = 1.0 / x_res;
  float dy = 1.0 / y_res;

  if (x_res < y_res)
    dx *= (float)x_res / y_res;
  if (x_res > y_res)
    dy *= (float)y_res / x_res;

  glBegin(GL_LINES);
  for (int x = 0; x < field.rows() + 1; x++) {
    glVertex3f(x * dx, 0, 1);
    glVertex3f(x * dx, 1, 1);
  }
  for (int y = 0; y < field.cols() + 1; y++) {
    glVertex3f(0, y * dy, 1);
    glVertex3f(1, y * dy, 1);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float half_zoom = zoom * 0.5;

  glOrtho(-half_zoom, half_zoom, -half_zoom, half_zoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eye_center[0], eye_center[1], 1, // eye
            eye_center[0], eye_center[1], 0, // center
            0, 1, 0);                      // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  float xLength = 1.0;
  float yLength = 1.0;

  if (x_res < y_res)
    xLength = (float)x_res / y_res;
  if (y_res < x_res)
    yLength = (float)y_res / x_res;

  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(0.0, yLength, 0.0);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(xLength, yLength, 0.0);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(xLength, 0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  // draw the grid, but only if the user wants
  if (draw_grid)
    drawGrid();

  // if there's a valid field index, print it
  if (x_field >= 0 && y_field >= 0 && x_field < field.rows() && y_field < field.cols()) {
    glLoadIdentity();

    // must set color before setting raster position, otherwise it won't take
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

    // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
    float half_zoom = 0.5 * zoom;
    glRasterPos3f(-half_zoom * 0.95, -half_zoom * 0.95, 0);

    // build the field value string
    char buffer[256];
    string fieldValue("(");
    sprintf(buffer, "%i", x_field);
    fieldValue = fieldValue + string(buffer);
    sprintf(buffer, "%i", y_field);
    fieldValue = fieldValue + string(", ") + string(buffer) + string(") = ");
    sprintf(buffer, "%f", field(x_field, y_field));
    fieldValue = fieldValue + string(buffer);

    // draw the grid, but only if the user wants
    if (draw_values)
      printGlString(fieldValue);
  }

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void printCommands()
{
  cout << "=============================================================== " << endl;
  cout << " Field viewer code CPSC 679" << endl;
  cout << "=============================================================== " << endl;
  cout << " q           - quit" << endl;
  cout << " v           - type the value of the cell under the mouse" << endl;
  cout << " g           - throw a grid over everything" << endl;
  cout << " w           - write out a PPM file " << endl;
  cout << " left mouse  - pan around" << endl;
  cout << " right mouse - zoom in and out " << endl;
  cout << " shift left mouse - draw on the grid " << endl;
}

///////////////////////////////////////////////////////////////////////
// Map the arrow keys to something here
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    break;
  case GLUT_KEY_RIGHT:
    break;
  case GLUT_KEY_UP:
    break;
  case GLUT_KEY_DOWN:
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'a':
    animate = !animate;
    break;
  case 'g':
    draw_grid = !draw_grid;
    break;
  case '?':
    printCommands();
    break;
  case 'v':
    draw_values = !draw_values;
    break;
  case 'w': {
    static int count = 0;
    char buffer[256];
    sprintf(buffer, "output_%i.ppm", count);
    /* field.writePPM(buffer); */
    count++;
  } break;
  case 'q':
    movie.streamWriteMovie("simulator.mov");
    exit(0);
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouse_button = button;
  mouse_state = state;
  mouse_modifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // zero out a 10x10 square of chemical b
    for (int sx = max(0, x_field - 5); sx < min(x_res, x_field + 5); sx++)
      for (int sy = max(0, y_field - 5); sy < min(y_res, y_field + 5); sy++)
      {
        simulator->M()(sx,sy) = 1.0f;
        simulator->pigments()[0]->g(sx,sy) = 0.3f;
      }

    // make sure nothing else is called
    return;
  }

  x_mouse = x;
  y_mouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked and moving
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  if (mouse_button == GLUT_LEFT_BUTTON && mouse_state == GLUT_DOWN && mouse_modifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // make sure nothing else is called
    return;
  }

  float x_diff = x - x_mouse;
  float y_diff = y - y_mouse;
  float speed = 0.001;

  if (mouse_button == GLUT_LEFT_BUTTON) {
    eye_center[0] -= x_diff * speed;
    eye_center[1] += y_diff * speed;
  }
  if (mouse_button == GLUT_RIGHT_BUTTON)
    zoom -= y_diff * speed;

  x_mouse = x;
  y_mouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is not clicked and moving
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
  refreshMouseFieldIndex(x, y);
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) {
    runEverytime();

    movie.addFrameGL();
  }

  updateTexture(field);
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(x_screen_res, y_screen_res);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei)x_screen_res, (GLsizei)y_screen_res);

  // set the background color to gray
  glClearColor(0.1, 0.1, 0.1, 0);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // In case the field is rectangular, make sure to center the eye
  if (x_res < y_res) {
    float xLength = (float)x_res / y_res;
    eye_center[0] = xLength * 0.5;
  }
  if (y_res < x_res) {
    float yLength = (float)y_res / x_res;
    eye_center[1] = yLength * 0.5;
  }

  printCommands();

  runOnce();

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glvuWindow();
  return 1;
}

///////////////////////////////////////////////////////////////////////
// This function is called every frame -- do something interesting
// here.
///////////////////////////////////////////////////////////////////////
void runEverytime()
{
  simulator->step();
  /* field = simulator->pigments()[0]->d + simulator.pigments()[0]->g; */
  field = simulator->pigments()[0]->d;
}

///////////////////////////////////////////////////////////////////////
// This is called once at the beginning so you can precache
// something here
///////////////////////////////////////////////////////////////////////
void runOnce()
{
  // larger wet area mask
  for (int y = 0.10 * y_res; y < 0.4 * y_res; y++)
    for (int x = 0.10 * x_res; x < 0.4 * x_res; x++)
      simulator->M()(x, y) = 1.0;
  for (int y = 0.60 * y_res; y < 0.9 * y_res; y++)
    for (int x = 0.60 * x_res; x < 0.9 * x_res; x++)
      simulator->M()(x, y) = 1.0;

  // pigment in center square
  for (int y = 0.20 * y_res; y < 0.30 * y_res; y++)
    for (int x = 0.20 * x_res; x < 0.30 * x_res; x++)
      simulator->pigments()[0]->g(x, y) = 0.05;
  for (int y = 0.70 * y_res; y < 0.80 * y_res; y++)
    for (int x = 0.70 * x_res; x < 0.80 * x_res; x++)
      simulator->pigments()[0]->g(x, y) = 0.05;

  /* float* paper = new float[3 * x_res * y_res]; */
  float* paper = NULL;
  int x, y;
  readPPM("../data/h_1.ppm", x, y, paper);
  simulator->setPaper(paper);
}
