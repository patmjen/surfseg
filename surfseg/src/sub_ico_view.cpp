/*
 * Most of this is from the OBJViewer demo from GEL
 */

#include <iostream>
#include <string>

#include <GL/glew.h>
#include <GL/glut.h>
#include <GEL/GLGraphics/QuatTrackBall.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Vec2i.h>

#include "subdivided_icosahedron.h"

using namespace CGLA;

constexpr int winSizeX = 800, winSizeY = 800;
constexpr int spinTimer = 20;

bool wireframe = false;
bool vertexNormals = true;
bool redoList = true;

GLGraphics::QuatTrackBall* trackBall;

SubdividedIcosahedron mesh(1.0f);

void mouse_motion(int x, int y)
{
    trackBall->roll_ball(Vec2i(x, winSizeY - y));
}

void spin(int x)
{
    trackBall->do_spin();
    glutTimerFunc(spinTimer, spin, 0);
    glutPostRedisplay();
}

void mouse(int btn, int state, int x, int y)
{
    y = winSizeY - y;
    if (state == GLUT_DOWN)
    {
        if (btn == GLUT_LEFT_BUTTON)
            trackBall->grab_ball(GLGraphics::ROTATE_ACTION, Vec2i(x, y));
        else if (btn == GLUT_MIDDLE_BUTTON)
            trackBall->grab_ball(GLGraphics::ZOOM_ACTION, Vec2i(x, y));
        else if (btn == GLUT_RIGHT_BUTTON)
            trackBall->grab_ball(GLGraphics::PAN_ACTION, Vec2i(x, y));
    } else if (state == GLUT_UP)
        trackBall->release_ball();
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
    case 'w':
        std::cout << "Saving to testMesh.obj\n";
        mesh.saveToObj("testMesh.obj");
        break;
    case 'f':
        vertexNormals = !vertexNormals;
        break;
    case 'e':
        mesh.flipEdge(0);
        mesh.computeVertexNormals();
        break;
    case 's':
        mesh.subdivide();
        mesh.computeVertexNormals();
        break;
	case 'r':
		mesh.removeEdge(0);
		mesh.removeInvalid();
		mesh.computeVertexNormals();
		break;
    case 'p':
        wireframe = !wireframe;
        if (wireframe)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    case '+':
        mesh.rescale(mesh.r()*1.1);
        break;
    case '-':
        mesh.rescale(mesh.r() / 1.1);
        break;
    }
    redoList = true;
}

void specialKeyboard(int key, int x, int y)
{
    switch (key) {
    case GLUT_KEY_UP:
        mesh.move(mesh.center() + Vec3f(0, 0.1, 0));
        break;
    case GLUT_KEY_DOWN:
        mesh.move(mesh.center() + Vec3f(0, -0.1, 0));
        break;
    case GLUT_KEY_RIGHT:
        mesh.move(mesh.center() + Vec3f(0.1, 0, 0));
        break;
    case GLUT_KEY_LEFT:
        mesh.move(mesh.center() + Vec3f(-0.1, 0, 0));
        break;
    }
    redoList = true;
}

void render()
{
    static unsigned int list;

    if (redoList) {
        list = glGenLists(1);
        glNewList(list, GL_COMPILE);

        mesh.render(vertexNormals);

        glEndList();
        redoList = false;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    trackBall->set_gl_modelview();

    glCallList(list);
    glutSwapBuffers();
}

int main(int argc, char *argv[])
{
    std::cout << "Setting up OpenGL window...\n";
    // GLUT Init
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(winSizeX, winSizeY);
    int main_window = glutCreateWindow(__FILE__);
    glutDisplayFunc(render);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialKeyboard);

    glutMotionFunc(mouse_motion);
    glutMouseFunc(mouse);

    glewInit();

    // GL Init
    glClearColor(.8f, 0.9f, 1.0f, 0.f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glShadeModel(GL_SMOOTH);

    std::cout << "Computing face normals...\n";
    mesh.computeFaceNormals();
    std::cout << "Computing vertex normals...\n";
    mesh.computeVertexNormals(true);

    std::cout << "Computing bounding sphere...\n";
    Vec3f center;
    float r = mesh.computeBoundingSphere(center);
    r *= 2;
    trackBall = new GLGraphics::QuatTrackBall(center, r, winSizeX, winSizeY);
    glutTimerFunc(spinTimer, spin, 0);

    std::cout << "Starting rendering loop...\n";
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(70, 1.0f, r / 100.0, r*3.0);
    glMatrixMode(GL_MODELVIEW);

    // Pass control to GLUT
    glutMainLoop();

    return 0;
}
