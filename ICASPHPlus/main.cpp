#include "ICASPHPlus.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QCoreApplication::setAttribute(Qt::AA_UseDesktopOpenGL);

	QApplication a(argc, argv);

	QSurfaceFormat format;
	format.setRenderableType(QSurfaceFormat::OpenGL);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setVersion(3, 3);
	//format.setOption(QSurfaceFormat::DebugContext);
    #ifdef    GL_DEBUG
	format.setOption(QSurfaceFormat::DebugContext);
    #endif // GL_DEBUG
	ICASPHPlus w;
	w.show();
	return a.exec();
}
