/********************************************************************************
** Form generated from reading UI file 'ICASPHPlus.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ICASPHPLUS_H
#define UI_ICASPHPLUS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ICASPHPlusClass
{
public:
    QAction *actionStart;
    QAction *actionStop;
    QAction *actionSavePic;
    QAction *actionshowBoundaryPar;
    QAction *actionIISPH;
    QAction *actionWCSPH;
    QAction *actionSingleBrak;
    QAction *actionDoubleDamBreak;
    QAction *actionStaticRigidInteraction;
    QAction *actionSurfaceSimulation;
    QAction *actionOn;
    QAction *actionOff;
    QWidget *centralWidget;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *BaseLayout;
    QLabel *SPHTypeName;
    QFrame *line;
    QMenuBar *menuBar;
    QMenu *menuSPH;
    QMenu *menuSPH_Interaction;
    QMenu *menuDistanceField;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ICASPHPlusClass)
    {
        if (ICASPHPlusClass->objectName().isEmpty())
            ICASPHPlusClass->setObjectName(QStringLiteral("ICASPHPlusClass"));
        ICASPHPlusClass->resize(1060, 795);
        actionStart = new QAction(ICASPHPlusClass);
        actionStart->setObjectName(QStringLiteral("actionStart"));
        actionStart->setCheckable(true);
        QIcon icon;
        icon.addFile(QStringLiteral("Resources/start.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionStart->setIcon(icon);
        actionStop = new QAction(ICASPHPlusClass);
        actionStop->setObjectName(QStringLiteral("actionStop"));
        actionStop->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QStringLiteral("Resources/stop.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionStop->setIcon(icon1);
        actionSavePic = new QAction(ICASPHPlusClass);
        actionSavePic->setObjectName(QStringLiteral("actionSavePic"));
        actionSavePic->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QStringLiteral("Resources/savepic.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSavePic->setIcon(icon2);
        actionshowBoundaryPar = new QAction(ICASPHPlusClass);
        actionshowBoundaryPar->setObjectName(QStringLiteral("actionshowBoundaryPar"));
        actionshowBoundaryPar->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QStringLiteral("Resources/showBoundaryP.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionshowBoundaryPar->setIcon(icon3);
        actionIISPH = new QAction(ICASPHPlusClass);
        actionIISPH->setObjectName(QStringLiteral("actionIISPH"));
        actionWCSPH = new QAction(ICASPHPlusClass);
        actionWCSPH->setObjectName(QStringLiteral("actionWCSPH"));
        actionSingleBrak = new QAction(ICASPHPlusClass);
        actionSingleBrak->setObjectName(QStringLiteral("actionSingleBrak"));
        actionDoubleDamBreak = new QAction(ICASPHPlusClass);
        actionDoubleDamBreak->setObjectName(QStringLiteral("actionDoubleDamBreak"));
        actionStaticRigidInteraction = new QAction(ICASPHPlusClass);
        actionStaticRigidInteraction->setObjectName(QStringLiteral("actionStaticRigidInteraction"));
        actionSurfaceSimulation = new QAction(ICASPHPlusClass);
        actionSurfaceSimulation->setObjectName(QStringLiteral("actionSurfaceSimulation"));
        actionOn = new QAction(ICASPHPlusClass);
        actionOn->setObjectName(QStringLiteral("actionOn"));
        actionOff = new QAction(ICASPHPlusClass);
        actionOff->setObjectName(QStringLiteral("actionOff"));
        centralWidget = new QWidget(ICASPHPlusClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QStringLiteral("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(0, 80, 1061, 661));
        BaseLayout = new QVBoxLayout(verticalLayoutWidget);
        BaseLayout->setSpacing(6);
        BaseLayout->setContentsMargins(11, 11, 11, 11);
        BaseLayout->setObjectName(QStringLiteral("BaseLayout"));
        BaseLayout->setContentsMargins(0, 0, 0, 0);
        SPHTypeName = new QLabel(centralWidget);
        SPHTypeName->setObjectName(QStringLiteral("SPHTypeName"));
        SPHTypeName->setGeometry(QRect(20, 10, 431, 51));
        QFont font;
        font.setFamily(QStringLiteral("Segoe Print"));
        font.setPointSize(22);
        SPHTypeName->setFont(font);
        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setGeometry(QRect(470, 0, 20, 81));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        ICASPHPlusClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ICASPHPlusClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1060, 23));
        menuSPH = new QMenu(menuBar);
        menuSPH->setObjectName(QStringLiteral("menuSPH"));
        menuSPH_Interaction = new QMenu(menuBar);
        menuSPH_Interaction->setObjectName(QStringLiteral("menuSPH_Interaction"));
        menuDistanceField = new QMenu(menuBar);
        menuDistanceField->setObjectName(QStringLiteral("menuDistanceField"));
        ICASPHPlusClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ICASPHPlusClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ICASPHPlusClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ICASPHPlusClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ICASPHPlusClass->setStatusBar(statusBar);

        menuBar->addAction(menuSPH->menuAction());
        menuBar->addAction(menuSPH_Interaction->menuAction());
        menuBar->addAction(menuDistanceField->menuAction());
        menuSPH->addAction(actionIISPH);
        menuSPH->addAction(actionWCSPH);
        menuSPH_Interaction->addAction(actionSingleBrak);
        menuSPH_Interaction->addAction(actionDoubleDamBreak);
        menuSPH_Interaction->addAction(actionStaticRigidInteraction);
        menuSPH_Interaction->addAction(actionSurfaceSimulation);
        menuDistanceField->addAction(actionOn);
        menuDistanceField->addAction(actionOff);
        mainToolBar->addAction(actionStart);
        mainToolBar->addAction(actionStop);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionSavePic);
        mainToolBar->addAction(actionshowBoundaryPar);

        retranslateUi(ICASPHPlusClass);
        QObject::connect(actionStart, SIGNAL(triggered()), ICASPHPlusClass, SLOT(onStart()));
        QObject::connect(actionStop, SIGNAL(triggered()), ICASPHPlusClass, SLOT(onStop()));
        QObject::connect(actionSavePic, SIGNAL(triggered()), ICASPHPlusClass, SLOT(isSavePic()));
        QObject::connect(actionshowBoundaryPar, SIGNAL(triggered()), ICASPHPlusClass, SLOT(isShowBoundaryPartilces()));
        QObject::connect(actionIISPH, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeIISPH()));
        QObject::connect(actionWCSPH, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeWCSPH()));
        QObject::connect(actionDoubleDamBreak, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeDoubleBreak()));
        QObject::connect(actionStaticRigidInteraction, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeStaticRigid()));
        QObject::connect(actionSurfaceSimulation, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeSurface()));
        QObject::connect(actionSingleBrak, SIGNAL(triggered()), ICASPHPlusClass, SLOT(changeSingleBreak()));
        QObject::connect(actionOn, SIGNAL(triggered()), ICASPHPlusClass, SLOT(onDrawDistanceField()));
        QObject::connect(actionOff, SIGNAL(triggered()), ICASPHPlusClass, SLOT(offDrawDistanceField()));

        QMetaObject::connectSlotsByName(ICASPHPlusClass);
    } // setupUi

    void retranslateUi(QMainWindow *ICASPHPlusClass)
    {
        ICASPHPlusClass->setWindowTitle(QApplication::translate("ICASPHPlusClass", "ICASPHPlus", Q_NULLPTR));
        actionStart->setText(QApplication::translate("ICASPHPlusClass", "Start", Q_NULLPTR));
        actionStop->setText(QApplication::translate("ICASPHPlusClass", "Stop", Q_NULLPTR));
        actionSavePic->setText(QApplication::translate("ICASPHPlusClass", "SavePic", Q_NULLPTR));
        actionshowBoundaryPar->setText(QApplication::translate("ICASPHPlusClass", "showBoundaryPar", Q_NULLPTR));
        actionIISPH->setText(QApplication::translate("ICASPHPlusClass", "IISPH", Q_NULLPTR));
        actionWCSPH->setText(QApplication::translate("ICASPHPlusClass", "WCSPH", Q_NULLPTR));
        actionSingleBrak->setText(QApplication::translate("ICASPHPlusClass", "SingleDamBreak", Q_NULLPTR));
        actionDoubleDamBreak->setText(QApplication::translate("ICASPHPlusClass", "DoubleDamBreak", Q_NULLPTR));
        actionStaticRigidInteraction->setText(QApplication::translate("ICASPHPlusClass", "StaticRigidInteraction", Q_NULLPTR));
        actionSurfaceSimulation->setText(QApplication::translate("ICASPHPlusClass", "SurfaceSimulation", Q_NULLPTR));
        actionOn->setText(QApplication::translate("ICASPHPlusClass", "On", Q_NULLPTR));
        actionOff->setText(QApplication::translate("ICASPHPlusClass", "Off", Q_NULLPTR));
        SPHTypeName->setText(QApplication::translate("ICASPHPlusClass", "IISPH DamBreak", Q_NULLPTR));
        menuSPH->setTitle(QApplication::translate("ICASPHPlusClass", "SPH Type", Q_NULLPTR));
        menuSPH_Interaction->setTitle(QApplication::translate("ICASPHPlusClass", "SPH Interaction", Q_NULLPTR));
        menuDistanceField->setTitle(QApplication::translate("ICASPHPlusClass", "DistanceField", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ICASPHPlusClass: public Ui_ICASPHPlusClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ICASPHPLUS_H
