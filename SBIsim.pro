QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

QMAKE_LFLAGS += -Wl,-rpath,/home/somayeh/xraylib-master/install/lib
LIBS += -L/home/somayeh/xraylib-master/install/lib -lxrl
# ITK_LIBS = $$system(find /usr/lib -name 'libITK*-5.2.so' -exec basename {} .so \; | sed 's/^/ -l/' | tr '\n' ' ')
# LIBS += -L/usr/lib $$ITK_LIBS
LIBS += -L/usr/lib -lITKCommon-5.2 -lITKIOImageBase-5.2 -lITKImageFeature-5.2 -lITKIOPNG-5.2 -lITKIOJPEG-5.2 \
        -lITKStatistics-5.2 -lITKTransform-5.2 -lITKLabelMap-5.2 -litkvnl-5.2 -litkvnl_algo-5.2 #lITKTransformFactory-5.2
INCLUDEPATH +=/usr/include/ITK-5.2

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    Diffuser.cpp \
    DiffuserAndObject.cpp \
    Object.cpp \
    Setup.cpp \
    Simulation.cpp \
    SourceAndDetector.cpp \
    main.cpp \
    MainSBI.cpp

HEADERS += \
    Diffuser.h \
    DiffuserAndObject.h \
    ImageExporter.h \
    ImageExporter.tpp \
    MainSBI.h \
    Materials.h \
    Object.h \
    ParseChemicalFormula.h \
    PhotonAttenuation.h \
    Setup.h \
    Simulation.h \
    SourceAndDetector.h

FORMS += \
    Diffuser.ui \
    DiffuserAndObject.ui \
    MainSBI.ui \
    Setup.ui \
    Simulation.ui \
    SourceAndDetector.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
