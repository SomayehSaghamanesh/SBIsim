QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

QMAKE_LFLAGS += -Wl,-rpath,/home/somayeh/xraylib-master/install/lib
LIBS += -L/home/somayeh/xraylib-master/install/lib -lxrl

# ITK_LIBS = $$system(find /usr/lib -name 'libITK*-5.2.so' -exec basename {} .so \; | sed 's/^/ -l/' | tr '\n' ' ')
# LIBS += -L/usr/lib $$ITK_LIBS
LIBS += -L/usr/local/lib -lITKCommon-6.0 -lITKIOImageBase-6.0 -lITKImageFeature-6.0 -lITKIONIFTI-6.0 -lITKIOTIFF-6.0\
        -lITKStatistics-6.0 -lITKTransform-6.0 -lITKLabelMap-6.0 -litkvnl-6.0 -litkvnl_algo-6.0 -litksys-6.0
INCLUDEPATH +=/usr/local/include/ITK-6.0

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
    DiffuserAndObject.ui \
    MainSBI.ui \
    Setup.ui \
    Simulation.ui \
    SourceAndDetector.ui

# Default rules for deployment.
# qnx: target.path = /tmp/$${TARGET}/bin
# else: unix:!android: target.path = /opt/$${TARGET}/bin
# !isEmpty(target.path): INSTALLS += target

# unix:!macx: LIBS += -L$$PWD/../../../../usr/local/lib/ -lITKCommon-6.0 -lITKIOImageBase-6.0 -lITKImageFeature-6.0 -lITKIONIFTI-6.0 \
# -lITKStatistics-6.0 -lITKTransform-6.0 -lITKLabelMap-6.0 -litkvnl-6.0 -litkvnl_algo-6.0

# INCLUDEPATH += $$PWD/../../../../usr/local/include/ITK-6.0
# DEPENDPATH += $$PWD/../../../../usr/local/include/ITK-6.0

DISTFILES += \
    build/Desktop_Qt_6_7_0-Debug/SBIsim_icon.tif \
    img/SBIsim_icon.tif

RESOURCES += \
    Resources.qrc
