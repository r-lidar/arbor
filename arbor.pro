TEMPLATE = lib
CONFIG += shared

QMAKE_CXXFLAGS += -std=c++20
QMAKE_CXXFLAGS += -O2 # compile O2 even in debug mode

DESTDIR     = $$OUT_PWD/libarbor
OBJECTS_DIR = $$OUT_PWD/libarbor/.obj
MOC_DIR     = $$OUT_PWD/libarbor/.moc

# OpenMP Support
unix {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS   += -fopenmp
}
win32-msvc* {
    QMAKE_CXXFLAGS += /openmp
}
win32-g++ {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS   += -fopenmp
}

DEFINES += NOGDAL

SOURCES += \
    $$files(src/api/*.cpp) \
    $$files(src/pointcloud/*.cpp) \
    $$files(src/fitting/*.cpp) \
    $$files(src/qsm/*.cpp) \
    $$files(src/dtm/*.cpp) \
    $$files(src/qsf/*.cpp) \
    $$files(src/seed/*.cpp) \
    $$files(src/segment/*.cpp) \
    $$files(src/utils/*.cpp) \
    $$files(src/vendor/dbscan/*.cpp) \
    $$files(src/vendor/hporro/*.cpp) \
    $$files(src/vendor/libqsm/*.cpp) \
    $$files(src/vendor/ptd/*.cpp)

INCLUDEPATH += \
    src/api \
    src/dtm \
    src/fitting \
    src/pointcloud \
    src/nanoflann \
    src/qsf \
    src/qsm \
    src/seed \
    src/segment \
    src/utils \
    src/vendor/dbscan \
    src/vendor \
    src/vendor/ptd \
    src/vendor/libqsm \
    src
