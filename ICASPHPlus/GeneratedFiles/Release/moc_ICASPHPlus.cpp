/****************************************************************************
** Meta object code from reading C++ file 'ICASPHPlus.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../ICASPHPlus.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ICASPHPlus.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ICASPHPlus_t {
    QByteArrayData data[14];
    char stringdata0[194];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ICASPHPlus_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ICASPHPlus_t qt_meta_stringdata_ICASPHPlus = {
    {
QT_MOC_LITERAL(0, 0, 10), // "ICASPHPlus"
QT_MOC_LITERAL(1, 11, 7), // "onStart"
QT_MOC_LITERAL(2, 19, 0), // ""
QT_MOC_LITERAL(3, 20, 6), // "onStop"
QT_MOC_LITERAL(4, 27, 9), // "isSavePic"
QT_MOC_LITERAL(5, 37, 19), // "onDrawDistanceField"
QT_MOC_LITERAL(6, 57, 20), // "offDrawDistanceField"
QT_MOC_LITERAL(7, 78, 11), // "changeIISPH"
QT_MOC_LITERAL(8, 90, 17), // "changeSingleBreak"
QT_MOC_LITERAL(9, 108, 17), // "changeDoubleBreak"
QT_MOC_LITERAL(10, 126, 17), // "changeStaticRigid"
QT_MOC_LITERAL(11, 144, 13), // "changeSurface"
QT_MOC_LITERAL(12, 158, 11), // "changeWCSPH"
QT_MOC_LITERAL(13, 170, 23) // "isShowBoundaryPartilces"

    },
    "ICASPHPlus\0onStart\0\0onStop\0isSavePic\0"
    "onDrawDistanceField\0offDrawDistanceField\0"
    "changeIISPH\0changeSingleBreak\0"
    "changeDoubleBreak\0changeStaticRigid\0"
    "changeSurface\0changeWCSPH\0"
    "isShowBoundaryPartilces"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ICASPHPlus[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   74,    2, 0x08 /* Private */,
       3,    0,   75,    2, 0x08 /* Private */,
       4,    0,   76,    2, 0x08 /* Private */,
       5,    0,   77,    2, 0x08 /* Private */,
       6,    0,   78,    2, 0x08 /* Private */,
       7,    0,   79,    2, 0x08 /* Private */,
       8,    0,   80,    2, 0x08 /* Private */,
       9,    0,   81,    2, 0x08 /* Private */,
      10,    0,   82,    2, 0x08 /* Private */,
      11,    0,   83,    2, 0x08 /* Private */,
      12,    0,   84,    2, 0x08 /* Private */,
      13,    0,   85,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void ICASPHPlus::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ICASPHPlus *_t = static_cast<ICASPHPlus *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->onStart(); break;
        case 1: _t->onStop(); break;
        case 2: _t->isSavePic(); break;
        case 3: _t->onDrawDistanceField(); break;
        case 4: _t->offDrawDistanceField(); break;
        case 5: _t->changeIISPH(); break;
        case 6: _t->changeSingleBreak(); break;
        case 7: _t->changeDoubleBreak(); break;
        case 8: _t->changeStaticRigid(); break;
        case 9: _t->changeSurface(); break;
        case 10: _t->changeWCSPH(); break;
        case 11: _t->isShowBoundaryPartilces(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject ICASPHPlus::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_ICASPHPlus.data,
      qt_meta_data_ICASPHPlus,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ICASPHPlus::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ICASPHPlus::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ICASPHPlus.stringdata0))
        return static_cast<void*>(const_cast< ICASPHPlus*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int ICASPHPlus::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE