#ifndef CLPARAMS_H
#define CLPARAMS_H

#include <QMap>
#include <QVariant>

class CLParams {
public:
    CLParams(int argc, char* argv[]);
    QMap<QString, QVariant> named_parameters;
    QList<QString> unnamed_parameters;
    bool success;
    QString error_message;
};

#endif // CLPARAMS_H

