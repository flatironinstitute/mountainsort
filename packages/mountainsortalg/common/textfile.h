#ifndef TEXTFILE_H
#define TEXTFILE_H

#include <QString>

namespace TextFile {
QString read(const QString& fname, QTextCodec* codec = 0);
bool write(const QString& fname, const QString& txt, QTextCodec* codec = 0);
bool write_single_try(const QString& fname, const QString& txt, QTextCodec* codec = 0);
};

#endif // TEXTFILE_H

