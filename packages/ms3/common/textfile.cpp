/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "textfile.h"
#include "mlutil.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>

QString TextFile::read(const QString& fname, QTextCodec* codec)
{
    QFile file(fname);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return QString();
    }
    QTextStream ts(&file);
    if (codec != 0)
        ts.setCodec(codec);
    QString ret = ts.readAll();
    file.close();
    return ret;
}

bool TextFile::write_single_try(const QString& fname, const QString& txt, QTextCodec* codec)
{
    /*
     * Modification on 5/23/16 by jfm
     * We don't want an program to try to read this while we have only partially completed writing the file.
     * Therefore we now create a temporary file and then copy it over
     */

    QString tmp_fname = fname + ".tf." + MLUtil::makeRandomId(6) + ".tmp";

    //if a file with this name already exists, we need to remove it
    //(should we really do this before testing whether writing is successful? I think yes)
    if (QFile::exists(fname)) {
        if (!QFile::remove(fname)) {
            if (QFile::exists(fname)) {
                qWarning() << "Problem in TextFile::write. Could not remove file even though it exists" << fname;
                return false;
            }
        }
    }

    //write text to temporary file
    QFile file(tmp_fname);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qWarning() << "Problem in TextFile::write. Could not open for writing... " << tmp_fname;
        return false;
    }
    QTextStream ts(&file);
    if (codec != 0) {
        ts.setAutoDetectUnicode(false);
        ts.setCodec(codec);
    }
    ts << txt;
    ts.flush();
    file.close();

    //check the contents of the file (is this overkill?)
    QString txt_test = TextFile::read(tmp_fname, codec);
    if (txt_test != txt) {
        QFile::remove(tmp_fname);
        qWarning() << "Problem in TextFile::write. The contents of the file do not match what was expected." << tmp_fname << txt_test.count() << txt.count();
        return false;
    }

    //check again if we need to remove the file
    if (QFile::exists(fname)) {
        if (!QFile::remove(fname)) {
            if (QFile::exists(fname)) {
                qWarning() << "Problem in TextFile::write. Could not remove file even though it exists (**)" << fname;
                QFile::remove(tmp_fname);
                return false;
            }
        }
    }

    //finally, rename the file
    if (!QFile::rename(tmp_fname, fname)) {
        qWarning() << "Problem in TextFile::write. Unable to rename file at the end of the write command" << fname << "Src/dst files exist?:" << QFile::exists(tmp_fname) << QFile::exists(fname);
        QFile::remove(tmp_fname);
        return false;
    }

    return true;
}

bool TextFile::write(const QString& fname, const QString& txt, QTextCodec* codec)
{
    int num_tries = 2;
    /*
    QFileInfo finfo(fname);
    if (!finfo.isWritable()) {
        qWarning() << "Problem in TextFile::write. File is not writable" << fname;
    }
    */
    for (int i = 0; i < num_tries; i++) {
        if (TextFile::write_single_try(fname, txt, codec)) {
            return true;
        }
    }
    qWarning() << "Problem in TextFile -- unable to write file: " << fname;
    return false;
}
