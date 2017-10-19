#include "clparams.h"

QVariant clp_string_to_variant(const QString& str);

CLParams::CLParams(int argc, char* argv[])
{
    this->success = true; //let's be optimistic!

    //find the named and unnamed parameters checking for errors along the way
    for (int i = 1; i < argc; i++) {
        QString str = QString(argv[i]);
        if (str.startsWith("--")) {
            int ind2 = str.indexOf("=");
            QString name = str.mid(2, ind2 - 2);
            QString val = "";
            if (ind2 >= 0)
                val = str.mid(ind2 + 1);
            if (name.isEmpty()) {
                this->success = false;
                this->error_message = "Problem with parameter: " + str;
                return;
            }
            QVariant val2 = clp_string_to_variant(val);
            if (this->named_parameters.contains(name)) {
                QVariant tmp = this->named_parameters[name];
                QVariantList list;
                if (tmp.type() == QVariant::List) {
                    list = tmp.toList();
                }
                else {
                    list.append(tmp);
                }
                if (val2.type() == QVariant::List)
                    list.append(val2.toList());
                else
                    list.append(val2);
                this->named_parameters[name] = list;
            }
            else {
                this->named_parameters[name] = val2;
            }
        }
        else {
            this->unnamed_parameters << str;
        }
    }
}

bool clp_is_long(const QString& str)
{
    bool ok;
    str.toLongLong(&ok);
    return ok;
}

bool clp_is_int(const QString& str)
{
    bool ok;
    str.toInt(&ok);
    return ok;
}

bool clp_is_float(const QString& str)
{
    bool ok;
    str.toFloat(&ok);
    return ok;
}

QVariant clp_string_to_variant(const QString& str)
{
    if (clp_is_long(str))
        return str.toLongLong();
    if (clp_is_int(str))
        return str.toInt();
    if (clp_is_float(str))
        return str.toFloat();
    if ((str.startsWith("[")) && (str.endsWith("]"))) {
        QString str2 = str.mid(1, str.count() - 2);
        QStringList list = str2.split("][");
        if (list.count() == 1) {
            return list[0];
        }
        else {
            QVariantList ret;
            foreach (QString tmp, list)
                ret.append(tmp);
            return ret;
        }
    }
    return str;
}
