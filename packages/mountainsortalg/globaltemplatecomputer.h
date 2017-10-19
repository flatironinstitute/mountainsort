#ifndef GLOBALTEMPLATECOMPUTER_H
#define GLOBALTEMPLATECOMPUTER_H

#include "mda32.h"

class GlobalTemplateComputerPrivate;
class GlobalTemplateComputer {
public:
    friend class GlobalTemplateComputerPrivate;
    GlobalTemplateComputer();
    virtual ~GlobalTemplateComputer();

    void setNumThreads(int num_threads);
    void setClipSize(int clip_size);
    void setTimesLabels(const QVector<double>& times, const QVector<int>& labels);
    void processTimeChunk(bigint t, const Mda32& X, bigint padding_left, bigint padding_right);

    Mda32 templates() const; //the output

private:
    GlobalTemplateComputerPrivate* d;
};

#endif // GLOBALTEMPLATECOMPUTER_H
