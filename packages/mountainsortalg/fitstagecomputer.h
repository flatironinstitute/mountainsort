#ifndef FITSTAGECOMPUTER_H
#define FITSTAGECOMPUTER_H

#include "mda32.h"

class FitStageComputerPrivate;
class FitStageComputer {
public:
    friend class FitStageComputerPrivate;
    FitStageComputer();
    virtual ~FitStageComputer();

    void setTemplates(const Mda32& templates);
    void setTimesLabels(const QVector<double>& times, const QVector<int>& labels);
    void processTimeChunk(bigint t, const Mda32& X, bigint padding_left, bigint padding_right);

    void finalize();
    QVector<bigint> eventIndicesToUse() const; //the output

private:
    FitStageComputerPrivate* d;
};

#endif // FitStageCOMPUTER_H
