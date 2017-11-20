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
