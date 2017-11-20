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
#include "p_create_firings.h"

#include "mda.h"

#include <mda32.h>

bool p_create_firings(QString event_times, QString labels, QString amplitudes, QString firings_out, int central_channel)
{
    Mda ET(event_times);
    Mda LA(labels);
    Mda32 AM;

    bigint L = ET.totalSize();

    int R = 3;
    if (!amplitudes.isEmpty()) {
        R = 4;
        AM.read(amplitudes);
    }
    Mda firings(R, L);
    for (bigint i = 0; i < L; i++) {
        firings.setValue(central_channel, 0, i);
        firings.setValue(ET.value(i), 1, i);
        firings.setValue(LA.value(i), 2, i);
        if (R >= 4)
            firings.setValue(AM.value(i), 3, i);
    }

    return firings.write64(firings_out);
}
