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

#include "p_apply_timestamp_offset.h"
#include <mda.h>

bool p_apply_timestamp_offset(QString firings_path, QString firings_out_path, double timestamp_offset)
{
    Mda firings(firings_path);
    if ((firings.N1() == 1) || (firings.N2() == 1)) {
        //event times
        for (bigint i = 0; i < firings.totalSize(); i++) {
            firings.set(firings.get(i) + timestamp_offset, i);
        }
    }
    else {
        for (bigint i = 0; i < firings.N2(); i++) {
            firings.setValue(firings.value(1, i) + timestamp_offset, 1, i);
        }
    }
    return firings.write64(firings_out_path);
}
