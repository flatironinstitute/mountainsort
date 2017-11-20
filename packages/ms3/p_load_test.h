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
#ifndef P_LOAD_TEST_H
#define P_LOAD_TEST_H

#include "mlutil.h"

struct P_load_test_opts {
    bigint num_cpu_ops = 0;
    bigint num_read_bytes = 0;
    bigint num_write_bytes = 0;
};

bool p_load_test(QString stats_out, P_load_test_opts opts);

bool p_misc_test(QString dir, QString info_out);

#endif // P_LOAD_TEST_H
